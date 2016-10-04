######
# Helper functions


getVars <- function(model.obj) {
  xvars <- AlteryxPredictive::getXVars2(model.obj)
  yvar <- as.character(formula(model.obj))[2]
  list(yvar = yvar, xvars = xvars)
}

# Remove variables from a list of factor levels that have a value of zero (e.g.,
# numeric and integer predictors)
noZeroLevels <- function(ll) {
  for(i in names(ll)) {
    if(is.numeric(ll[[i]]))
      ll[[i]] <- NULL
  }
  ll
}

# Remove variables from a list of factor levels that have a NULL value
noNullLevels <- function(ll) {
  for(i in names(ll)) {
    if(is.null(ll[[i]]))
      ll[[i]] <- NULL
  }
  ll
}

# Obtain the levels for each factor variable used in a model
getXlevels <- function(model.obj) {
  model.class <- class(model.obj)[1]
  if (model.class %in% c("glm", "svyglm", "rxLogit", "rxDTree", "rxDForest", "earth", "naiveBayes", "svm.formula", "nnet.formula")) {
    xlevels <- model.obj$xlevels
  } else {
    if (model.class == "rpart") {
      xlevels <- attr(model.obj, "xlevels")
    } else {
      if (model.class == "randomForest.formula") {
        xlevels <- noZeroLevels(model.obj$forest$xlevels)
      } else {
        if (model.class == "gbm") {
          xlevels <- model.obj$var.levels
          names(xlevels) <- getVars(model.obj)$xvars
          xlevels <- noZeroLevels(xlevels)
        } else {
          stop.Alteryx("There is a model of an unexpected class included. Please remove it and re-run the module.")
        }
      }
    }
  }
  xlevels
}

# Obtain the levels of a model's target variable
getYlevels <- function(model.obj) {
  model.class <- class(model.obj)[1]
  if (model.class %in% c("glm", "svyglm")) {
    #    y.name <- unlist(strsplit(as.character(model.obj$call)[2], " ~ "))[1]
    y.name <- as.character(model.obj$formula)[2]
    y.var <- eval(parse(text = paste("model.obj$data$", y.name, sep = "")))
    ylevels <- levels(y.var)
  } else {
    if (model.class == "rpart") {
      ylevels <- attributes(model.obj)$ylevels
    } else {
      if (model.class == "randomForest.formula") {
        ylevels <- attr(model.obj$confusion, "dimnames")[[1]]
      } else {
        if (model.class %in% c("rxLogit", "rxDTree", "rxDForest")) {
          ylevels <- model.obj$yinfo$levels
        } else {
          if (model.class == "earth" || model.class == "naiveBayes" || model.class == "svm.formula" ) {
            ylevels <- model.obj$levels
          } else {
            if (model.class == "nnet.formula") {
              ylevels <- model.obj$lev
            } else { # gbm
              ylevels <- model.obj$target.levels
            }
          }
        }
      }
    }
  }
  ylevels
}

# matchLevels coerces the levels in new data factors to exactly match the levels
# of factors in the original data, and is needed for Revo ScaleR models
matchLevels <- function(nd, ol) {
  # Address the non-standard way randomForest returns xlevel values
  check.ol <- sapply(ol, is.numeric)
  ol[check.ol] <- NULL
  # if the model does appear to have factors then sort out the levels
  if (!(is.null(ol) || length(ol) == 0)) {
    factor.test <- sapply(nd, is.factor)
    the.factors <- names(nd)[factor.test]
    if (!all(names(ol) %in% the.factors))
      stop("There are factor variables in the model that are not present in the data to be scored.")
    
    # The function to use with sapply to determine which factors have different
    # levels in the new data versus the levels used in model estimation.
    checkLevels <- function(z) {
      current.levels <- levels(nd[[z]])
      desired.levels <- ol[[z]]
      !all(current.levels == desired.levels)
    }
    the.factors <- the.factors[sapply(the.factors, checkLevels)]
    
    # The function to use with lapply to relevel a factor
    relevelFac <- function(z) {
      orig.factor <- nd[[z]]
      these.levels <- ol[[z]]
      new.factor <- factor(orig.factor, levels = these.levels)
      new.factor
    }
    
    # Relevel the factors that differ from their levels in model estimation
    new.factor.list <- lapply(the.factors, relevelFac)
    names(new.factor.list) <- the.factors
    nd[the.factors] <- new.factor.list
  }
  
  nd
}

# modelProbs gives the predictecd probability of the desired target level
modelProbs <- function(model.object, y.levels, targ.level, new.data) {
  the.type <- class(model.object)[1]
  if(the.type == "glm") {
    pred.prob <- predict(model.object, newdata = new.data, type = "response")
  } else {
    if (the.type == "randomForest.formula") {
      pred.prob  <- predict(model.object, newdata = matchLevels(new.data, model.object$forest$xlevels), type="prob")[ , targ.level]
    } else {
      if (the.type == "gbm") {
        pred <- predict(model.object, newdata = new.data, n.trees = model.object$best.trees, type = "response")
        if (model.object$distribution == "multinomial") {
          pred.prob <- pred[,1,1]
          if (model.object$target.levels[2] == targ.level) {
            pred.prob <- pred[,2,1]
          }
        } else {
          if (y.levels[1] == targ.level) {
            pred.prob <- 1 - pred
          } else {
            pred.prob <- pred
          }
        }
      } else {
        if (the.type == "rpart") {
          pred.prob <- predict(model.object, newdata = new.data)[ , targ.level]
        } else {
          if (the.type == "rxLogit") {
            if (model.object$yinfo$level[1] == targ.level) {
              pred.prob <- 1 - rxPredict(model.object, data = matchLevels(new.data, model.object$xlevels), type = "response")
            } else {
              pred.prob <- rxPredict(model.object, data = matchLevels(new.data, model.object$xlevels), type = "response")
            }
          } else {
            if (the.type == "rxDTree" || the.type == "rxDForest") {
              pred.probs <- rxPredict(model.object, data = matchLevels(new.data, model.object$xlevels), type = "prob")
              if (model.object$yinfo$level[1] == targ.level) {
                pred.prob <- pred.probs[[1]]
              } else {
                pred.prob <- pred.probs[[2]]
              }
            } else {
              if (the.type == "naiveBayes" ) {
                pred.prob <- predict(model.object, newdata = new.data, type="raw")[ , targ.level]
              } else {
                if (the.type == "svm.formula") {
                  pred.prob <- attr(predict(model.object, newdata = new.data, decision.values = TRUE, probability = TRUE),"probabilities")[ , targ.level]
                } else {
                  if (the.type == "nnet.formula") {
                    pred.prob <- predict(model.object, newdata = new.data, type = "raw")
                  } else { # earth
                    if (the.type != "earth")
                      stop.Alteryx("One of the provided models is of an unknown type")
                    if (model.object$glm.list[[1]]$family$family != "binomial")
                      stop.Alteryx("A Spline Model that is not a binary classifier is included in the set of models, only binary classifiers are appropriate for constructing lift charts.")
                    if (model.object$levels[1] == targ.level) {
                      pred.prob <- 1 - as.numeric(predict(model.object, newdata = new.data, type = "response"))
                    } else {
                      pred.prob <- as.numeric(predict(model.object, newdata = new.data, type = "response"))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  pred.prob
}

# The function giniCoef calculates the Gini coefficient of the distribution
# of ordered values. The order will typcially be based on a predicted score.
# Optionally, the function can produce a Lorenz curve plot of the data as a
# side effect. 
giniCoef <- function (y, wts = rep(1, length(y)), order.var = NULL, target.level = NULL, decending = TRUE, Lorenz = FALSE){
  if (is.factor(y)) {
    if (!is.null(target.level)) {
      y <- as.numeric(as.character(y) == target.level)
    } else {
      y <- as.numeric(as.character(y) == levels(y)[2])
    }
  }
  if (is.null(order.var))
    order.var <- y
  y <- y[order(order.var, decreasing = decending)]
  # Normalize the weights to sum to one
  wts <- wts[order(order.var, decreasing = decending)]/sum(wts)
  n <- length(y) + 1
  if (decending) {
    cum.wts <- c(0, cumsum(wts))
    cum.wt.y <- c(0, cumsum(wts*y)/sum(wts*y))
    gini <- sum((cum.wts[-1] - cum.wts[-n])*(cum.wt.y[-1] + cum.wt.y[-n])) - 1
  } else {
    cum.wts <- c(0, cumsum(wts))
    cum.wt.y <- c(0, cumsum(wts*y)/sum(wts*y))
    gini <- 1 - sum((cum.wts[-1] - cum.wts[-n])*(cum.wt.y[-1] + cum.wt.y[-n]))
  }
  # The optional plot of the Lorenz curve
  if (Lorenz) {
    plot(cum.wts, cum.wt.y, type = "l", lwd = 2, col = "forestgreen")
    lines(c(0,1), c(0,1), lty = 2, lwd = 2, col = "gray")
  }
  gini
}

########
# The main part of the macro

nRecordLimit <- %Question.num.records%
  
  # Get the data and models sorted out by use of the meta data, and then get the
  # first chunk of data and determine the path of an XDF file if the input is XDF
  is.XDF <- FALSE
for (i in 1:2) {
  meta.data <- read.AlteryxMetaInfo(name = paste("#", i, sep = ""))
  names.field <- as.character(meta.data$Name)
  if(names.field[1] == "Name" & names.field[2] == "Object" & length(names.field) == 2) {
    model.stream <- paste("#", i, sep = "")
    the.models <- read.Alteryx(model.stream)
  } else {
    data.stream <- paste("#", i, sep = "")
    the.source <- as.character(meta.data$Source)
    if (all(substr(the.source, 3, 9) == "Context")) {
      library(rjson)
      meta.list <- fromJSON(the.source[1])
      if (meta.list$Context == "XDF") {
        is.XDF <- TRUE
        xdf.path <- meta.list$File.Loc
      } else {
        stop.Alteryx("At this time only scoring with XDF files is supported")
      }
    }
    the.data <- read.Alteryx.First(data.stream, nRecordLimit = nRecordLimit)
  }
}

# The names in the.data
data.names <- names(the.data)

# Get the user provided parameters
cumulative <- '%Question.cumulative%'
incremental <- '%Question.incremental%'
true.resp <- as.numeric(trim.blanks('%Question.true_resp%'))
targ.level <- trim.blanks('%Question.targ_level%')
subtitle <- '%Question.samp_name%'

# Address data appearing in an Alteryx data stream first, then move on to
# compute contexts
if (!is.XDF) {
  # The allowed model types
  allowed.types <- c("glm", "svyglm", "randomForest.formula", "rpart", "gbm", "rxLogit", "rxDTree", "rxDForest", "earth", "naiveBayes", "svm.formula", "nnet.formula")
  
  # Put the models into a list and figure out what is in them and do error
  # checking
  yvar <- rep(NA, nrow(the.models))
  mod.class <- rep(NA, nrow(the.models))
  nmodels <- nrow(the.models)
  mod.list <- vector(mode = "list", length = nmodels)
  mod.names <- as.character(the.models$Name)
  the.models$Object <- as.character(the.models$Object)
  
  # Get the full set of variables actually used across all the models, which
  # will allow for a common set of records for comparsion purposes after
  # missing values have been addressed.
  all.xvars <- character(0)
  
  print("nmodels: ")
  print(nmodels)
  for (i in 1:nrow(the.models)) {
    mod.obj <- unserializeObject(the.models$Object[i]) 	
    mod.class[i] <- class(mod.obj)[1]
    print(paste("Model type:", mod.class[1]))
    if (!(mod.class[i] %in% allowed.types))
      stop.Alteryx(paste("Model", i, "is not one of the allowed types."))
    if (mod.class[i] %in% c("rxLogit", "rxGlm", "rxLinMod", "rxDTree", "rxDForest") && !("RevoScaleR" %in% row.names(installed.packages())))
      stop.Alteryx(paste("One of the provided models is a RevoScaleR model which is not supported on this machine."))
    if (!(mod.class[i] %in% c("naiveBayes", "svm.formula", "nnet.formula"))){
      the.vars <- getVars(mod.obj)
      if (!(the.vars$yvar %in% data.names) || !(the.vars$xvars %in% data.names))
        stop.Alteryx(paste("Not all of Model ", i, "'s variables are in the data.", sep=""))
      yvar[i] <- the.vars$yvar
      all.xvars <- c(all.xvars, the.vars$xvars)
      mod.list[[i]] <- mod.obj
    }
    if (mod.class[i] %in% c("naiveBayes", "svm.formula")) {
      yvar[i] <- mod.obj$yvars
      all.xvars <- c(all.xvars, mod.obj$xvars)
      mod.list[[i]] <- mod.obj
    }
    if (mod.class[i] == "nnet.formula") {
      these.vars <- names(attr(mod.obj$terms, "dataClasses"))
      yvar[i] <- these.vars[1]
      all.xvars <- c(all.xvars, these.vars[-1])
      mod.list[[i]] <- mod.obj
    }
  }
  
  saveRDS(mod.list, "C:\\Users\\kliu\\Documents\\Work\\201608_Logistic_Regression\\test\\modlist.RDS")
  # Do error checking on the target variables of the supplied models
  if (length(unique(yvar)) > 1)
    stop.Alteryx("Not all the models have the same target variable")
  y.var <- eval(parse(text=paste("the.data$", yvar[1], sep="")))
  if (!is.factor(y.var))
    stop.Alteryx("The target variable must be a factor")
  y.levels <- getYlevels(mod.list[[1]])
  print("y Levels is: ")
  print(y.levels)
  if (length(y.levels) != 2) 
    stop.Alteryx("Lift charts are only defined for binary target variables")
  if (!(targ.level %in% y.levels))
    stop("The target value given does not match any of the target variable values")
  
  # Load any needed libraries based on the classes of the input models
  if ("rpart" %in% mod.class)
    library(rpart)
  if ("randomForest.formula" %in% mod.class)
    library(randomForest)
  if ("gbm" %in% mod.class)
    library(gbm)
  if ("earth" %in% mod.class)
    library(earth)
  if ("naiveBayes" %in% mod.class)
    library(e1071)
  if ("svm.formula" %in% mod.class)
    library(e1071)
  if ("nnet.formula" %in% mod.class)
    library(nnet)
  
  # Determine which of the predictors is a factor
  unique.xvars <- unique(all.xvars)
  all.xfac <- unique.xvars[sapply(unique.xvars, function(x) is.factor(the.data[[x]]))]
  
  # Determine the set of levels for the predictor factor variables in each model
  all.mod.levels <- lapply(mod.list, getXlevels)
  names(all.mod.levels) <- mod.names
  
  # Determine the full set of factor levels across the factor predictors used
  # in all models. This is used to cull records 
  all.levels <- lapply(all.xfac, function(x) NULL)
  names(all.levels) <- all.xfac
  for (i in mod.names) {
    this.model.levels <- all.mod.levels[[i]]
    for (j in all.xfac) {
      all.levels[[j]] <- c(all.levels[[j]], this.model.levels[[j]])
    }
  }
  all.levels <- lapply(all.levels, unique)
  # Add in the levels of the target
  eval(parse(text = paste("all.levels$", yvar[1], " <- y.levels", sep = "")))
  
  # Determine the unique set of variables that appear across the models to be
  # examined, including the target variables
  all.vars <- unique(c(yvar[1], all.xvars))
  # Deal with possible new line characters that can happen with a lot of fields
  all.vars <- trim.blanks(gsub("\n", "", all.vars))
  
  ## PREPARE THE FIRST CHUNK OF DATA TO GET THE ORDERED BINARY VECTOR
  #
  # Remove any variables not needed for creating the lift charts
  rel.data <- the.data[, all.vars]
  # Remove records from the data with missing values
  rel.data <- na.omit(rel.data)
  
  # Delete the original data to free up memory
  rm(the.data)
  gc()
  
  # Remove records that have factors variables with levels not in the model
  #
  # The levels of the fields in the relevant data
  rd.levels <- noNullLevels(sapply(rel.data, function(x) if (class(x) == "factor") levels(x)))
  # Address the possibility that all factors have the same number of levels
  if (class(rd.levels) == "matrix") {
    tf <- attr(rd.levels, "dimnames")[[2]]
    rd.levels2 <- list()
    for (i in 1:length(tf))
      rd.levels2[[tf[i]]] <- rd.levels[,i]
    rd.levels <- rd.levels2
  }
  # Remove the records with unknown factor levels and re-do the factors to
  # remove the missing levels
  for (i in names(all.levels)) {
    if(!all(rd.levels[[i]] %in% all.levels[[i]])) {
      extra <- rd.levels[[i]][!(rd.levels[[i]] %in% all.levels[[i]])]
      if (length(extra) == 1)
        the.message <- paste("The category", extra, "in the field", i, "was not present in the data used to estimate the models. Records with this category value will be omitted.")
      else
        the.message <- paste("The categories", paste(extra, collapse = ", "), "in the field", i, "were not present in the data used to estimate the models. Records with these category values will be omitted.")
      AlteryxMessage(the.message, iType = 2, iPriority = 3)
      rel.data <- rel.data[!(as.character(rel.data[[i]]) %in% extra),]
      new.factor <- as.factor(as.character(rel.data[[i]]))
      rel.data[[i]] <- NULL
      rel.data[[i]] <- new.factor
      rm(new.factor)
    }
  }
  
  # Shuffle the data into random order. This addresses oddities with tree
  # models when the original data is structured with all responses of one type
  # coming first in the data set. This is only done for the first data chunk
  # since it is an issue expected to be encountered mostly with smaller data
  # sets, which will be read as a single chunk.
  set.seed(1)
  rel.data <- rel.data[order(runif(nrow(rel.data))), ]
  
  # Get the 0/1 version of the target values in this data chunk
  y.vals <- eval(parse(text=paste("rel.data$", yvar[1], sep="")))
  y.binary <- as.numeric(as.character(y.vals) == targ.level)
  
  # Intialize the predicted probabilities for the different models as a list
  model.probs <- list()
  for (i in 1:nmodels) {
    model.probs[[i]] <- modelProbs(mod.list[[i]], y.levels, targ.level, rel.data)
  }
  
  # PREPARE AND GET THE ORDERED 0/1 VALUES FOR THE DIFFERENT MODELS IN THE
  # REMAINING DATA CHUNKS
  #
  while (!is.null(the.data <- read.Alteryx.Next(data.stream))) {
    # Remove any variables not needed for creating the lift charts
    rel.data <- the.data[, all.vars]
    # Remove records from the data with missing values
    rel.data <- na.omit(rel.data)
    
    # Delete the original data to free up memory
    rm(the.data)
    gc()
    
    # Remove records that have factors variables with levels not in the model
    #
    # The levels of the fields in the relevant data
    rd.levels <- noNullLevels(sapply(rel.data, function(x) if (class(x) == "factor") levels(x)))
    # Address the possibility that all factors have the same number of levels
    if (class(rd.levels) == "matrix") {
      tf <- attr(rd.levels, "dimnames")[[2]]
      rd.levels2 <- list()
      for (i in 1:length(tf))
        rd.levels2[[tf[i]]] <- rd.levels[,i]
      rd.levels <- rd.levels2
    }
    # Remove the records with unknown factor levels and re-do the factors to
    # remove the missing levels
    for (i in names(all.levels)) {
      if(!all(rd.levels[[i]] %in% all.levels[[i]])) {
        extra <- rd.levels[[i]][!(rd.levels[[i]] %in% all.levels[[i]])]
        if (length(extra) == 1)
          the.message <- paste("The category", extra, "in the field", i, "was not present in the data used to estimate the models. Records with this category value will be omitted.")
        else
          the.message <- paste("The categories", paste(extra, collapse = ", "), "in the field", i, "were not present in the data used to estimate the models. Records with these category values will be omitted.")
        AlteryxMessage(the.message, iType = 2, iPriority = 3)
        rel.data <- rel.data[!(as.character(rel.data[[i]]) %in% extra),]
        new.factor <- as.factor(as.character(rel.data[[i]]))
        rel.data[[i]] <- NULL
        rel.data[[i]] <- new.factor
        rm(new.factor)
      }
    }
    
    # Get the 0/1 version of the target values in this data chunk and append it
    y.vals <- eval(parse(text=paste("rel.data$", yvar[1], sep="")))
    y.binary <- c(y.binary, as.numeric(as.character(y.vals) == targ.level))
    
    # Get the model probabilities for this chunk of data
    for (i in 1:nmodels) {
      model.probs[[i]] <- c(model.probs[[i]], modelProbs(mod.list[[i]], y.levels, targ.level, rel.data))
    }
  }
  names(model.probs) <- mod.names
  
  # Remove rel.data to minmize memory usage
  rm(rel.data)
  gc()
  
  # Information needed to construct the charts and tables. All the needed data
  # obtained first since the data could be chunked
  samp.resp <- sum(y.binary)/length(y.binary)
  samp.wt <- (samp.resp*(1 - true.resp))/(true.resp*(1 - samp.resp))
  
  # Construct the charts and the gains/response table
  the.tbl <- data.frame(Decile = as.character(seq(0, 100, 10)))
  
  # The function graphWHR takes information about the desired graph size (in
  # inches or centemeters), the desired graph resolution, and the base font
  # size used to determine the render size of the objects on the graph
  whr <- graphWHR(inches = '%Question.inches%', in.w = '%Question.in.w%', in.h = '%Question.in.h%', cm.w = '%Question.cm.w%', cm.h = '%Question.cm.h%', resolution = '%Question.graph.resolution%', print.high = FALSE)
  AlteryxGraph(2, width = whr[1], height = whr[2], res = whr[3], pointsize = '%Question.pointsize%')
  colr <- rep(palette(), ceiling(nmodels/length(palette())))
  
  # CREATE THE PLOTS, THE GAINS TABLE, THE AREA UNDER THE CURVE, AND THE GINI
  area.under.curve <- NULL
  gini <- NULL
  if (cumulative == "True") {
    iter.num <- 1
    plot(seq(0.1, 1, 0.1), seq(0.1, 1, 0.1), main="Weighted Cumulative Response Captured", sub=subtitle, xlab="Sample Proportion", ylab="Percent of Total Response Captured", type="l", lwd=2, xaxs="i", yaxs="i")
    for (i in mod.names) {
      var1 <- y.binary[order(model.probs[[i]], decreasing = TRUE)]
      var.ind1 <- rep(1, length(var1))
      var.ind1[var1 == 0] <- samp.wt
      var.ind <- cut(cumsum(var.ind1)/sum(var.ind1), seq(0, 1, 0.1), include.lowest = TRUE)
      var2 <- as.vector(by(var1, var.ind, sum))
      lines(seq(0.1, 1, 0.1), cumsum(var2)/sum(var2), col = colr[iter.num], lwd = 2)
      points(seq(0.1, 1, 0.1), cumsum(var2)/sum(var2), col = colr[iter.num], pch = iter.num)
      Type <- "Gains"
      cum.pct <- c(0,cumsum(var2)/sum(var2))
      area.under.curve <- c(area.under.curve, sum(0.05*(cum.pct[2:11] + cum.pct[1:10])))
      gini <- c(gini, giniCoef(y.binary, wts = var.ind1, order.var = model.probs[[i]]))
      gain <- c(0, 10*(cum.pct[2:11] - cum.pct[1:10]))
      eval(parse(text = paste("the.tbl$", i, "_Cum_Pct <- round(100*cum.pct, 3)", sep = "")))
      eval(parse(text = paste("the.tbl$", i, "_Lift <- round(gain, 3)", sep = "")))
      iter.num <- iter.num + 1
    }
    legend("bottomright", legend=mod.names, col=colr[1:nmodels], pch=1:nmodels, lty=1, lwd=2)
  } else {
    # Calculate the response rates. This needs to be done in a two-step
    # process since the maximum incremental response rate needs to be
    # determined in order to determine the extents of the y-axis in the
    # plots
    iter.num <- 1
    resp.matrix <- matrix(NA, nrow=10, ncol=nmodels)
    for (i in mod.names) {
      var1 <- y.binary[order(model.probs[[i]], decreasing = TRUE)]
      var.ind1 <- rep(1, length(var1))
      var.ind1[var1 == 0] <- samp.wt
      var.ind <- cut(cumsum(var.ind1)/sum(var.ind1), seq(0, 1, 0.1), include.lowest = TRUE)
      var2 <- as.vector(by(var1, var.ind, sum))
      var3 <- as.vector(by(var.ind1, var.ind, sum))
      resp.matrix[, iter.num] <- var2/var3
      iter.num <- iter.num + 1
    }
    # Create the base plot canvas
    max.resp <- max(resp.matrix)
    plot(seq(0.1, 1, 0.1), seq(0, max.resp, length=10), type = "n", main = "Weighted Incremental Response Rate", sub = subtitle, xlab = "Sample Percentile", ylab = "Resposne Rate", lwd = 2, xaxs = "i", yaxs = "i")
    lines(seq(0.1, 1, 0.1), rep(true.resp, 10), lwd=2)
    
    # The "grp" label of the table output
    Type <- "Resp Rate"
    
    # Add the incremental response rates to the plot for each model
    for (j in 1:nmodels) {
      lines(seq(0.1, 1, 0.1), as.vector(resp.matrix[, j]), col=colr[j], lwd=2)
      points(seq(0.1, 1, 0.1), as.vector(resp.matrix[, j]), col=colr[j], pch=j)
      eval(parse(text = paste("the.tbl$", mod.names[j], " <- c(0,as.vector(resp.matrix[, j]))", sep = "")))
    }
    legend("topright", legend=mod.names, col=colr[1:nmodels], pch=1:nmodels, lty=1, lwd=2)
  }
  
  # The case of data in an XDF file
} else {
  # The allowed model types
  allowed.types <- c("rxLogit", "rxDTree", "rxDForest")
  
  # Put the models into a list and figure out what is in them and do error
  # checking
  yvar <- rep(NA, nrow(the.models))
  mod.class <- rep(NA, nrow(the.models))
  nmodels <- nrow(the.models)
  mod.list <- vector(mode = "list", length = nmodels)
  mod.names <- as.character(the.models$Name)
  the.models$Object <- as.character(the.models$Object)
  
  # Get the full set of variables actually used across all the models, which
  # will allow for a common set of records for comparsion purposes after
  # missing values have been addressed.
  all.xvars <- character(0)
  for (i in 1:nrow(the.models)) {
    mod.obj <- unserializeObject(the.models$Object[i])
    mod.class[i] <- class(mod.obj)[1]
    if (!(mod.class[i] %in% allowed.types))
      stop.Alteryx(paste("Model", i, "is not one of the allowed types."))
    if (mod.class[i] %in% c("rxLogit", "rxGlm", "rxLinMod", "rxDTree", "rxDForest") && !("RevoScaleR" %in% row.names(installed.packages())))
      stop.Alteryx(paste("One of the provided models is a RevoScaleR model which is not supported on this machine."))
    the.vars <- getVars(mod.obj)
    if (!(the.vars$yvar %in% data.names) || !(the.vars$xvars %in% data.names))
      stop.Alteryx(paste("Not all of Model ", i, "'s variables are in the data.", sep=""))
    yvar[i] <- the.vars$yvar
    all.xvars <- c(all.xvars, the.vars$xvars)
    mod.list[[i]] <- mod.obj
  }
  
  # Do error checking on the target variables of the supplied models
  if (length(unique(yvar)) > 1)
    stop.Alteryx("Not all the models have the same target variable")
  y.var <- eval(parse(text=paste("the.data$", yvar[1], sep="")))
  if (!is.factor(y.var))
    stop.Alteryx("The target variable must be a factor")
  y.levels <- getYlevels(mod.list[[1]])
  if (length(y.levels) != 2) 
    stop.Alteryx("Lift charts are only defined for binary target variables")
  if (!(targ.level %in% y.levels))
    stop("The target value given does not match any of the target variable values")
  
  # Determine which of the predictors is a factor
  unique.xvars <- unique(all.xvars)
  all.xfac <- unique.xvars[sapply(unique.xvars, function(x) is.factor(the.data[[x]]))]
  
  # Determine the set of levels for the predictor factor variables in each model
  all.mod.levels <- lapply(mod.list, getXlevels)
  names(all.mod.levels) <- mod.names
  
  # Determine the full set of factor levels across the factor predictors used
  # in all models. This is used to cull records 
  all.levels <- lapply(all.xfac, function(x) NULL)
  names(all.levels) <- all.xfac
  for (i in mod.names) {
    this.model.levels <- all.mod.levels[[i]]
    for (j in all.xfac) {
      all.levels[[j]] <- c(all.levels[[j]], this.model.levels[[j]])
    }
  }
  all.levels <- lapply(all.levels, unique)
  # Add in the levels of the target
  eval(parse(text = paste("all.levels$", yvar[1], " <- y.levels", sep = "")))
  
  # Determine the unique set of variables that appear across the models to be
  # examined, including the target variables
  all.vars <- unique(c(yvar[1], all.xvars))
  
  ## PREPARE THE FIRST CHUNK OF DATA 	TO GET THE ORDERED BINARY VECTOR
  #
  # Setting up the read of the XDF file
  xdf.data <- RxXdfData(xdf.path, blocksPerRead = 1)
  rxOpen(xdf.data)
  this.df <- rxReadNext(xdf.data)
  # Remove any variables not needed for creating the lift charts
  rel.data <- this.df[, all.vars]
  # Remove records from the data with missing values
  rel.data <- na.omit(rel.data)
  
  # Delete this.df to free up memory
  rm(this.df)
  gc()
  
  # Shuffle the data into random order. This addresses oddities with tree
  # models when the original data is structured with all responses of one type
  # coming first in the data set. This is only done for the first data chunk
  # since it is an issue expected to be encountered mostly with smaller data
  # sets, which will be read as a single chunk.
  set.seed(1)
  rel.data <- rel.data[order(runif(nrow(rel.data))), ]
  
  # Remove records that have factors variables with levels not in the model
  #
  # The levels of the fields in the relevant data
  rd.levels <- noNullLevels(sapply(rel.data, function(x) if (class(x) == "factor") levels(x)))
  # Address the possibility that all factors have the same number of levels
  if (class(rd.levels) == "matrix") {
    tf <- attr(rd.levels, "dimnames")[[2]]
    rd.levels2 <- list()
    for (i in 1:length(tf))
      rd.levels2[[tf[i]]] <- rd.levels[,i]
    rd.levels <- rd.levels2
  }
  # Remove the records with unknown factor levels and re-do the factors to
  # remove the missing levels
  for (i in names(all.levels)) {
    if(!all(rd.levels[[i]] %in% all.levels[[i]])) {
      extra <- rd.levels[[i]][!(rd.levels[[i]] %in% all.levels[[i]])]
      if (length(extra) == 1)
        the.message <- paste("The category", extra, "in the field", i, "was not present in the data used to estimate the models. Records with this category value will be omitted.")
      else
        the.message <- paste("The categories", paste(extra, collapse = ", "), "in the field", i, "were not present in the data used to estimate the models. Records with these category values will be omitted.")
      AlteryxMessage(the.message, iType = 2, iPriority = 3)
      rel.data <- rel.data[!(as.character(rel.data[[i]]) %in% extra),]
      new.factor <- as.factor(as.character(rel.data[[i]]))
      rel.data[[i]] <- NULL
      rel.data[[i]] <- new.factor
      rm(new.factor)
    }
  }
  
  # Get the 0/1 version of the target values in this data chunk
  y.vals <- eval(parse(text=paste("rel.data$", yvar[1], sep="")))
  y.binary <- as.numeric(as.character(y.vals) == targ.level)
  
  # Intialize the predicted probabilities for the different models as a list
  model.probs <- list()
  for (i in 1:nmodels) {
    model.probs[[i]] <- modelProbs(mod.list[[i]], y.levels, targ.level, rel.data)
  }
  
  # PREPARE AND GET THE ORDERED 0/1 VALUES FOR THE DIFFERENT MODELS IN THE
  # REMAINING DATA CHUNKS
  #
  repeat {
    this.df <- rxReadNext(xdf.data)
    if (length(this.df) == 0) {
      break
    } else {
      # Remove any variables not needed for creating the lift charts
      rel.data <- this.df[, all.vars]
      # Remove records from the data with missing values
      rel.data <- na.omit(rel.data)
      
      # Delete this.df to free up memory
      rm(this.df)
      gc()
      
      # Remove records that have factors variables with levels not in the model
      #
      # The levels of the fields in the relevant data
      rd.levels <- noNullLevels(sapply(rel.data, function(x) if (class(x) == "factor") levels(x)))
      # Address the possibility that all factors have the same number of levels
      if (class(rd.levels) == "matrix") {
        tf <- attr(rd.levels, "dimnames")[[2]]
        rd.levels2 <- list()
        for (i in 1:length(tf))
          rd.levels2[[tf[i]]] <- rd.levels[,i]
        rd.levels <- rd.levels2
      }
      # Remove the records with unknown factor levels and re-do the factors to
      # remove the missing levels
      for (i in names(all.levels)) {
        if(!all(rd.levels[[i]] %in% all.levels[[i]])) {
          extra <- rd.levels[[i]][!(rd.levels[[i]] %in% all.levels[[i]])]
          if (length(extra) == 1)
            the.message <- paste("The category", extra, "in the field", i, "was not present in the data used to estimate the models. Records with this category value will be omitted.")
          else
            the.message <- paste("The categories", paste(extra, collapse = ", "), "in the field", i, "were not present in the data used to estimate the models. Records with these category values will be omitted.")
          AlteryxMessage(the.message, iType = 2, iPriority = 3)
          rel.data <- rel.data[!(as.character(rel.data[[i]]) %in% extra),]
          new.factor <- as.factor(as.character(rel.data[[i]]))
          rel.data[[i]] <- NULL
          rel.data[[i]] <- new.factor
          rm(new.factor)
        }
      }
      
      # Get the 0/1 version of the target values in this data chunk
      y.vals <- eval(parse(text=paste("rel.data$", yvar[1], sep="")))
      y.binary <- as.numeric(y.binary, as.character(y.vals) == targ.level)
      
      # Intialize the predicted probabilities for the different models as a list
      for (i in 1:nmodels) {
        model.probs[[i]] <- c(model.probs[[i]], modelProbs(mod.list[[i]], y.levels, targ.level, rel.data))
      }
    }
  }
  names(model.probs) <- mod.names
  
  # Remove rel.data to minmize memory usage
  rm(rel.data)
  gc()
  
  # Information needed to construct the charts and tables. All the needed data
  # obtained first since the data could be chunked
  samp.resp <- sum(y.binary)/length(y.binary)
  samp.wt <- (samp.resp*(1 - true.resp))/(true.resp*(1 - samp.resp))
  
  # Construct the charts and the gains/response table
  the.tbl <- data.frame(Decile = as.character(seq(0, 100, 10)))
  
  # The function graphWHR takes information about the desired graph size (in
  # inches or centemeters), the desired graph resolution, and the base font
  # size used to determine the render size of the objects on the graph
  whr <- graphWHR(inches = '%Question.inches%', in.w = '%Question.in.w%', in.h = '%Question.in.h%', cm.w = '%Question.cm.w%', cm.h = '%Question.cm.h%', resolution = '%Question.graph.resolution%', print.high = FALSE)
  AlteryxGraph(2, width = whr[1], height = whr[2], res = whr[3], pointsize = '%Question.pointsize%')
  colr <- rep(palette(), ceiling(nmodels/length(palette())))
  
  # CREATE THE PLOTS, THE GAINS TABLE, AND THE AREA UNDER THE CURVE
  area.under.curve <- NULL
  gini <- NULL
  if (cumulative == "True") {
    iter.num <- 1
    plot(seq(0.1, 1, 0.1), seq(0.1, 1, 0.1), main="Weighted Cumulative Response Captured", sub=subtitle, xlab="Sample Proportion", ylab="Percent of Total Response Captured", type="l", lwd=2, xaxs="i", yaxs="i")
    for (i in mod.names) {
      var1 <- y.binary[order(model.probs[[i]], decreasing = TRUE)]
      var.ind1 <- rep(1, length(var1))
      var.ind1[var1 == 0] <- samp.wt
      var.ind <- cut(cumsum(var.ind1)/sum(var.ind1), seq(0, 1, 0.1), include.lowest = TRUE)
      var2 <- as.vector(by(var1, var.ind, sum))
      lines(seq(0.1, 1, 0.1), cumsum(var2)/sum(var2), col = colr[iter.num], lwd = 2)
      points(seq(0.1, 1, 0.1), cumsum(var2)/sum(var2), col = colr[iter.num], pch = iter.num)
      Type <- "Gains"
      cum.pct <- c(0,cumsum(var2)/sum(var2))
      area.under.curve <- c(area.under.curve, sum(0.05*(cum.pct[2:11] + cum.pct[1:10])))
      gini <- c(gini, giniCoef(y.binary, wts = var.ind1, order.var = model.probs[[i]]))
      gain <- c(0, 10*(cum.pct[2:11] - cum.pct[1:10]))
      eval(parse(text = paste("the.tbl$", i, "_Cum_Pct <- round(100*cum.pct, 3)", sep = "")))
      eval(parse(text = paste("the.tbl$", i, "_Lift <- round(gain, 3)", sep = "")))
      iter.num <- iter.num + 1
    }
    legend("bottomright", legend=mod.names, col=colr[1:nmodels], pch=1:nmodels, lty=1, lwd=2)
  } else {
    # Calculate the response rates. This needs to be done in a two-step
    # process since the maximum incremental response rate needs to be
    # determined in order to determine the extents of the y-axis in the
    # plots
    iter.num <- 1
    resp.matrix <- matrix(NA, nrow=10, ncol=nmodels)
    for (i in mod.names) {
      var1 <- y.binary[order(model.probs[[i]], decreasing = TRUE)]
      var.ind1 <- rep(1, length(var1))
      var.ind1[var1 == 0] <- samp.wt
      var.ind <- cut(cumsum(var.ind1)/sum(var.ind1), seq(0, 1, 0.1), include.lowest = TRUE)
      var2 <- as.vector(by(var1, var.ind, sum))
      var3 <- as.vector(by(var.ind1, var.ind, sum))
      resp.matrix[, iter.num] <- var2/var3
      iter.num <- iter.num + 1
    }
    # Create the base plot canvas
    max.resp <- max(resp.matrix)
    plot(seq(0.1, 1, 0.1), seq(0, max.resp, length=10), type = "n", main = "Weighted Incremental Response Rate", sub = subtitle, xlab = "Sample Percentile", ylab = "Resposne Rate", lwd = 2, xaxs = "i", yaxs = "i")
    lines(seq(0.1, 1, 0.1), rep(true.resp, 10), lwd=2)
    
    # The "grp" label of the table output
    Type <- "Resp Rate"
    
    # Add the incremental response rates to the plot for each model
    for (j in 1:nmodels) {
      lines(seq(0.1, 1, 0.1), as.vector(resp.matrix[, j]), col=colr[j], lwd=2)
      points(seq(0.1, 1, 0.1), as.vector(resp.matrix[, j]), col=colr[j], pch=j)
      eval(parse(text = paste("the.tbl$", mod.names[j], " <- c(0,as.vector(resp.matrix[, j]))", sep = "")))
    }
    legend("topright", legend=mod.names, col=colr[1:nmodels], pch=1:nmodels, lty=1, lwd=2)
  }
}

# Close the AlteryxGraph device
invisible(dev.off())
# Write out the lift table
the.tbl2 <- wrapTable(format(the.tbl, trim = TRUE, digits = 3, nsmall = 1), width = 5)
the.tbl2 <- sub("  ", "Decile", the.tbl2)
grp.out <- data.frame(grp = c("type", rep("table", length(the.tbl2))), out = c(Type, the.tbl2))
if (cumulative == "True") {
  grp.out <- rbind(grp.out, data.frame(grp = rep("area", length(mod.names)), out = paste(mod.names, area.under.curve, gini, sep = "|")))
}
write.Alteryx(grp.out)
