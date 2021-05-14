#' Mixed Logistic Support Vector Machine for Binary Data
#'
#' @param Y The outcome variable.
#' @param x A formula string contains the predictors.
#' @param random A string in lme4 format indicates the random effect model.
#' @param data The data set as a data frame.
#' @param initialRandomEffects The initial values for random effects.
#' @param ErrorTolerance The tolerance for log-likelihood.
#' @param ErrorTolerance0 The tolerance for eta (penalized quasi-likelihood, PQL).
#' @param MaxIterations The maximum iteration times for each run of PQL.
#' @param MaxIterations0 The maximum iteration times for PQL.
#' @param verbose The option to monitor each run of PQL or not.
#'
#' @return A list contains the random forest, mixed model, and random effects.
#' See the example below for the usage. A predict() function is also available below.
#' @examples
#'

# one hot encoding
one_hot <- function(data, columns){
  data[, columns] <- sapply(data[, columns], as.factor)
  nominal <- dummyVars("~.", data=data[,columns], levelsOnly = TRUE)
  oheed <- predict(nominal, newdata = data[,columns])
  # drop old
  data <- data[ , !(names(data) %in% columns)]
  # add new
  newdata <- cbind(data, new_col = oheed)
  return(newdata)
}

# Factorial type to numerical type
fac_to_num <- function(data, columns){
  num_ed <- sapply(data[, columns], as.numeric)
  num_ed <- num_ed-1
  # drop old
  data <- data[ , !(names(data) %in% columns)]
  # add new
  newdata <- cbind(data, num_ed)
  return(newdata)
}

Mix_SVM <- function(Y, features, random, data, initialRandomEffects=0,
                    kernel, C,epsilon,nu,gamma,type,
                    ErrorTolerance=0.001, MaxIterations=200,
                    ErrorTolerance0=0.001, MaxIterations0=15, verbose=FALSE) {
  
  # Condition that indicates the loop has not converged or run out of iterations
  ContinueCondition0 <- TRUE
  iterations0 = 0
  
  # Get initial values
  
  mu = rep(mean(Y),length(Y))
  eta = log(mu/(1-mu))
  y = eta + (Y-mu)/(mu*(1-mu))
  # weights = mu*(1-mu)
  
  AdjustedTarget <- y - initialRandomEffects
  
  f0 = as.formula(paste0('resi ~ -1 + ', random))
  
  oldLogLik = oldEta = -Inf
  newy = y
  
  # PQL
  while(ContinueCondition0) {
    
    iterations0 <- iterations0 + 1
    
    iterations = 0
    ContinueCondition <- TRUE
    
    # random forest + lmer
    # Expectation Maximization esque
    while(ContinueCondition) {
      
      iterations <- iterations + 1
      
      # random Forest
      #df = data.frame(x = data[names(data) %in% x], y = AdjustedTarget)
      # data$AdjustedTarget = AdjustedTarget
      # model = svm(AdjustedTarget ~ ., data=data, kernel = kernel, C=C, gamma=gamma, scale=FALSE,
      #           type = type, epsilon=epsilon, nu=nu)
      predictors = data[names(data) %in% features]
      mod = svm(y= AdjustedTarget, x= predictors, kernel = kernel, C=C, gamma=gamma, scale=TRUE,
                type = type, epsilon=epsilon, nu=nu, na.action = na.omit)
      
      # y - X*beta 
      pred = predict(mod, newdata=predictors)
      #resi = y - pred
      resi = newy - pred
      
      
      ## Estimate New Random Effects and Errors using lmer
      lmefit <- lmer(f0, data=data)
      
      # check convergence
      LogLik <- as.numeric(logLik(lmefit))
      
      ContinueCondition <- (abs(LogLik-oldLogLik)>ErrorTolerance & iterations < MaxIterations)
      oldLogLik <- LogLik
      
      # Extract (the only) random effects part (Zb) to make the new adjusted target
      AllEffects <- predict(lmefit)
      
      #  y-Zb
      #AdjustedTarget <- y - AllEffects
      AdjustedTarget <- newy - AllEffects
      
      # monitor the change the of logLikelihood
      if(verbose) print(c(iterations0,iterations,LogLik))
    }
    
    neweta  <-pred + AllEffects
    newmu <- 1/(1+exp(-neweta))
    newY <- (newmu>=0.35)
    newy <- eta + (newY-newmu)/(newmu*(1-newmu))
    AdjustedTarget <- newy - AllEffects
    
    print(c(iter = iterations0, maxEtaChange=max(abs(neweta-oldEta))))
    
    ContinueCondition0 <- (max(abs(neweta-oldEta))>ErrorTolerance0 & iterations0 < MaxIterations0)
    oldEta <- eta
  }
  
  result <- list(SVR=mod, MixedModel=lmefit, RandomEffects=ranef(lmefit),
                 IterationsUsed=iterations0)
  
  return(result)
}


# predict 
predict.Mix_SVM <- function(object, newdata, EstimateRE=TRUE){
  
  svmPrediction <- predict(object$SVR,newdata=newdata)
  
  # If not estimate random effects, just use the forest for prediction.
  if(!EstimateRE){
    return(svmPrediction)
  }
  
  RandomEffects <- predict(object$MixedModel, newdata=newdata, allow.new.levels=T)
  
  completePrediction = svmPrediction + RandomEffects
  
  return(completePrediction)
}

# Over all prediction
# mean_RE=TRUE will get the average Random Effects in training samples
# add in the final prediction
# FALSE is default due to no literature recommendations
predicting <- function(model, data, mean_RE=FALSE){
  y_prob = predict.Mix_SVM(model, data, EstimateRE=FALSE)
  if(!mean_RE){
    return(y_prob)
  }
  random_effects = model$RandomEffects$empid
  random_effects = mean(random_effects$`(Intercept)`)
  y_prob_mRE = y_prob+random_effects
  return(y_prob_mRE)
}
