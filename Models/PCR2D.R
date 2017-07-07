#' Fit a Scalar-on-Image Regression Model based on Principal Components
#' 
#' This function fits a scalar-on-image regression model based on a principal 
#' component analysis. The optimal number of components to use is found by 
#' k-fold cross-validation. Optionally, confidence intervals can be returned for
#' the coefficient image (pointwise) and for potential scalar covariates. 
#' 
#' Requirements: This function requires the packages \pkg{funData} and \pkg{MFPCA}.
#' 
#' @param scores The principal component scores for \code{N} images, organized 
#'   in an \code{N x K} matrix.
#' @param y A vector of length \code{N} containing the response values.
#' @param FixEf An optional matrix containing the intercept and additional 
#'   covariates in row-wise manner
#' @param functions the principal component functions as a \code{funData} object
#'   with \code{K} observations.
#' @param nGroups An integer, specifying the number of folds for the 
#'   cross-validation.
#' @param kVals A vector containing the maximal number of principal components 
#'   to test in the cross-validation (must all be smaller or equal to \code{K}.)
#' @param CI Logical. If \code{TRUE}, pointwise bootstrap confidence intervals 
#'   are calculated for the coefficient image and for the effects of the scalar 
#'   covariates. Defaults to \code{FALSE}.
#' @param CIit An integer giving the number of bootstrap iterations for the
#'   calculation of the confidence intervals. Defaults to \code{200}.
#'   
#' @return \item{kBest}{The optimal number of principal components to use found by cross-validation.} \item{alphaHat}{The estimated coefficients for the scalar 
#'   covariates.} \item{betaHat}{The estimated coefficient image.} \item{fitted.values}{The vector of fitted response values.}
#'   \item{CI}{The pointwise bootstrap confidence intervals for the coefficient image saved in an array containing the lower and upper bounds.}
#'   \item{CI_fixEf}{The bootstrap confidence intervals for the scalar coefficients saved in an array containing the lower and upper bounds.}
PCRCV <- function(scores, y, FixEf, functions, nGroups, kVals, CI = FALSE, CIit = 200)
{
  predTest <-rep(NA, nGroups)
  predTotal <- rep(NA, length(kVals))
  
  
  for(k in 1:length(kVals))
  {
    # sample group indices
    groupInd <- sample(rep(1:nGroups, length.out = length(y)))
    
    
    for(g in 1:nGroups)
    {
      trainInd <- which(groupInd != g)
      
      LM <- lm(y ~ FixEf + scores - 1, data = list(y = y[trainInd], FixEf = FixEf[trainInd, , drop = FALSE], scores = scores[trainInd, 1:kVals[k], drop = FALSE] )) # include intercept in fixed effects
      
      fitTest <- predict(LM, newdata = list(FixEf = FixEf[-trainInd, , drop = FALSE], scores = scores[-trainInd, 1:kVals[k], drop = FALSE]))
      
      predTest[g] <- mean((y[-trainInd] - fitTest)^2)
    }
    
    predTotal[k] <- mean(predTest)
  }
  
  # find optimal k and refit model
  kBest <- kVals[which.min(predTotal)]
  LM <- lm(y ~ FixEf+ scores[, 1:kBest] - 1) # include intercept in fixed effects
  
  if(CI)
  {
    bootEst <- array(NA, dim = c(CIit, nObsPoints(functions)))
    bootFixEf <- array(NA, dim = c(CIit, dim(FixEf)[2]))
    
    for(i in 1:CIit)
    {
      ind <- sample(length(y), replace = TRUE)
      LMboot <- lm(y[ind] ~ FixEf[ind, , drop = FALSE]+ scores[ind, 1:kBest, drop = FALSE] - 1)
      bootEst[i,,] <- MFPCA::expandBasisFunction(scores = matrix(LMboot$coef[-(1:dim(FixEf)[2])], nrow = 1), functions =extractObs(functions, 1:kBest))@X[1,,]
      bootFixEf[i,] <- LMboot$coef[1:dim(FixEf)[2]]
    }
    
    CI_beta <- apply(bootEst,2:3, quantile, c(0.025, 0.975))
    dimnames(CI_beta) <- list(c("lower", "upper"), NULL, NULL) 
    
    CI_fixEf <- apply(bootFixEf, 2,  quantile, c(0.025, 0.975))
    dimnames(CI_fixEf) <- list(c("lower", "upper"), NULL) 
  }
  else
  {
    CI_beta <- NULL
    CI_fixEf <- NULL
  }
  
  return(list(kBest = kBest,
              alphaHat = LM$coef[1:dim(FixEf)[2]],
              betaHat = MFPCA::expandBasisFunction(scores = matrix(LM$coef[-(1:dim(FixEf)[2])], nrow = 1), functions =extractObs(functions, 1:kBest)),
              fitted.values = LM$fitted.values,
              CI = CI_beta,
              CI_fixEf= CI_fixEf))
}