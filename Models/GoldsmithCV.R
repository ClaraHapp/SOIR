
#' Algorithm of Goldsmith et al (2014) for Scalar-on-Image regression
#' @param hyper A list of hyperparameters with names sigEps, sigBeta, a, b
#' @param nGroups The number of groups to use for crossvalidation
#' @param X An array of image covariates
#' @param FixEf A matrix of fixed effects 
#' @param y A vector of response values
#' @param delta Matrix of neighbourhood structure, see function findNeigh
#' @param K The number of MCMC iterations for each run
#' @param nMCMC Thinning of MCMC chains, i.e. if nMCMC = 10, save only each 10th value
#' @param burnin The burnin length for each MCMC run
#' @param verbose Should the function produce an output to document the current state?
GoldsmithCV <- function(hyper, nGroups,
                        X, FixEf, y,
                        delta, 
                        K, nMCMC, burnin,
                        verbose = FALSE)
{
  # Define auxiliary parameters
  nPat <- dim(X)[2]
  nVox <- prod(dim(X)[1])
  nFixEf <- dim(FixEf)[2]
  FixEfInvRoot <- chol(solve(t(FixEf)%*%FixEf))
  
  allHyper <- expand.grid(hyper)
  allPred <- rep(NA, length(allHyper))
  
  for(i in 1:dim(allHyper)[1])
  {
    if(verbose)
      cat("Hyper: sigEps = ", allHyper$sigEps[i], ", sigBeta = ", allHyper$sigBeta[i], ", a = ", allHyper$a[i], ", b = ", allHyper$b[i], "\n")
    # groups for cross-validation
    groupInd <- sample(rep(1:nGroups, length.out = nPat))
    
    testPred <- rep(NA, nGroups)
    
    for(g in 1:nGroups)
    {
      trainInd <- which(groupInd != g)
      FixEfInvRoot <- chol(solve(t(FixEf[trainInd,])%*%FixEf[trainInd,]))
      
      fit <- .C("allGibbs",
                K= as.integer(K), nPat = as.integer(length(trainInd)), nVox = as.integer(nVox), nFixEf = as.integer(nFixEf), nDelta = as.integer(dim(delta)[2]), 
                nMCMC = as.integer(nMCMC), burnin = as.integer(burnin),
                X = as.double(X[,trainInd]), FixEf = as.double(FixEf[trainInd,]), Y = as.double(y[trainInd]),
                delta=as.integer(delta), FixEfInvRoot = as.double(FixEfInvRoot),
                sigEps = as.double(allHyper$sigEps[i]),
                sigBeta = as.double(allHyper$sigBeta[i]),
                a = as.double(allHyper$a[i]), b = as.double(allHyper$b[i]), 
                calcCI = as.integer(0), # no quantile CI calculation in CV
                alphaRes  = as.double(rep(0, nFixEf)), betaRes = as.double(rep(0, nVox)), gammaRes = as.double(rep(0, nVox)),
                betaLower = as.double(0), betaUpper = as.double(0),
                alphaLower = as.double(0), alphaUpper = as.double(0),
                NAOK = TRUE)
      fitTest <- FixEf[-trainInd,, drop = FALSE] %*% fit$alphaRes + apply(X[,-trainInd],2, function(x){sum(x * fit$betaRes)})
      
      testPred[g] <- sum((y[-trainInd] - fitTest)^2)
    } 
    
    allPred[i] <- mean(testPred)
  }
  
  
  bestHyper <- allHyper[which.min(allPred),]
  FixEfInvRoot <- chol(solve(t(FixEf)%*%FixEf))
  
  # refit model on full data using the best parameter combination and calculate CI
  fit <- .C("allGibbs",
            K= as.integer(K), nPat = as.integer(nPat), nVox = as.integer(nVox), nFixEf = as.integer(nFixEf), nDelta = as.integer(dim(delta)[2]), 
            nMCMC = as.integer(nMCMC), burnin = as.integer(burnin),
            X = as.double(X), FixEf = as.double(FixEf), Y = as.double(y),
            delta=as.integer(delta), FixEfInvRoot = as.double(FixEfInvRoot),
            sigEps = as.double(bestHyper$sigEps),
            sigBeta = as.double(bestHyper$sigBeta),
            a = as.double(bestHyper$a), b = as.double(bestHyper$b),    
            calcCI = as.integer(1), # now quantile CI calculation
            alphaRes  = as.double(rep(0, nFixEf)), betaRes = as.double(rep(0, nVox)), gammaRes = as.double(rep(0, nVox)),
            betaLower = as.double(rep(0, nVox)), betaUpper = as.double(rep(0, nVox)),
            alphaLower = as.double(rep(0, nFixEf)), alphaUpper = as.double(rep(0, nFixEf)),
            NAOK = TRUE)
  
  # results: should be self-explanatory... ;)
  return(list(hyper = bestHyper, 
              alpha = fit$alphaRes, 
              beta = fit$betaRes, 
              gamma = fit$gammaRes,
              fitted.values = FixEf %*% fit$alphaRes + apply(X,2, function(x){sum(x * fit$betaRes)}),
              betaCI  = list(lower = fit$betaLower, upper = fit$betaUpper),
              alphaCI = list(lower = fit$alphaLower, upper = fit$alphaUpper)
  ))
}