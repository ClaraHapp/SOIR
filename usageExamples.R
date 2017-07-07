#### Usage examples ####

#' Assume the demeaned image covariates to be given in an array XnoMean of 
#' dimension N x S x S, where N is the number of observations and the sidelength
#' S is a power of 2.
#' 
#' The response vector y is a numeric vector of length N.
#' 
#' Potential scalar covariates without an intercept are stored in the N x (p-1) 
#' matrix covt. The matrix FixEf contains covt with an additional column for the
#' intercept.

#### Source all relevant files ####
source("Models/*.R")


#### Spline Regression ####
SplineReg <- splineRegression(XnoMean = XnoMean, y = y, covt = covt, # data
                              m = 2, # cubic B-spline with 2nd order penalty
                              k = 15) # number of basis functions for each direction


#### FPCR ####
FPCR <- refund::fpcr(y = y, xfuncs = XnoMean, covt = covt, # data
                     ncomp = c(5, 10, 25, 50, 100, 150),  # number of principal components to use (optimal value is chosen in CV)
                     nbasis = list(v1 = 15, v2 = 15), # number of B-spline basis functions in each direction
                     family = "gaussian", # response is normally distributed
                     method = "REML") # smoothing parameter is found via REML


#### PCR 2D ####
# calculate PCA first
PCA2D <- MFPCA::fcptpaBasis(funData::funData(argvals = list(1:dim(XnoMean)[2],1:dim(XnoMean)[3]), 
                                             X = XnoMean), # create funData object containing the images
                            npc = 25, # calculate 25 PCs
                            smoothingDegree = c(2,2), # penalize 2nd differences in both directions
                            alphaRange = list(v = c(10e-4, 1e2), w = c(1e-4, 1e2))) # choose optimal smoothing parameters in (1e-4, 1e2) for both directions 


# Scalar-on-Image regression
PCR2D <- PCRCV(scores = PCA2D$scores, # scores from the PCA
               y = y, FixEf = FixEf, # further data
               functions = PCA2D$functions, # principal component functions
               nGroups = 5, # number of groups for cross-validation
               kVals = c(1,5,10,15,20,25), # number of principal components to use (optimal value is chosen in CV)
               CI = FALSE) # do not calculate confidence bands


#### WCR ####
WCR <- refund.wave::wcr(y, xfuncs = XnoMean, covt = covt, # data
                        min.scale = 3, # resolution level
                        nfeatures = c(10, 25, 50, 100, 250, 500, 1000), # number of wavelet-coefficients to retain (optimal value is chosen in CV)
                        ncomp = c(5, 10, 15, 25, 50, 75), # number of principal components to use (optimal value is chosen in CV)
                        method = "pcr")


#### WPLS ####
WPLS <- refund.wave::wcr(y, xfuncs = XnoMean, covt = covt, # data
                         min.scale = 3, # resolution level
                         nfeatures = c(10, 25, 50, 100, 250, 500, 1000), # number of wavelet-coefficients to retain (optimal value is chosen in CV)
                         ncomp = c(5, 10, 15, 25, 50, 75), # number of partial least squares components to use (optimal value is chosen in CV)
                         method = "pls")


#### WNET ####
WNET <- refund.wave::wnet(y, xfuncs = XnoMean, covt = covt, # data
                          min.scale = 3,  # resolution level
                          nfeatures = c(10, 25, 50, 100, 250, 500, 1000),  # number of wavelet-coefficients to retain (optimal value is chosen in CV)
                          alpha = c(0, 0.25, 0.5, 0.75, 1), # mixing parameter, 0 = RIDGE, 1 = LASSO (optimal value is chosen in CV)
                          store.glmnet = TRUE) # save the whole glmnet result


#### Extra variables for GMRF models ####
# For sparse GMRF
dyn.load("C/mainGibbs_HyperparamsFixed.so")
# For GMRF
dyn.load("C/mainGibbs_GMRF.so")


nPix <- prod(dim(XnoMean)[2:3])

# Find indices of neighbouring pixels
delta <- t(apply(matrix(1:nPix), 1, findNeigh, # apply function findNeigh to each pixel
                 maskInd = 1:nPix, match=1:nPix, maskDim = dim(XnoMean)[2:3], dim = "2D")) # parameters of findNeigh

# Reshape images to matrix that contains the vectorized image covariates in a row-wise format
Xreshaped <- aperm(XnoMean, c(2,3,1))
dim(Xreshaped) <- c(prod(sideLengths), dim(XnoMean)[1])

####  SparseGMRF ####
SparseGMRF <- GoldsmithCV(hyper = list(sigEps = c(1e-5, 1e-3, 1e-1), # hyperparameter values to chose in CV for variance of error term ...
                                       sigBeta = c(1e-5, 1e-3, 1e-1), # ... variance of beta ...
                                       a = c(-4, -2, -0.5), # ... Ising parameter a ...
                                       b = c(0.1, 0.5, 1.5)), # ... and Ising parameter b.
                          nGroups = 5, # number of groups to use in cross-validation
                          X = Xreshaped, FixEf = FixEf, y = y, # data with reshaped X
                          delta = delta, # neighbourhood information
                          K = 250, # number of MCMC iterations
                          nMCMC = 1, # no thinning, i.e. save every iteration
                          burnin = 100) # first 100 iterations are discarded as burnin


#### GMRF ####  
GMRF <- .C("allGibbsGMRF", # call C function
           K = as.integer(5000), # number of MCMC iterations
           nPat = as.integer(length(y)), # number of observations
           nVox = as.integer(nPix), # number of pixels/voxels per image
           nFixEf = as.integer(dim(FixEf)[2]), # number of scalar covariates (incl. intercept)
           nDelta = as.integer(dim(delta)[2]), # maximal number of neighbours
           nMCMC = as.integer(20), # thinning, i.e. save only each 20. iteration
           burnin = as.integer(1000), # discard first 1000 iterations as burn-in
           X = as.double(Xreshaped), # matrix of reshaped data
           FixEf = as.double(FixEf), # scalar covariates (incl. intercept)
           Y = as.double(y), # response
           delta = as.integer(delta), # neighbourhood information
           FixEfInvRoot = as.double(chol(solve(t(FixEf)%*%FixEf))), # auxiliary parameter
           shapeEps = as.double(1), # prior shape for the variance parameter sigEps
           scaleEps = as.double(1), # prior scale for the variance parameter sigEps
           shapeBeta = as.double(1), # prior shape for the variance parameter sigBeta
           scaleBeta = as.double(1), # prior scale for the variance parameter sigBeta
           calcCI = as.integer(0), # do not calculate CIs
           alphaRes  = as.double(rep(0, dim(FixEf)[2])), # allocate result vector for scalar coefficients alpha
           betaRes = as.double(rep(0, nPix)), # allocate result vector for image coefficient beta
           sigEpsRes = as.double(0), # result variable for the variance parameter sigEps
           sigBetaRes = as.double(0), # result variable for the variance parameter sigBeta
           betaLower = as.double(0), # variables for confidence bands, not used here
           betaUpper = as.double(0),
           alphaLower = as.double(0), 
           alphaUpper = as.double(0),
           NAOK = TRUE) # allow for NA values, that are included in neighbourhood structure delta

# construct fitted values
GMRF$fitted.values <- FixEf%*%GMRF$alphaRes + apply(XnoMean,1,function(x){sum(x * GMRF$betaRes)})