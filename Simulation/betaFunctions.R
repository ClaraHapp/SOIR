#' Simulate smooth beta function
#' 
#' @param x1,x2 Vectors forming the rectangle, where the new beta function should be evaluated
#' @param random A logical variable, indicating wheter the weights should be 
#'   generated randomly (\code{random = TRUE}) or following a fixed scheme 
#'   (\code{random = FALSE}). Defaults to \code{FALSE}.
#'   
#' @return A new simulated image beta.
smoothBeta <- function(x1,x2, random = FALSE)
{
  grid <- expand.grid(x1, x2)
  sideLengths <- c(length(x1), length(x2))
  
  means <- {if(random) {matrix(runif(6), nrow = 3)} else {rbind(c(.15, .5), c(.7, .2), c(0.8, 0.7))}}
  
  beta <- array(-0.05+ 0.09*mvtnorm::dmvnorm(x = grid, mean = means[1,], sigma=diag(c(0.1, 0.1))) -
                  0.01*mvtnorm::dmvnorm(x = grid, mean = means[2,], sigma=diag(c(0.02, 0.03))) +
                  0.05*mvtnorm::dmvnorm(x = grid, mean = means[3,], sigma = matrix(c(0.05,0.02, 0.02, 0.1), nrow = 2)),
                dim = sideLengths)
  return(beta)
}


#' Simulate sparse beta function
#' 
#' @param x1,x2 Vectors forming the rectangle, where the new beta function should be evaluated
#' @param random A logical variable, indicating wheter the weights should be 
#'   generated randomly (\code{random = TRUE}) or following a fixed scheme 
#'   (\code{random = FALSE}). Defaults to \code{FALSE}.
#'   
#' @return A new simulated image beta.
sparseBeta <- function(x1,x2, random = FALSE)
{
  grid <- expand.grid(x1, x2)
  sideLengths <- c(length(x1), length(x2))
  
  means <- {if(random) {matrix(runif(4), nrow = 2)} else {rbind(c(.2, .3), c(0.4, 0.8))}}
  
  beta <- array(.006 * mvtnorm::dmvnorm(x = grid, mean = means[1,], sigma=diag(c(.0025, 0.0015))) -
                  0.003 * mvtnorm::dmvnorm(x = grid, mean = means[2,], sigma = matrix(c(0.002,-0.001, -0.001, 0.001), nrow = 2, ncol = 2)),
                dim =sideLengths)
  beta[abs(beta) < 0.25] <- 0
  return(beta/2)
}

#' Beta function based on functional data
#' 
#' @param functions A funData object, representing the (principal component)
#'   functions to be used for generating beta
#' @param random A logical variable, indicating wheter the weights should be
#'   generated randomly (\code{random = TRUE}) or following a fixed scheme
#'   (\code{random = FALSE}). Defaults to \code{FALSE}.
#'   
#' @return A new simulated image beta.
pcaBeta <- function(functions, random = FALSE)
{
  if(random)
  {
    K <- sample(nObs(functions),1)
    weights <- runif(K, min = -1, max = 1)
  }
  else
  {
    K <- 5
    weights <- (-1)^(1:K)*exp(-(1:K)/K)
  }
  
  
  return(MFPCA::ttv(functions@X[1:K,,,drop = FALSE], list(v = weights), dim = 1))
}


#' Bumpy beta function
#' 
#' The authors thank P. T. Reiss for providing the code as in Reiss et al. 2015
#' 
#' @param seed A random seed that is set before simulating the bumps.
#' @param n Sidelength of the images. Should be a power of 2. Defaults to 64.
#'   
#' @return A new simulated image beta.
bumpsBeta <- function(seed, n = 64)
{
  # bummpy wave
  xcol <- matrix(rep(1:n, n), ncol=n)
  xrow <- matrix(rep(1:n, n), ncol=n, byrow=TRUE)
  
  set.seed(seed)
  t <- floor(matrix(runif(22), ncol=2)*n)
  h2 <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
  w <- seq(1,50,length=length(h2))
  bumps <- array(0, dim = c(n,n))
  for (i in 1 : length(h2)){
    bumps <- bumps + h2[i] * pmax(0, 1-abs(sqrt((xcol - t[i,1])^2 + (xrow - t[i,2])^2)/w[i]))^4
  }
  
  return(bumps)
}


# Simulate beta functions
x1 <- seq(0, 1, length.out = 64)
x2 <- seq(0, 1, length.out = 64)

beta <- list(
  smooth = smoothBeta(x1, x2, random = FALSE),
  sparse = sparseBeta(x1, x2, random = FALSE),
  bumpy = bumpsBeta(seed = 24, n = 64)/50,
  pca =  pcaBeta(<PCA2D$functions>, random = FALSE) # note that PCA2D depends on the ADNI data and cannot be published.
)
