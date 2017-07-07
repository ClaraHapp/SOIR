#' Measure the smoothness via a smoothness penalty
#' 
#' This function calculates the smoothness penalty for given differences, i.e. 
#' it vectorizes the image to a vector x and calculates x' P x for a smoothness 
#' penalty matrix based on differences. According to Rayleighs theorem, 
#' lambda_min(P) <= x' P x /(x'x)  <= lambda_max(P) with equalities for the
#' corresponding eigenvectors. Since P is a difference matrix, it is positive
#' semidefinite and has no full rank, i.e. lambda_min = 0. Hence
#' 0 <= x' P x /lambda_max x' x) <= 1
#' can be seen as a smoothness measure.
#' Values close to 0 indicate smoothness, values close to 1 indicate maximal non-smoothness. 
smooth <- function(img, differences, maxEig)
{
  return( (sum(diff(img, differences = differences)^2)+ sum(diff(t(img), differences = differences)^2))/(maxEig*sum(img^2)))
}


#' Calculate  a measure of sparsity for an image based on the gini index
#' 
#' The result is between 0 and 1. In the case of sparsity (extreme inequality), the result approaches 0 and 1 in the case of no sparsity, i.e. extreme equality
sparsity <- function(img)
{
  return(  2*sum(sort( abs(img))*( length(img) - 1: length(img) + 1/2))/( length(img)*sum( abs(img))))
}

#' Transform an image to the wavelet domain and calculate the measure of sparsity of the Wavelet coefficients
sparsityWave <- function(img)
{
  return(sparsity(unlist(wavethresh::imwd(img)[-(1:6)]))) # extract all wavelet coefficients and calculate sparsity measure for them
}


#' Project an image on a basis
#' 
#' This function projects an image \code{img} on some vectorized basis images and rescales by the norm of the image.
#' If the image lies completely in the span of the basis images, the measure is equal to 0. If the image is orthogonal to all basis images, the measure equals 1.
#' 
#'  @param img An image of dimensions dxd
#'  @param basis A system of K vectorized basis functions, supplied as an array of dimension K x d^2 
basisProjection <- function(img, basis)
{
  return(sum(lm(as.vector(img) ~ basis)$residuals^2) / sum(img^2))
}

#' Project image in wavelet space and calculate relative residuals
waveletProj <- function(img)
{
  return(sum((wavethresh::imwr(wavethresh::imwd(img)) - img)^2) / sum(img^2))
}

# Project image on spline space and calculate relative residuals
splineProj <- function(img, m, k)
{
  d <- dim(img)
  
  xx <- rep(1:d[1], times = d[2])
  yy <- rep(1:d[2], each = d[1])
  
  g <- mgcv::gam(as.vector(img) ~ te(xx,yy, bs = "ps", m = m, k = k), method = "REML")
  
  return(sum(g$residuals^2)/sum(img^2))
}

#' Calculate a measure of prior variability in the GMRF model based on the Kullback-Leibler divergence
#' 
#' @param aPri,bPri The parameters of the inverse gamma prior
#' @param aPost,bPost The parameters of the full conditional
#' 
#' @return A measure of the prior variability
prior <- function(aPri, bPri, aPost, bPost)
{
  D <- aPost * log(bPost) - aPri * log(bPri)  + lgamma(aPri) - lgamma(aPost) + ( aPri - aPost) * (log(bPost) - digamma(aPost)) + (bPri - bPost) * aPost/bPost
  
  return(1 - exp(-D/10))
}
