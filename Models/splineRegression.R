#' Estimate Scalar-on-Image Regression using a Penalized Spline Representation
#' 
#' This function fits a scalar-on-image regression model based on a penalized 
#' spline representation of the unknown coefficient image.
#' 
#' Requirements: This function requires the package \pkg{mgcv}.
#' 
#' @param XnoMean The demeaned iamge data, an array of dimension \code{N x S_1 x
#'   S_2}, containing \code{N} images of dimension \code{S_1 x S_2}.
#' @param y A vector of length \code{N} containing the response values.
#' @param covt An optional matrix containing the intercept and additional 
#'   covariates in row-wise manner.
#' @param m,k Further parameters passed to \pkg[mgcv]{gam}, specifying the order
#'   of the penalty and the number of basis functions used.
#'   
#' @return \item{alphaHat}{The estimated coefficients for the scalar 
#'   covariates.} \item{betaHat}{The estimated coefficient image.} 
#'   \item{betaCoef}{The spline coefficients that yield the coefficient image.} 
#'   \item{fitted.values}{The vector of fitted response values.}
#'   \item{seAlpha}{Standard errors for the scalar coefficients.}
#'   \item{seBeta}{A matrix containing pointwise standard errors for the
#'   coefficient image.}
splineRegression <- function(XnoMean, y, covt = NULL, m, k)
{
  d <- dim(XnoMean) # original dimensions
  
  # reshape X
  dim(XnoMean) <- c(d[1], prod(d[-1]))
  
  # x and y coordinates of images
  coord <- expand.grid(x = 1:d[2], y = 1:d[3])
  xx <- matrix(coord$x, nrow = d[1], ncol = prod(d[-1]), byrow = TRUE)
  yy <- matrix(coord$y, nrow = d[1], ncol = prod(d[-1]), byrow = TRUE)
  
  intercept <- matrix(1, nrow = d[1], ncol = 1)
  
  if(is.null(covt))
    covt <- intercept
  else
    covt <- cbind(intercept, covt)
  
  fit <- mgcv::gam(y ~ covt -1 + te(xx, yy, by = XnoMean, bs = "ps", m = m, k = k), method = "REML")
  
  betaHat <- predict.gam(fit, newdata = list(xx = coord$x, yy = coord$y, XnoMean = rep(1,prod(d[-1])), covt = matrix(0, nrow = prod(d[-1]), ncol = dim(covt)[2])), type = "terms", se.fit = TRUE)
  
  return(list(alphaHat = fit$coef[1:dim(covt)[2]], betaHat = matrix(betaHat$fit[,2], nrow = d[2], ncol = d[3]),
              betaCoef = fit$coef[-(1:dim(covt)[2])], fitted.values = fit$fitted.values,
              seAlpha = summary(fit)$se[1:dim(covt)[2]],
              seBeta = matrix(betaHat$se.fit[,2], nrow = d[2], ncol = d[3])))
}