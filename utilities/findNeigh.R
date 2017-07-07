#' Find neighbouring pixels/voxels for a given index, possibly within a mask
#' 
#' @param l The index for which the neighbours should be found
#' @param maskInd # the indices of the mask
#' @param match # the indices in the full image
#' @param maskDim Dimensions of the mask
#' @param dim A character string, specifying the dimension. Can be either
#'   \code{"2D"} (the default) or \code{"3D"}
#'   
#' @return The indices of neighbouring pixels/voxels in maskInd
findNeigh <- function(l, maskInd, match, maskDim, dim = "2D")
{
  neigh <- switch(dim,
                  "2D" =  c(-1,+1,-maskDim[1],+maskDim[1]),
                  "3D" =  c(-1,+1 , -maskDim[1],+maskDim[1],-maskDim[1]*maskDim[2], +maskDim[1]*maskDim[2] ),
                  stop("Function findNeigh: dim must be '2D' or '3D'!"))
  
  # theoretical neighbourhood of voxel l in vectorized(!) full mask
  deltaFull = maskInd[l] +  neigh
  
  # set values ouside mask to NA
  deltaFull[deltaFull <= 0] <- NA
  deltaFull[deltaFull > max(maskInd)] <- NA
  
  # check borders
  deltaFull[1:2][((deltaFull[1:2]-1)%%(maskDim[1]*maskDim[2]))%/%maskDim[1] != ((maskInd[l]-1)%%(maskDim[1]*maskDim[2]))%/%maskDim[1]] <- NA
  deltaFull[3:4][ (deltaFull[3:4]-1)%/%(maskDim[1]*maskDim[2])  != (maskInd[l]-1)%/%(maskDim[1]*maskDim[2])] <- NA
  
  if(dim == "3D")
    deltaFull[5:6][(deltaFull[5:6]-1)%%(maskDim[1]*maskDim[2]) != (maskInd[l]-1)%%(maskDim[1]*maskDim[2])] <- NA
  
  return(match[deltaFull]) # indices of neighbourhood voxels in maskInd
}