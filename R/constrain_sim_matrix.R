#' Hard constraining of similarity matrix outside of global fit window
#'
#' This function takes in timeVectors from both runs, a global-fit object and
#' sample-length of window of no constraining. Outside of window, all elements
#' of matrix are equally weighted.
#' @param tRunAVec A numeric vector.
#' @param tRunBVec A numeric vector.
#' @param Fit An object of class "loess".
#' @param noBeef A parameter which defines the window size of no constraining.
#' @export
calcNoBeefMaskGlobal <- function(tRunAVec, tRunBVec, Fit, noBeef = 50){
  A1 <- tRunAVec[1]; A2 <- tRunAVec[length(tRunAVec)]
  B1 <- tRunBVec[1]; B2 <- tRunBVec[length(tRunBVec)]
  B1.p <- predict(Fit, A1)
  B2.p <- predict(Fit, A2)
  if(is.na(B1.p + B2.p)){mask <- matrix(0, nrow = length(tRunAVec), ncol = length(tRunBVec))}
  else{
    deltaTime <- (A2-A1)/(length(tRunAVec)-1)
    idx <- floor((B1.p - B1)/deltaTime) + 1
    leftTri1 <- c(idx-noBeef,0); rightTri1 <- c(idx+noBeef,0)
    idx <- ceiling((B2.p - B1)/deltaTime)
    leftTri2 <- c(idx-noBeef,length(tRunAVec)-1); rightTri2 <- c(idx+noBeef,length(tRunAVec)-1) # Considering our square starts at (0,0)
    mask <- matrix(1, nrow = length(tRunAVec), ncol = length(tRunBVec))
    for(x in 0:(length(tRunBVec)-1)){
      for(y in 0:(length(tRunAVec)-1)){
        lessThanLeftTri <- y-leftTri1[2]-
          (leftTri2[2]-leftTri1[2])*(x-leftTri1[1])/(leftTri2[1]-leftTri1[1])
        greaterThanRightTri <- y-rightTri1[2]-
          (rightTri2[2]-rightTri1[2])*(x-rightTri1[1])/(rightTri2[1]-rightTri1[1])
        if(lessThanLeftTri<=0 & greaterThanRightTri>=0){mask[y+1, x+1]<-0}
      }
    }
  }
  return(mask)
}

#' Soft constraining of similarity matrix outside of global fit window
#'
#' This function takes in datafile, precursor id, names of run pair. Based on
#' similarity measure, it calculates the similarity matrix. Outside of window,
#' elements of matrix are weighted proportional to distance from window-boundry.
#' @param tRunAVec A numeric vector.
#' @param tRunBVec A numeric vector.
#' @param Fit An object of class "loess".
#' @param noBeef A parameter which defines the window size of no constraining.
#' @export
calcNoBeefMaskGlobalWSlope <- function(tRunAVec, tRunBVec, Fit, noBeef =50){
  A1 <- tRunAVec[1]; A2 <- tRunAVec[length(tRunAVec)]
  B1 <- tRunBVec[1]; B2 <- tRunBVec[length(tRunBVec)]
  B1.p <- predict(Fit, A1)
  B2.p <- predict(Fit, A2)
  if(is.na(B1.p + B2.p)){mask <- matrix(0, nrow = length(tRunAVec), ncol = length(tRunBVec))}
  else{
    deltaTime <- (A2-A1)/(length(tRunAVec)-1)
    idx <- floor((B1.p - B1)/deltaTime) + 1
    leftTri1 <- c(idx-noBeef,0); rightTri1 <- c(idx+noBeef,0)
    idx <- ceiling((B2.p - B1)/deltaTime)
    leftTri2 <- c(idx-noBeef,length(tRunAVec)-1); rightTri2 <- c(idx+noBeef,length(tRunAVec)-1) # Considering our square starts at (0,0)
    mask <- matrix(NA, nrow = length(tRunAVec), ncol = length(tRunBVec))
    for(x in 0:(length(tRunBVec)-1)){
      for(y in 0:(length(tRunAVec)-1)){
        # distance(P1, P2, (X0, Y0)) <- -((Y2-Y1)X0 - (X2-X1)Y0 + X2*Y1 -
        # Y2*X1)/sqrt((Y2-Y1)^2 +(X2-X1)^2)
        distFromLeftBndry <- -((leftTri2[2]-leftTri1[2])*x -
                                 (leftTri2[1]-leftTri1[1])*y + leftTri2[1]*leftTri1[2] -
                                 leftTri2[2]*leftTri1[1])/sqrt((leftTri2[2]-leftTri1[2])^2 +
                                                                 (leftTri2[1]-leftTri1[1])^2)
        distFromRightBndry <- -((rightTri2[2]-rightTri1[2])*x -
                                  (rightTri2[1]-rightTri1[1])*y + rightTri2[1]*rightTri1[2] -
                                  rightTri2[2]*rightTri1[1])/sqrt((rightTri2[2]-rightTri1[2])^2 +
                                                                    (rightTri2[1]-rightTri1[1])^2)
        if(distFromLeftBndry<=0 & distFromRightBndry>=0){mask[y+1, x+1] <- 0
        } else if(distFromLeftBndry > 0) { mask[y+1, x+1] <- distFromLeftBndry
        } else if(distFromRightBndry < 0) {mask[y+1, x+1] <- -distFromRightBndry}
      }
    }
    # Above eqution is just from a parallel line: slope = (A2-A1)/(B2-B1) =
    # (y-A1)/(x- B1) Also, we are considering a very simple linear fit with
    # slope ~= 1. This will allow a fitting line to cross square only at
    # adjacent sides. If the slope comes out to be deviating too much from 1,
    # then this equation will fail and also in that case better to fit
    # non-linear boundaries for No-Beef region.
  }
  return(mask)
}

#' Hard constraining of similarity matrix outside of global fit window
#'
#' This function takes in timeVectors from both runs, a global-fit object and
#' sample-length of window of no constraining. Outside of window, all elements
#' of matrix are equally weighted.
#' @param s A numeric similarity matrix.
#' @param mask A numeric matrix with weights for constraining.
#' @param unitConstrain Constraining value for unit weight.
#' @export
constrainSimilarity <- function(s, mask, unitConstrain){
  return(s + unitConstrain*mask)
}
