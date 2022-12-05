#' SDE plugin estimator solver
#'
#' @description Calculates recursive estimator for given hazard estimates, integrand function and gradients.
#'
#' @param n Total number of indiivduals
#' @param hazMatrix Matrix consisting of hazards(rows) and their increments(columns) along the same time scale
#' @param F_fun Integrand function F = (F_1,F_2,...) for the differential equation system
#' @param JacobianList The Jacobian matrices of F_1, F_2, ... organized in a \code{list}
#' @param X0 Matrix containing initial values for the parameter
#' @param V0 Matrix containing initial values for the variance
#' @param isLebesgue (Optional, to improve efficientcy) Provide the index of \code{A} (e.g. 2) that is a regular \code{dt} integral, i.e. not a cumulative hazard.
#'
#' @author Pål Christie Ryalen \email{p.c.ryalen@@medisin.uio.no}
#'
#' @return \code{list} containing the parameter estimate \code{X}, and its covariance estimates \code{V/n}
#'
#'
#' @example inst/examples/pluginEstimate_example.R
#'
#' @references Ryalen, P.C., Stensrud, M.J., Røysland, K.: \emph{Transforming cumulative hazards}, arXiv, to appear in Biometrika 2018.
#'
#' @export



# Euler solver
pluginEstimate <- function(n,hazMatrix,F_fun,JacobianList,X0,V0,isLebesgue = NULL){
        if(is.null(dim(hazMatrix)))
                hazMatrix <- matrix(hazMatrix,nrow=1)

        if( !( "matrix" %in% class(V0) & "matrix" %in% class(F_fun(X0)) &
                 "matrix" %in% class(JacobianList[[1]](V0)) & "matrix" %in%  class(hazMatrix) ) )
        stop("Please provide all input on matrix form")

        numIncrements <- ncol(hazMatrix)
        X <- matrix(0,nrow=length(X0),ncol=numIncrements)
        # V <- matrix(0,nrow=prod(dim(V0)),ncol=numIncrements)

        V <- array(0,dim = c(dim(V0),numIncrements))
        X[,1] <- X0
        V[,,1] <- V0

        XRows <- length(X0)
        hazRows <- nrow(hazMatrix)

        numHaz <- length(JacobianList)

        for(i in 2:numIncrements){
                X_last <- X[,i-1]

                # Euler step for parameter X
                FX_mat <- F_fun(X_last)
                X[,i] <- X_last + FX_mat %*% hazMatrix[,i]

                Vtemp <- matrix(0,nrow=XRows,ncol=XRows)
                Vn <- V[,,i-1]

                # Euler step for covariance V
                for(j in 1:numHaz){
                        Vtemp <- Vtemp + (Vn %*% t(JacobianList[[j]](X_last)) + JacobianList[[j]](X_last) %*% Vn) * hazMatrix[j,i]
                }

                # If 'dt' integrals are included, efficientcy can be improved by setting values equal to 0
                dB <- hazMatrix[,i]
                dB[isLebesgue] <- 0

                # The outer product of the increments
                dW <- dB %o% dB

                Vn <- Vn + Vtemp + n*FX_mat %*% dW %*% t(FX_mat)
                V[,,i] <- Vn


        }

        retList <- list(X=X,covariance=V/n)
        return(retList)
}
