#' SDE plugin estimator solver
#' 
#' @description Calculates recursive estimator for given hazard estimates, integrand function and gradients. Assumes there are no tied event times.
#' 
#' @param n Total number of indiivduals
#' @param hazMatrix Matrix consisting of hazards(rows) and their increments(columns) along the same time scale
#' @param F_fun Integrand function for the differential equation system
#' @param gradientList \code{list} containing derivatives of \code{F_fun}
#' @param X0 Vector containing initial values for the parameter
#' @param V0 Matrix containing initial values for the variance
#' 
#' @author Pål Christie Ryalen \email{p.c.ryalen@@medisin.uio.no}
#' 
#' @return \code{list} containing the parameter estimate \code{X}, and its covariance estimates
#' 
#' 
#' @examples 
#' 
#' ###############################################################
#' #########################  Survival  ##########################
#' ###############################################################
#' library(timereg)
#' 
#' n <- 100
#' dfr <- data.frame(to = rexp(n,1),from=0,event=1)
#' 
#' fit <- survfit(Surv(from,to,event==1)~1,data=dfr)
#' 
#' times <- fit$time
#' dN <- fit$n.event
#' Y <- fit$n.risk
#' Y[1] <- n
#' 
#' dA <- dN/Y
#' 
#' # Function specification
#' F_fun_Survival <- function(x)-x
#' gradientListSurvival <- list(function(x)-1)
#' X0_Survival <- 1
#' V0_Survival <- 0
#' 
#' paramEst_survival <- pluginEstimate(n,dA,F_fun_Survival,gradientListSurvival,X0_Survival,V0_Survival)
#' 
#' KM <- cumprod(1 - dA)
#' Greenwood <- KM^2 * cumsum(dA^2)
#' 
#' plot(times,paramEst_survival$X,type="s",main="SDE plugin survival estimates",ylab="",xlab="time")
#' lines(times,paramEst_survival$X + 1.96*sqrt(paramEst_survival$covariance),type="s")
#' lines(times,paramEst_survival$X - 1.96*sqrt(paramEst_survival$covariance),type="s")
#' lines(seq(0,10,length.out=100),exp(-seq(0,10,length.out=100)),col=2)
#' legend("topright",c("SDE plugin estimates","Exact"),lty=1,col=c(1,2),bty="n")
#' 
#' plot(times,paramEst_survival$covariance,type="s",main="SDE plugin variance vs Greenwood variance",ylab="",xlab="time")
#' lines(times,Greenwood,type="s",col=4,lty=1)
#' legend("topright",c("SDE plugin estimates","Greenwood"),lty=1,col=c(1,4),bty="n")
#' 
#' #################################################################################
#' ############ Competing risks and Cumulative incidence(two states) ###############
#' #################################################################################
#' 
#' n <- 200
#' x1 <- rexp(n,1)
#' x2 <- rexp(n,1.3)
#' to.states <- ifelse(x1<x2,0,1)
#' 
#' dfr <- data.frame(from=0,to=c(x1[to.states==0],x2[to.states==1]),to.state=to.states)
#' dfr <- dfr[order(dfr$to),]
#' 
#' nrisk <- c(n,n:1)
#' dA1 <- c(0,1*(dfr$to.state==0))/nrisk
#' dA2 <- c(0,1*(dfr$to.state==1))/nrisk
#' 
#' hazMatrix <- rbind(dA1,dA2)
#' 
#' F_fun_cuminc <- function(X)rbind(c(X[2],0),c(-X[2],-X[2]))
#' gradientList_cuminc <- list( function(X)matrix(c(0,0,1,0),nrow=2),
#'                              function(X)matrix(c(0,0,-1,-1),nrow=2) )
#'                              
#' X0_cuminc <- c(0,1)
#' V0_cuminc <- matrix(0,nrow=2,ncol=2)
#' 
#' paramEst_cuminc <- pluginEstimate(n,hazMatrix,F_fun_cuminc,gradientList_cuminc,X0_cuminc,V0_cuminc)
#' 
#' times <- c(0,dfr$to)
#' 
#' plot(times,paramEst_cuminc$X[1,],type="s",ylab="",xlab="time",main="SDE plugin cumulative incidence estimate",ylim=c(0,0.7))
#' lines(times,paramEst_cuminc$X[1,] + 1.96*sqrt(paramEst_cuminc$covariance[1,]),type="s")
#' lines(times,paramEst_cuminc$X[1,] - 1.96*sqrt(paramEst_cuminc$covariance[1,]),type="s")
#' lines(seq(0,10,length.out = 1000),10/1000*cumsum(exp(-seq(0,10,length.out = 1000)*(1 + 1.3))),col=2)
#' legend("topright",c("SDE plugin estimates","Exact"),lty=1,col=c(1,2),bty="n")
#' 
#' 
#' ##########################################################################################
#' ####################  Relative survival(two different populations)  ######################
#' ##########################################################################################
#' 
#' n <- 300
#' t1 <- sort(rexp(n,1))
#' t2 <- sort(rexp(n,1.3))
#' times <- sort(c(0,t1,t2))
#' nrisk <- 300:1
#' dA1 <- 1*(t1 %in% times)/nrisk
#' dA2 <- 1*(t2 %in% times)/nrisk
#'
#' tmatch1 <- match(t1,times)
#' tmatch2 <- match(t2,times)
#'
#' hazMatrix <- matrix(0,nrow=2,ncol=length(times))
#' hazMatrix[1,tmatch1] <- dA1
#' hazMatrix[2,tmatch2] <- dA2
#' 
#' F_fun_RelSurv <- function(X)rbind(c(-X[1],X[1]),c(-X[2],0),c(0,-X[3]))
#' gradientList_RelSurv <- list(function(X)rbind(c(-1,0,0),c(0,-1,0),c(0,0,0)),
#'                              function(X)rbind( c(1,0,0), c(0,0,0), c(0,0,-1) ))
#'                              
#' X0_RelSurv <- c(1,1,1)
#' V0_RelSurv <- matrix(0,nrow=3,ncol=3)
#' 
#' 
#' paramEst_relsurv <- pluginEstimate(n,hazMatrix,F_fun_RelSurv,gradientList_RelSurv,X0_RelSurv,V0_RelSurv)
#' 
#' 
#' plot(times,paramEst_relsurv$X[1,],type="s",ylab="",xlab="time",main="SDE plugin relative survival estimate",ylim=c(-1,5.8),xlim=c(0,4))
#' lines(times,paramEst_relsurv$X[1,] + 1.96*sqrt(paramEst_relsurv$covariance[1,]),type="s")
#' lines(times,paramEst_relsurv$X[1,] - 1.96*sqrt(paramEst_relsurv$covariance[1,]),type="s")
#' lines(seq(0,10,length.out = 100),exp(seq(0,10,length.out = 100)*(1.3-1)),col=2)
#' legend("topright",c("SDE plugin estimates","Exact"),lty=1,col=c(1,2),bty="n")
#' 
#' @references Ryalen, P.C., Stensrud, M.J., Røysland, K.: \emph{Transforming cumulative hazards}, arXiv:1710.07422, to appear in Biometrika 2018.
#' 
#' @export



# Euler solver  ~~  Built-in lebesgue
pluginEstimate <- function(n,hazMatrix,F_fun,gradientList,X0,V0){
        if(is.null(dim(hazMatrix)))
                hazMatrix <- matrix(hazMatrix,nrow=1)
        
        numIncrements <- ncol(hazMatrix)
        X <- matrix(0,nrow=length(X0),ncol=numIncrements)
        V <- matrix(0,nrow=prod(dim(V0)),ncol=numIncrements)
        X[,1] <- X0
        V[,1] <- V0
        
        XRows <- length(X0)
        hazRows <- nrow(hazMatrix)
        
        numHaz <- length(gradientList)
        
        for(i in 2:numIncrements){
                X_last <- X[,i-1]
                # Hack for dealing with NA
                X_last[is.na(X_last)] <- X0[is.na(X_last)]
                
                FX_mat <- F_fun(X_last)
                X[,i] <- X_last + FX_mat %*% hazMatrix[,i]
                
                Vtemp <- matrix(0,nrow=XRows,ncol=XRows)
                Vn <- matrix(V[,i-1],nrow=XRows,ncol=XRows)
                
                # Hack for dealing with NA
                Vn[is.na(Vn)] <- V0[is.na(Vn)]
                for(j in 1:numHaz){
                        Vtemp <- Vtemp + (Vn %*% t(gradientList[[j]](X_last)) + gradientList[[j]](X_last) %*% Vn) * hazMatrix[j,i]
                }
                
                dW <- diag(hazRows) * hazMatrix[,i]^2
                
                Vn <- Vn + Vtemp + n*FX_mat %*% dW %*% t(FX_mat)
                V[,i] <- as.vector(Vn)
                
                
        }
        
        retList <- list(X=X,covariance=V/n)
        return(retList)
}
