
  ## tseriesEntropy
  ## Entropy based tests of serial dependence and nonlinearity
  #
  #  The authors of this software is
  #
  #  Simone Giannerini, Copyright (c) 2009.
  #
  #  Permission to use, copy, modify, and distribute this software for any
  #  purpose without fee is hereby granted, provided that this entire notice
  #  is included in all copies of any software which is or includes a copy
  #  or modification of this software and in all copies of the supporting
  #  documentation for such software.
  #
  #  This program is free software; you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation; either version 2 of the License, or
  #  (at your option) any later version.
  #
  #  This program is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #  GNU General Public License for more details.
  #
  #  A copy of the GNU General Public License is available at
  #  http://www.r-project.org/Licenses/

 ## ***************************************************************************************************

 setClass("Srho.ts", contains="Srho",
         representation(method ="character", bandwidth="character"))

## ***************************************************************************************************
setMethod ("show" , "Srho.ts",
    function(object){
    out <- object@.Data
    names(out) <- object@lags;
    n <- length(out);
    lag.max <- object@lags[n]
    cat (" Srho computed on", lag.max, "lags \n")
    cat (" ------------------------------------------------ \n")
    print(out)
    cat (" ------------------------------------------------ \n")
    cat (" Data type               :" , object@data.type , "\n")
    cat (" Stationary version      :" , object@stationary , "\n")
    cat (" Computation method used :" , object@method , "\n")
    cat (" Bandwidth method used   :" , object@bandwidth , "\n")
    if(length(object@notes)>0){
    cat (" Additional notes        :" , object@notes, "\n");
    }
    cat ("\n")
    }
)

## ***************************************************************************************************

Srho.ts <- function(x,y, lag.max = 10, bw = c("reference",
"mlcv", "lscv", "scv", "pi"), bdiag=TRUE, method =c("integral","summation"), plot = TRUE, tol=1e-03,...){

    n <- length(x);
    if (missing(lag.max)) lag.max = round(n/4);
    if((lag.max >= n)||(lag.max < 1)) stop('incorrect number of lags')

    bw     <- match.arg(bw)
    method <- match.arg(method)
#    if(stationary=TRUE) stop("Stationary version not yet implemented for continuous data");
    if(!is.numeric(x)) stop('input series must be numeric')
    if (missing(y)||is.null(y)){ # RICORDARSI DI AGGIORNARE GLI HELP *********
        return(Srho.ts.uni(x,lag.max=lag.max,bw=bw,bdiag=bdiag,method=method, plot=plot, tol=tol,...))
    } else {
        if(!is.numeric(y)) stop('input series must be numeric')
        return(Srho.ts.biv(x,y,lag.max=lag.max,bw=bw,bdiag=bdiag,method=method, plot=plot, tol=tol,...))
    }
}

## **************************************************************************************************

# $Id: Entropy_lib.R,v 2.14 2002/06/23 17:14:05 jracine Exp jracine $'
# Written by J. Racine email: jracine@maxwell.syr.edu

# R code defining Srho described in Granger, Maasoumi, and Racine JTSA
# "A Dependence Metric for Possibly Nonlinear Processes".
# Finally, univariate and bivariate kernel density routines
# and cross-validation routines for data-driven bandwidth selection.

## ***************************************************************
## Extensively modified by Simone Giannerini
## email: simone.giannerini@unibo.it
## August 2005 -
## *******************************************************************
## All the bottleneck parts of the program related to Srho
## have been rewritten in F90 and are dynamically loaded
## Some R parts have been rewritten in a cleaner and more portable way
## *******************************************************************
## May 2011:
## the program uses the package "cubature" instead of "adapt" for
## multidimensional integration
## *******************************************************************

# maxpts <- 1e06 # Maximum number of function evaluations
# tol <- 1e-03   # Relative error '

##     ************************************************************************

Srho.func <- function(x1, x2, h1,h2, Hbiv, method = c("integral", "summation"),
tol=1e-03,...) {

    if(length(x1) != length(x2))  stop("Input vectors differ in length.\n")
    if((h1 <= 0) || (h2 <= 0))    stop("Bandwidths must be positive.\n")

    method <- match.arg(method)

    if(method == "integral") {
        # For range of integration, use min - range and max + range for each
        # variable
        lo.default <- c(3*min(x1)-2*max(x1),3*min(x2)-2*max(x2)) # Lower vertices (min -        2range)
        up.default <- c(3*max(x1)-2*min(x1),3*max(x2)-2*min(x2)) # Upper vertices (max +        2range)
        xb                 <- cbind(x1,x2)
        N                  <- nrow(xb)
        storage.mode(xb)   <- 'double'
        Hinv               <- solve(Hbiv)
        storage.mode(Hinv) <- 'double'

         Srho.integrand2 <- function(x) {
            storage.mode(x) <- 'double'
            nx              <- NCOL(x)
           # SUBROUTINE SRhointegrand2(X,nX,xb,N,h1,h2,Hinv,SINT)
            res <- .Fortran("srhointegrand2",x,as.integer(nx),xb,as.integer(N),
            as.double(h1),as.double(h2),Hinv,SINT=double(nx),PACKAGE='tseriesEntropy')$SINT;
            return(matrix(res,ncol=nx))
         }
         return(0.5*hcubature(f=Srho.integrand2, lowerLimit=lo.default,
         upperLimit=up.default, tol = tol,vectorInterface = TRUE,
         fDim = 1, absError=0,...)$integral)
    } else if(method == "summation") {
        h1.biv <- Hbiv[1,1]
        h2.biv <- Hbiv[2,2]
        Srho.sum <-
        .Fortran("srhosum", as.double(x1),as.double(x2), as.integer(length(x1)),
            as.double(h1),as.double(h2), as.double(h1.biv),as.double(h2.biv),S=double(1),PACKAGE='tseriesEntropy')$S;
        return(Srho.sum)
   }
}

## *****************************************************************************

Srho.func.old <- function(x1, x2, h1, h2, Hbiv, method = c("integral",
"summation"), tol=1e-03,...) {

    if(length(x1) != length(x2))  stop("Input vectors differ in length.\n")
    if((h1 <= 0) || (h2 <= 0))    stop("Bandwidths must be positive.\n")

    method <- match.arg(method)
    h1.biv <- Hbiv[1,1]
    h2.biv <- Hbiv[2,2]
    if(method == "integral") {
        # For range of integration, use min - range and max + range for each
        # variable
        lo.default <- c(3*min(x1)-2*max(x1),3*min(x2)-2*max(x2)) # Lower vertices (min - 2range)
        up.default <- c(3*max(x1)-2*min(x1),3*max(x2)-2*min(x2)) # Upper vertices (max + 2range)
        # Here we define a function which is the integrand of Srho when the
        # densities are estimated using a Gaussian kernel
         Srho.integrand <- function(x) {
            Srho.integrand <- .Fortran("srhointegrand",as.double(x),as.double(x1),as.double(x2),as.integer(length(x1)),
            as.double(h1),as.double(h2),as.double(h1.biv),as.double(h2.biv),SINT=double(1),PACKAGE='tseriesEntropy')$SINT;
            return(Srho.integrand)
         }
         return(0.5*hcubature(f=Srho.integrand, lowerLimit=lo.default, upperLimit=up.default, tol = tol,
         fDim = 1, absError=0,...)$integral)
    } else if(method == "summation") {

        Srho.sum <- .Fortran("srhosum",as.double(x1),as.double(x2),as.integer(length(x1)),
            as.double(h1),as.double(h2),as.double(h1.biv),as.double(h2.biv),S=double(1),PACKAGE='tseriesEntropy')$S;
        return(Srho.sum)
   }
}
## **************************************************************************************************

Srho.ts.uni <- function(x, lag.max = 10, bw = c("reference",
"mlcv", "lscv", "scv", "pi"),bdiag=TRUE, method =c("integral","summation"), plot = TRUE, tol=1e-03,...) {
    n <- length(x)
    S      <- rep(NA,lag.max)
    bw     <- match.arg(bw)
    method <- match.arg(method)
    bw.fun.uni <- paste(bw,"uni",sep=".")
    bw.fun.biv <- paste(bw,"biv",sep=".")
    h1   <- do.call(bw.fun.uni,list(x));
    h2   <- h1
    for(k in 1:lag.max){
        x.lag0 <- x[(1+k):n]
        x.lagk <- x[1:(n-k)]
        Hbiv   <- do.call(bw.fun.biv,list(x.lag0,x.lagk,bdiag));
        S[k]   <- Srho.func(x.lag0,x.lagk,h1,h2,Hbiv,method,tol,...)
    }
    out <- new("Srho.ts")
    out@.Data      <- S
    out@lags       <- 1:lag.max
    out@stationary <- FALSE
    out@data.type  <- "continuous"
    out@bandwidth  <- bw
    out@method     <- method
    if (plot){
        plot(out)
        return(invisible(out))
    }
    else return(out)
}

## **************************************************************************************************

Srho.ts.biv <- function(x, y, lag.max = 10, bw = c("reference",
"mlcv", "lscv", "scv", "pi"),bdiag=TRUE, method = c("integral","summation"),plot=TRUE, tol=1e-03,...){

    if(length(x) != length(y))  stop("Input vectors differ in length.\n")

    n <- length(x)
    bw     <- match.arg(bw)
    method <- match.arg(method)
    S      <- rep(NA,(2*lag.max+1))
    bw.fun.uni <- paste(bw,"uni",sep=".")
    bw.fun.biv <- paste(bw,"biv",sep=".")
    h1   <- do.call(bw.fun.uni,list(x));
    h2   <- do.call(bw.fun.uni,list(y));
    Hbiv <- do.call(bw.fun.biv,list(x,y,bdiag));

    S[lag.max+1] <- Srho.func(x,y,h1,h2,Hbiv,method,tol,...)

    for(k in 1:lag.max) {
        x1     <- x[1:(n-k)];
        y1     <- y[(k+1):n];
        Hbiv   <- do.call(bw.fun.biv,list(x1,y1,bdiag));
        S[lag.max+1+k] <- Srho.func(x1,y1,h1,h2,Hbiv,method,tol,...)
#
        x1     <- x[(k+1):n];
        y1     <- y[1:(n-k)];
        Hbiv   <- do.call(bw.fun.biv,list(x1,y1,bdiag));

        S[lag.max+1-k] <- Srho.func(x1,y1,h1,h2,Hbiv,method,tol,...)
    }
    out <- new("Srho.ts")
    out@.Data      <-  S
    out@lags       <- -lag.max:lag.max
    out@stationary <- FALSE
    out@data.type  <- "continuous"
    out@bandwidth  <- bw
    out@method     <- method
    if (plot) {
        plot(out)
        return(invisible(out))
    }
    else return(out)
}
## **************************************************************************************************

# Next we implement cross-validation for univariate and bivariate
# density estimation using second-order Gaussian kernels as is used in
# the metrics defined above.

# Kullback-Leibler cross-validation for kernel density estimation,
# second order Gaussian kernel (likelihood cross-validation)

mlcv.uni <- function(x) {
    DMACH <- rep(0,4)
    DMACH[1]<-.Machine$double.eps
    DMACH[2]<-.Machine$double.neg.eps
    DMACH[3]<-.Machine$double.xmin
    DMACH[4]<-.Machine$double.xmax
    x <- x
    n <- length(x)
    kdenest.mlcv <- function(h) {
       kdenest.mlcv <- .Fortran("kdenestmlcv",as.double(x),as.integer(length(x)),as.double(h),F=double(1),as.double(DMACH),PACKAGE='tseriesEntropy');
       return(kdenest.mlcv$F);
    }
    return(nlm(kdenest.mlcv, 1.06*(min(sd(x),IQR(x)/1.34))*n^{-1/5})$estimate)
}
## **************************************************************************************************

mlcv.biv <- function(x,y,bdiag) {
    DMACH <- rep(0,4)
    DMACH[1]<-.Machine$double.eps
    DMACH[2]<-.Machine$double.neg.eps
    DMACH[3]<-.Machine$double.xmin
    DMACH[4]<-.Machine$double.xmax
    x <- x
    y <- y
    n <- length(x);
    kdenest.mlcv <- function(h) {
        kdenest.mlcv <- .Fortran("kdenestmlcvb",as.double(x),as.double(y),as.integer(length(x)),as.double(h),F=double(1),as.double(DMACH),PACKAGE='tseriesEntropy');
        return(kdenest.mlcv$F);
    }
    res <- nlm(kdenest.mlcv, c(1.06*(min(sd(x),IQR(x)/1.34))*n^{-1/6},1.06*sd(y)*n^{-1/6}))$estimate
    return(diag(res))
}
## **************************************************************************************************
scv.uni <- function(x){
    return(hscv(x))
}
## **************************************************************************************************
scv.biv <- function(x,y,bdiag){
    if(bdiag){
        return(Hscv.diag(cbind(c(x),c(y)),pilot="samse"))
    }else{
        return(Hscv(cbind(c(x),c(y)),pilot="samse"))
    }
}
## **************************************************************************************************
pi.uni <- function(x){
    return(hpi(x))
}
## **************************************************************************************************
pi.biv <- function(x,y,bdiag){
    if(bdiag){
        return(Hpi.diag(cbind(c(x),c(y)),pilot="samse"))
    }else{
        return(Hpi(cbind(c(x),c(y)),pilot="samse"))
    }
}
## **************************************************************************************************
lscv.uni <- function(x){
    return(hlscv(x))
}
## **************************************************************************************************

lscv.biv <- function(x,y,bdiag){
    if(bdiag){
        return(Hlscv.diag(cbind(c(x),c(y))))
    }else{
        return(Hlscv(cbind(c(x),c(y))))
    }
}
## **************************************************************************************************

reference.uni <- function(x){
    #return( 1.06*(min(sd(x),IQR(x)/1.34))*length(x)^{-1/5})
    return(1.06*sd(x)*length(x)^{-1/5})
}
## **************************************************************************************************

reference.biv <- function(x,y,bdiag){
#        h1.biv <- 1.06*(min(sd(x),IQR(x)/1.34))*length(x)^{-1/6}
#        h2.biv <- 1.06*sd(y)*length(y)^{-1/6}
    h1 <- 1.06*sd(x)*length(x)^{-1/6}
    h2 <- 1.06*sd(y)*length(y)^{-1/6}
    return(diag(c(h1,h2)))
}
## **************************************************************************************************

surrogate.SA <- function(x,nlag=trunc(length(x)/4),nsurr,Te=0.0015,RT=0.9,eps.SA=0.05,nsuccmax=30,nmax=300,che=100000){

    ## Wrapper for the F90 routine SURROGATEACF
    ## Given a time series in input x computes nsurr surrogates through Simulated Annealing
    ##
    ## INPUT: ********************************************************************************
     N    <- length(x);     #  length of the series
    ##  nlag                #  minimization is performed w.r.t. to the first 'nlag' lags
    ##  Te                  #  starting value for the temperature (should be dependent on N)
    ##  RT                  #  reduction factor for the temperature
    ##  eps.SA              #  target tolerance
    ##  nsuccmax            #  Te is decreased after nsuccmax*N successes
    ##  nmax                #  Te is decreased after nmax*N iterations
    ##  nsurr               #  number of surrogates
    ##  che                 #  after check*2N global iterations the algorithm starts again
    ##
    ## OUTPUT:
    ##   A matrix with N rows and nsurr columns, in each column is stored a surrogate
    surrogate.acf <- matrix(0,N,nsurr);
    surr <- .Fortran("surrogateacf",as.double(x),as.integer(N),as.integer(nlag),as.double(Te),as.double(RT),as.double(eps.SA),
            as.integer(nsuccmax),as.integer(nmax),as.integer(nsurr),as.integer(che),surrogate.acf, PACKAGE='tseriesEntropy');
    surrogate.acf <- surr[[11]];
    return(list(surr=surrogate.acf,call=match.call()));
}

## **************************************************************************************************

Srho.test.SA <- function(x, y, lag.max = 10,  B = 100, plot = TRUE, quant = c(0.95, 0.99),
bw = c("reference", "mlcv", "lscv", "scv", "pi"), bdiag=TRUE, method =c("integral","summation"), tol=1e-03,
nlag=trunc(length(x)/4),Te=0.0015,RT=0.9,eps.SA=0.05,nsuccmax=30,nmax=300,che=100000,...)
{
    if(any(quant<=0|quant>=1)) stop("elements of quant must lie in ]0,1[");
    if(length(quant)==1){
        if(quant==0.99){
            quant <- c(0.95,quant)
        }else{
            quant <- c(quant,0.99)
        }
    }
    bw     <- match.arg(bw)
    method <- match.arg(method)
    if (missing(y)){
        S.x    <- Srho.ts(x,lag.max=lag.max,bw=bw,bdiag=bdiag,method=method, plot=FALSE, tol=tol,...)@.Data
        x.surr <- surrogate.SA(x=x,nlag=nlag,nsurr = (B+5),Te=Te,RT=RT,eps.SA=eps.SA,nsuccmax=nsuccmax,nmax=nmax,che=che)
        M      <-  safe.Srho(x.surr=x.surr$surr,B.good=B,lag.max=lag.max,bw=bw,bdiag=bdiag,method=method,tol=tol,...);
        M.95   <- apply(M,1,quantile,probs=c(quant[1]));
        M.99   <- apply(M,1,quantile,probs=c(quant[2]));
        names(S.x) <- 1:lag.max
        ind95  <- which(S.x>=M.95);
        ind99  <- which(S.x>=M.99);
        out            <- new("Srho.test")
        out@.Data      <- S.x@.Data
        out@lags       <- 1:lag.max
        out@stationary <- FALSE
        out@data.type  <- "continuous"
        out@test.type  <- "nonlinearity (SA)"
        out@call       <- match.call();
        out@call.h     <- x.surr$call;
        out@quantiles  <- cbind(M.95,M.99)
        q.names        <- paste("Q",as.character((quant*100)),"%",sep='')
        colnames(out@quantiles)     <- q.names
        rownames(out@quantiles)     <- 1:lag.max
        out@significant.lags        <- list(as.integer(names(ind95)),as.integer(names(ind99)))
        names(out@significant.lags) <- q.names
        out@p.value                 <- rowMeans(M >= S.x) # bootstrap p-value
        names(out@p.value)          <- 1:lag.max
        if (plot) {
            plot(out)
            return(invisible(out))
        }
        else return(out)

    } else {
        return(cat("Cross-Entropy testing for non linearity not yet implemented"))
    }
}

## **************************************************************************************************

surrogate.AR <- function(x, order.max=NULL, fit.method=c("yule-walker", "burg", "ols", "mle", "yw"), nsurr){
    fit.method <- match.arg(fit.method)
    n         <- length(x);
    x.ar      <- ar(x,method=fit.method,order.max = order.max, na.action=na.exclude);
    x.res.ar  <- x.ar$resid;
    ind       <- is.na(x.res.ar);
    x.res.ar[ind] <- median(x.res.ar,na.rm=T);
    x.res.ar  <- scale(x.res.ar,center=TRUE,scale=FALSE) # Recenters the residuals (NEW)
    x.par.ar  <- x.ar$ar;             # Estimated parameters AR
    le.ar     <- x.ar$order;          # Estimated order of the AR model (AIC)
    x.surr    <- matrix(0,nrow=n,ncol=nsurr) # Matrix of the bootstrap replications
    for(i in (1:nsurr)) {
        x.resb     <- sample(x.res.ar,replace=TRUE);
        x.surr[,i] <- arima.sim(n = n, list(ar = x.par.ar),innov=x.resb,n.start=le.ar+1);
    }
    return(list(surr=x.surr,call=match.call()));
}
## **************************************************************************************************

surrogate.ARs <- function(x, order.max=NULL, fit.method=c("yule-walker", "burg", "ols", "mle", "yw"), nsurr){
    ## SMOOTHED SIEVE BOOTSRAP
    fit.method <- match.arg(fit.method)
    n          <- length(x);
    x.ar       <- ar(x,method=fit.method,order.max = order.max, na.action=na.exclude);
    x.res.ar   <- x.ar$resid;
    ind        <- is.na(x.res.ar);
    x.res.ar[ind] <- median(x.res.ar,na.rm=TRUE);
    x.res.ar   <- scale(x.res.ar,center=TRUE,scale=FALSE) # Recenter the residuals
    x.par.ar   <- x.ar$ar;                 # Estimated parameters AR
    le.ar      <- x.ar$order;              # Estimated order of the AR model (AIC)
    fit        <- density(x.res.ar)        # Kernel density fit
    x.surr     <- matrix(0,nrow=n,ncol=nsurr) # Matrix of the bootstrap replications
    for(i in (1:nsurr)) {
        x.resb     <- rnorm(n, sample(x.res.ar, size = n, replace = TRUE), fit$bw) # draws from the density
        x.surr[,i] <- arima.sim(n = n, list(ar = x.par.ar),innov=x.resb,n.start=le.ar+1);
    }
    return(list(surr=x.surr,call=match.call()));
}
## **************************************************************************************************

Srho.test.AR <- function(x, y, lag.max = 10,  B = 100, plot = TRUE, quant = c(0.95, 0.99),
bw = c("reference", "mlcv", "lscv", "scv", "pi"), bdiag=TRUE, method =c("integral","summation"), tol=1e-03,
order.max=NULL, fit.method=c("yule-walker", "burg", "ols", "mle", "yw"),smoothed=TRUE,...)
{
    if(any(quant<=0|quant>=1)) stop("elements of quant must lie in ]0,1[");
    if(length(quant)==1){
        if(quant==0.99){
            quant <- c(0.95,quant)
        }else{
            quant <- c(quant,0.99)
        }
    }
    bw         <- match.arg(bw)
    method     <- match.arg(method)
    fit.method <- match.arg(fit.method)
    arg.s      <- list(x=x,order.max=order.max,fit.method=fit.method,nsurr = (B+5))
    fun        <- switch(smoothed+1,"surrogate.AR","surrogate.ARs")
    if (missing(y)){
        S.x    <- Srho.ts(x,lag.max=lag.max,bw=bw,bdiag=bdiag,method=method, plot=FALSE, tol=tol,...)@.Data
#        x.surr <- surrogate.AR(x=x,order.max=order.max,fit.method=fit.method,nsurr = (B+5))
        x.surr <-  do.call(eval(fun),args=arg.s)
        M      <-  safe.Srho(x.surr=x.surr$surr,B.good=B,lag.max=lag.max,bw=bw,bdiag=bdiag,method=method,tol=tol,...);
        M.95   <- apply(M,1,quantile,probs=c(quant[1]));
        M.99   <- apply(M,1,quantile,probs=c(quant[2]));
        names(S.x) <- 1:lag.max
        ind95  <- which(S.x>=M.95);
        ind99  <- which(S.x>=M.99);
        out            <- new("Srho.test")
        out@.Data      <- S.x@.Data
        out@lags       <- 1:lag.max
        out@stationary <- FALSE
        out@data.type  <- "continuous"
        out@test.type  <- "nonlinearity (AR)"
        out@call       <- match.call();
        out@call.h     <- x.surr$call
        out@quantiles  <- cbind(M.95,M.99)
        q.names        <- paste("Q",as.character((quant*100)),"%",sep='')
        colnames(out@quantiles)     <- q.names
        rownames(out@quantiles)     <- 1:lag.max
        out@significant.lags        <- list(as.integer(names(ind95)),as.integer(names(ind99)))
        names(out@significant.lags) <- q.names
        out@p.value                 <- rowMeans(M >= S.x) # bootstrap p-value
        names(out@p.value)          <- 1:lag.max
        if (plot) {
            plot(out)
            return(invisible(out))
        }
        else return(out)

    } else {
        return(cat("Cross-Entropy testing for non linearity not yet implemented"))
    }
}
## **************************************************************************************************

    safe.Srho <- function(x.surr,y.surr,B.good,lag.max,writeout=200,bw,bdiag,method,tol,...){

        ## Function to protect against crashes in computing Srho.ts on surrogates or bootstrap replicates

        ## INPUT:
        ##
        ##  [...]     same as Srho.ts
        ##  B.good       :  minimum number of replications requested.
        ##  x.surr, ysurr:  matrix of the surrogates, one column each.
        ##                  It has to be ncol(x.surr) > B.good.
        ##  writeout     :  writes every writeout replications.
        ## OUTPUT:
        ##  S.surr       :  nlag by B.good matrix with the results.

        n.surr <- ncol(x.surr);
        if(missing(y.surr)||is.null(y.surr)){
            y.surr <- NULL
            S.surr <- matrix(0,nrow=lag.max,ncol=B.good)
        }else{
            if(any(dim(x.surr)!=dim(y.surr))){stop("x.surr and y.surr do not match")}
            S.surr <- matrix(0,nrow=(2*lag.max+1),ncol=B.good)
        }
        if(n.surr<B.good) stop("The number of surrogates generated are less than the number of the results requested")
        if(n.surr==B.good) warning("The number of surrogates generated are equal to the number of the results requested")
        k <- 1;
        j <- 1;
       # cat("\r Replication number  1");
        while(k <= B.good){
            result <- try(Srho.ts(x=x.surr[,j],y=y.surr[,j],lag.max=lag.max,bw=bw,bdiag=bdiag,method=method,plot=FALSE,tol=tol,...)@.Data,TRUE)
            cond <- (class(result)=="try-error"|any(is.na(result))|any(is.nan(result)))
            if(cond){
                cat("\r ***************************************");
                cat("\r Crash encountered, surrogate skipped \n");
                cat("\r Surrogate number ", j);
                cat("\r ***************************************");
            } else {
                S.surr[,k] <- result
                if(k%%writeout==0) cat("\r Replication number ", k); ## Writes every writeout replications
                k <- k+1;
            }
            j <- j+1;
        }
        return(S.surr)
    }

## **************************************************************************************************

Srho.test.ts <- function(x, y, lag.max = 10,  B = 100, plot = TRUE, quant = c(0.95, 0.99),
bw = c("reference", "mlcv", "lscv", "scv", "pi"), bdiag=TRUE, method =c("integral","summation"), tol=1e-03,
ci.type = c("mbb","perm"),...)
{
    if(any(quant<=0|quant>=1)) stop("elements of quant must lie in ]0,1[");
    if(length(quant)==1){
        if(quant==0.99){
            quant <- c(0.95,quant)
        }else{
            quant <- c(quant,0.99)
        }
    }
    bw      <- match.arg(bw)
    method  <- match.arg(method)
    ci.type <- match.arg(ci.type)
    if(missing(y)||is.null(y)){
        y <- NULL
        y.surr <- NULL
        lag.names <- 1:lag.max
        ci.type <- "perm"
    }else{
      if(!exists('blag')){
        blag <- lag.max
      }
        X <- ts.intersect(as.ts(x), as.ts(y)) # time alignment
        x <- X[,1]
        y <- X[,2]
        y.surr  <- switch(ci.type,"mbb"=mbboot(x=y,B=B+5,l=blag),"perm"=boot.perm(x=y, B=B+5))
        lag.names <- -lag.max:lag.max
    }
    x.surr <- switch(ci.type,"mbb"=mbboot(x=x,B=B+5,l=blag),"perm"=boot.perm(x=x, B=B+5))
    S.x    <- Srho.ts(x=x,y=y,lag.max=lag.max,bw=bw,bdiag=bdiag,method=method, plot=FALSE,tol=tol,...)@.Data
    M      <-  safe.Srho(x.surr=x.surr$surr,y.surr=y.surr$surr,B.good=B,lag.max=lag.max,
                        bw=bw,bdiag=bdiag,method=method,tol=tol,...);
    M.95   <- apply(M,1,quantile,probs=c(quant[1]));
    M.99   <- apply(M,1,quantile,probs=c(quant[2]));
    names(S.x) <- lag.names
    ind95  <- which(S.x>=M.95);
    ind99  <- which(S.x>=M.99);
    out <- new("Srho.test")
    out@.Data      <- S.x@.Data
    out@lags       <- lag.names
    out@stationary <- FALSE
    out@data.type  <- "continuous"
    out@test.type  <- paste("independence",ci.type,sep="-")
    out@call       <- match.call();
    out@quantiles  <- cbind(M.95,M.99)
    q.names <- paste("Q",as.character((quant*100)),"%",sep='')
    colnames(out@quantiles)     <- q.names
    rownames(out@quantiles)     <- lag.names
    out@significant.lags        <- list(as.integer(names(ind95)),as.integer(names(ind99)))
    names(out@significant.lags) <- q.names
    out@p.value                 <- rowMeans(M >= S.x) # bootstrap p-value
    names(out@p.value)          <- lag.names
    if (plot) {
        plot(out)
        return(invisible(out))
    }
    else return(out)
}
## **************************************************************************************************

boot.perm <- function(x, B){
## generates B random permutations of x
# INPUT:  vector x of length n
# OUTPUT: n by B matrix; each column contains a random permutation of x
    n      <- length(x)
    x.surr <- matrix(rep(x,B), n,B)
    x.surr <- apply(x.surr,FUN=sample,MARGIN=2)
    return(list(surr=x.surr,call=match.call()));
}
## **************************************************************************************************

mbboot <- function(x, B, l){
# Moving Block Bootstrap
# INPUT  x: vector of length n
#        l: length of the MBB block
# OUTPUT x.surr : n by B matrix; each column contains a MBB resample from x
    n <- length(x)
    if ((l < 1) || (l >= n)) stop("l should be in [1,length(x)]")
    nblocks <- n%/%l+1 # max. number of blocks for each resample
    ind.mat   <- embed(1:n,l)[,l:1] # embedding of the indices
    ind.block <- sample.int(n=n-l+1,size=nblocks*B,replace=TRUE) # block indices
    ind.x     <- as.vector(t(ind.mat[ind.block,]))
    x.surr <- matrix(x[ind.x],nrow=nblocks*l,ncol=B)[1:n,]
 #   return(x.surr)
    return(list(surr=x.surr,call=match.call()));
}
## **************************************************************************************************

    safe.Trho <- function(x.surr,y.surr,B.good,lag.max,writeout=200,bw,bdiag,method,tol,...){

        ## Function to protect against crashes in computing Trho on surrogates or bootstrap replicates

        ## INPUT:
        ##
        ##  [...]     same as Srho.ts
        ##  B.good       :  minimum number of replications requested.
        ##  x.surr, ysurr:  matrix of the surrogates, one column each.
        ##                  It has to be ncol(x.surr) > B.good.
        ##  writeout     :  writes every writeout replications.
        ## OUTPUT:
        ##  S.surr       :  nlag by B.good matrix with the results.

        n.surr <- ncol(x.surr);
        if(missing(y.surr)||is.null(y.surr)){
            y.surr <- NULL
            S.surr <- matrix(0,nrow=lag.max,ncol=B.good)
        }else{
            if(any(dim(x.surr)!=dim(y.surr))){stop("x.surr and y.surr do not match")}
            S.surr <- matrix(0,nrow=(2*lag.max+1),ncol=B.good)
        }
        if(n.surr<B.good) stop("The number of surrogates generated are less than the number of the results requested")
        if(n.surr==B.good) warning("The number of surrogates generated are equal to the number of the results requested")
        k <- 1;
        j <- 1;
       # cat("\r Replication number  1");
        while(k <= B.good){
            result <- try((Srho.ts(x=x.surr[,j],y=y.surr[,j],lag.max=lag.max,bw=bw,bdiag=bdiag,method=method,plot=FALSE,tol=tol,...)@.Data-
             Srho.cor(x=x.surr[,j],lag.max=lag.max,plot=FALSE)@.Data)^2,TRUE)
            cond <- (class(result)=="try-error"|any(is.na(result))|any(is.nan(result)))
            if(cond){
                cat("\r ***************************************");
                cat("\r Crash encountered, surrogate skipped \n");
                cat("\r Surrogate number ", j);
                cat("\r ***************************************");
            } else {
                S.surr[,k] <- result
                if(k%%writeout==0) cat("\r Replication number ", k); ## Writes every writeout replications
                k <- k+1;
            }
            j <- j+1;
        }
        return(S.surr)
    }

## **************************************************************************************************

Trho.test.AR <- function(x, y, lag.max = 10,  B = 100, plot = TRUE, quant = c(0.95, 0.99),
bw = c("reference", "mlcv", "lscv", "scv", "pi"), bdiag=TRUE, method =c("integral","summation"), tol=1e-03,
order.max=NULL, fit.method=c("yule-walker", "burg", "ols", "mle", "yw"),smoothed=TRUE,...)
{
    if(any(quant<=0|quant>=1)) stop("elements of quant must lie in ]0,1[");
    if(length(quant)==1){
        if(quant==0.99){
            quant <- c(0.95,quant)
        }else{
            quant <- c(quant,0.99)
        }
    }
    bw         <- match.arg(bw)
    method     <- match.arg(method)
    fit.method <- match.arg(fit.method)
    arg.s      <- list(x=x,order.max=order.max,fit.method=fit.method,nsurr = (B+5))
    fun        <- switch(smoothed+1,"surrogate.AR","surrogate.ARs")
    if (missing(y)){
        S.x    <- (Srho.ts(x,lag.max=lag.max,bw=bw,bdiag=bdiag,method=method, plot=FALSE,
                                tol=tol,...)@.Data - Srho.cor(x,lag.max=lag.max,plot=FALSE)@.Data)^2
        x.surr <-  do.call(eval(fun),args=arg.s)
#        x.surr <- surrogate.AR(x=x,order.max=order.max,fit.method=fit.method,nsurr = (B+5))

        M      <-  safe.Trho(x.surr=x.surr$surr,B.good=B,lag.max=lag.max, bw=bw,bdiag=bdiag,method=method,tol=tol,...);
        M.95   <- apply(M,1,quantile,probs=c(quant[1]));
        M.99   <- apply(M,1,quantile,probs=c(quant[2]));
        names(S.x) <- 1:lag.max
        ind95      <- which(S.x>=M.95);
        ind99      <- which(S.x>=M.99);
        out        <- new("Srho.test")
        out@.Data  <- S.x
        out@lags   <- 1:lag.max
        out@stationary <- FALSE
        out@data.type  <- "continuous"
        out@test.type  <- "nonlinearity (AR) Nonparametric vs. Parametric"
        out@call       <- match.call();
        out@call.h     <- x.surr$call
        out@quantiles  <- cbind(M.95,M.99)
        q.names        <- paste("Q",as.character((quant*100)),"%",sep='')
        colnames(out@quantiles) <- q.names
        rownames(out@quantiles) <- 1:lag.max
        out@significant.lags    <- list(as.integer(names(ind95)),as.integer(names(ind99)))
        names(out@significant.lags) <- q.names
        out@p.value                 <- rowMeans(M >= S.x) # bootstrap p-value
        names(out@p.value)          <- 1:lag.max
        if (plot) {
            plot(out,ylab = "T")
            return(invisible(out))
        }
        else return(out)

    } else {
        return(cat("Cross-Entropy testing for non linearity not yet implemented"))
    }
}
## **************************************************************************************************

Trho.test.SA <- function(x, y, lag.max = 10,  B = 100, plot = TRUE, quant = c(0.95, 0.99),
bw = c("reference", "mlcv", "lscv", "scv", "pi"), bdiag=TRUE, method =c("integral","summation"), tol=1e-03,
nlag=trunc(length(x)/4),Te=0.0015,RT=0.9,eps.SA=0.05,nsuccmax=30,nmax=300,che=100000,...)
{
    if(any(quant<=0|quant>=1)) stop("elements of quant must lie in ]0,1[");
    if(length(quant)==1){
        if(quant==0.99){
            quant <- c(0.95,quant)
        }else{
            quant <- c(quant,0.99)
        }
    }
    bw         <- match.arg(bw)
    method     <- match.arg(method)
    if (missing(y)){
        S.x    <- (Srho.ts(x,lag.max=lag.max,bw=bw,bdiag=bdiag,method=method, plot=FALSE,
                                tol=tol,...)@.Data - Srho.cor(x,lag.max=lag.max,plot=FALSE)@.Data)^2
        x.surr <- surrogate.SA(x=x,nlag=nlag,nsurr = (B+5),Te=Te,RT=RT,eps.SA=eps.SA,nsuccmax=nsuccmax,nmax=nmax,che=che)
        M      <-  safe.Trho(x.surr=x.surr$surr,B.good=B,lag.max=lag.max, bw=bw,bdiag=bdiag,method=method,tol=tol,...);
        M.95  <- apply(M,1,quantile,probs=c(quant[1]));
        M.99  <- apply(M,1,quantile,probs=c(quant[2]));
        names(S.x) <- 1:lag.max
        ind95 <- which(S.x>=M.95);
        ind99 <- which(S.x>=M.99);
        out            <- new("Srho.test")
        out@.Data      <- S.x
        out@lags       <- 1:lag.max
        out@stationary <- FALSE
        out@data.type  <- "continuous"
        out@test.type  <- "nonlinearity (SA)"
        out@call       <- match.call();
        out@call.h     <- x.surr$call;
        out@quantiles  <- cbind(M.95,M.99)
        q.names <- paste("Q",as.character((quant*100)),"%",sep='')
        colnames(out@quantiles) <- q.names
        rownames(out@quantiles) <- 1:lag.max
        out@significant.lags    <- list(as.integer(names(ind95)),as.integer(names(ind99)))
        names(out@significant.lags) <- q.names
        out@p.value                 <- rowMeans(M >= S.x) # bootstrap p-value
        names(out@p.value)          <- 1:lag.max
        if (plot) {
            plot(out,ylab = "T")
            return(invisible(out))
        }
        else return(out)

    } else {
        return(cat("Cross-Entropy testing for non linearity not yet implemented"))
    }
}
## **************************************************************************************************

Srho.cor <- function(x,lag.max=10,plot=TRUE){
    # Parametric estimation of Srho based on the ACF
    x.cor <- acf(x,lag.max=lag.max,plot=FALSE,na.action=na.pass)$acf[2:(lag.max+1),,1]
    out <- new("Srho")
    out@.Data      <- cor2Srho(x.cor)
    out@lags       <- 1:lag.max
    out@data.type  <- "continuous"
    out@stationary <- TRUE
    if (plot) {
        plot(out)
        return(invisible(out))
    }
    else return(out)
}

## ****************************************************************************************************

cor2Srho <- function(rho){
    S <- 1- (2*(1-rho^2)^(1/4))/(4-rho^2)^(1/2);
    return(as.vector(S))
}
## ****************************************************************************************************
