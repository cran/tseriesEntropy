\name{Srho}
\alias{Srho}
\title{Entropy Measure Of Serial And Cross Dependence}
\description{
Entropy based measure of serial and cross dependence for integer or categorical data. Implements a normalized version
of the Hellinger/Matusita distance. As shown in the references the metric measure is a proper distance.}

\usage{
Srho(x, y, lag.max, stationary = TRUE, plot = TRUE, version = c("FORTRAN","R"),
 nor = FALSE)
}

\arguments{
  \item{x, y}{integer or factor time series objects or vectors. (\code{y} is missing in the univariate case).}
  \item{lag.max}{maximum lag at which to calculate Srho; default is \code{round(N/4)} where N is the number of observations.}
  \item{stationary}{logical. If \code{TRUE} assumes stationarity and computes marginal probabilities by using N observations. If \code{FALSE} uses N-k observations where k is the lag.}
  \item{plot}{logical. If \code{TRUE} (the default) Srho is plotted.}
  \item{version}{either \code{"FORTRAN"} or \code{"R"}. \code{FORTRAN} version is the default and is preferred over the pure \code{R} version which is considerably slower but is included in case of portability issues.}
  \item{nor}{logical. If \code{TRUE} normalizes Srho with respect to its attainable maximum. Defaults to \code{FALSE}.}
}

\details{
\describe{
\item{Univariate version: serial entropy}{\preformatted{Srho(x, lag.max,
stationary = TRUE, plot = TRUE, version = c("FORTRAN","R"), nor = FALSE)}}

\item{Bivariate version: cross entropy}{\preformatted{Srho(x, y, lag.max,
stationary = TRUE, plot = TRUE, version = c("FORTRAN","R"), nor = FALSE)}}}

This implementation of the measure is normalized to take values in [0, 1]. Normalization is
performed with respect to the maximum attainable value computed analytically.
This makes the results of Srho comparable among different series.}

\value{An object of S4 class "Srho", which is a list with the following elements:
    \item{.Data}{ vector of \code{lag.max} elements containing Srho computed at each lag.}
    \item{lags}{integer vector that contains the lags at which Srho is computed.}
    \item{stationary}{Object of class \code{"logical"}:  \code{TRUE} if the stationary version is computed.}
    \item{data.type}{Object of class \code{"character"}: contains the data type.}
    \item{notes}{Object of class \code{"character"}: additional notes.}
}
\references{
Granger C. W. J., Maasoumi E., Racine J., (2004) A dependence metric for possibly nonlinear processes.
\emph{Journal of Time Series Analysis}, \bold{25(5)}, 649--669.

Giannerini S., Maasoumi E., Bee Dagum E., (2015), Entropy testing
for nonlinear serial dependence in time series, \emph{Biometrika}, \bold{102(3)}, 661--675
\doi{10.1093/biomet/asv007}.

Maasoumi E., (1993) A compendium to information theory in economics and econometrics.
\emph{Econometric Reviews}, \bold{12(2)}, 137--181.

}
 \author{Simone Giannerini<simone.giannerini@unibo.it>}

\section{Warning}{
Unlike \code{ccf} the lag k value returned by \code{Srho(x,y)} estimates Srho between \code{x[t]} and \code{y[t+k]}.
The result is returned invisibly if plot is TRUE.
}

\seealso{See Also \code{\link{Srho.test}}. The function \code{\link{Srho.ts}} implements the same measure for numeric data.}
\examples{

## UNIVARIATE VERSION
x <- as.integer(rbinom(n=20,size=4,prob=0.5))
Srho(x,lag.max=4)

## BIVARIATE VERSION
y <- as.integer(rbinom(n=20,size=4,prob=0.5))
Srho(x,y,lag.max=4)

## EXAMPLE  1: the effect of normalization
## computes the maximum attainable value by correlating x with itself

set.seed(12)
K    <- 5           # number of categories
smax <- 1-1/sqrt(K) # theoretical maximum under the uniform distribution
x    <- as.integer(sample(1:K,size=1e3,replace=TRUE)) # generates the sequence
S    <- Srho(x,x,lag.max=2,nor=FALSE,plot=FALSE)

plot(S,lwd=2,col=4)
abline(h=smax,col=2,lty=2)
text(x=-1,y=0.54,labels=paste("theoretical maximum = ",round(smax,4),sep=""),col=2)
text(x=-1,y=0.45,labels=paste("estimated maximum = ",round(S[3],4),sep=""),col=4)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
