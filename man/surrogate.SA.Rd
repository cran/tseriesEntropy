\name{surrogate.SA}
\alias{surrogate.SA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Surrogate Time Series Through Simulated Annealing}
\description{
Starting from a time series \code{x} given as input, the function generates surrogate series through Simulated Annealing.
Each surrogate series is a constrained random permutation having the same autocorrelation function (up to \code{nlag} lags) of the original series \code{x}.
The surrogates can be used for testing for non linearity in time series.
}
\usage{
surrogate.SA(x, nlag, nsurr, Te = 0.0015, RT = 0.9, eps.SA = 0.05, nsuccmax = 30,
 nmax = 300, che = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a univariate numeric time series object or a numeric vector.}
  \item{nlag}{minimization is performed w.r.t. to the first \code{nlag} lags.}
  \item{nsurr}{number of surrogates.}
  \item{Te}{starting value for the temperature.}
  \item{RT}{reduction factor for the temperature \code{Te}.}
  \item{eps.SA}{target tolerance.}
  \item{nsuccmax}{Te is decreased after \code{nsuccmax*N} successes.}
  \item{nmax}{Te is decreased after \code{nmax*N} successes.}
  \item{che}{after \code{che*2N} global iterations the algorithm starts again.}
}
\details{
Let \code{N} be the length of the series \code{x}. Sensible (\code{N}-dependent) defaults are derived for the parameters of the
algorithm, there should not be the need to change them. In case, the user could try increasing \code{eps.SA}.
}
\value{
A list with the following elements:
    \item{surr}{a matrix with \code{N} rows and \code{nsurr} columns, in each column is stored a surrogate.}
    \item{call}{contains the call to the routine.}
}
\references{
Giannerini S., Maasoumi E., Bee Dagum E., (2015), Entropy testing
for nonlinear serial dependence in time series, \emph{Biometrika}, \bold{102(3)}, 661--675
\doi{10.1093/biomet/asv007}.

Schreiber T., Schmitz A.,(2000) Surrogate time series. \emph{Physica D}, \bold{142(3-4)}, 346--382.
}

\author{Simone Giannerini<simone.giannerini@unibo.it>}

\seealso{See Also \code{\link{Trho.test.SA}}, \code{\link{surrogate.AR}}, \code{\link{Trho.test.AR}}.}
\examples{
set.seed(1345)
# Generates a AR(1) series
x      <- arima.sim(n=120, model = list(ar=0.8));
x.surr <- surrogate.SA(x, nlag=10, nsurr=3);
plot.ts(x.surr$surr,col=4);


## Check that the surrogates have the same ACF of x
corig <- acf(x,10,plot=FALSE)$acf[,,1];
csurr <- acf(x.surr$surr[,1],10,plot=FALSE)$acf[,,1];
round(cbind(corig,csurr,"abs(difference)"=abs(corig-csurr)),3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
