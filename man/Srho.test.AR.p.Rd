\name{Srho.test.AR.p}
\alias{Srho.test.AR.p}
\alias{Srho.test.AR}
\title{Entropy Tests For Nonlinearity In Time Series - Parallel Version}

\description{Entropy test of nonlinearity for time series based
on \code{\link{Srho.ts}} and surrogate data obtained through
the sieve bootstrap. The parallel version requires \code{parallel}.}

\usage{
Srho.test.AR(x, y, lag.max = 10, B = 100, plot = TRUE, quant = c(0.95, 0.99),
 bw = c("reference", "mlcv", "lscv", "scv", "pi"), bdiag=TRUE,
 method = c("integral", "summation"), tol = 0.001, order.max = NULL,
 fit.method=c("yule-walker", "burg", "ols", "mle", "yw"),  smoothed = TRUE ,...)

## Parallel version

Srho.test.AR.p(x, y, lag.max = 10, B = 100, plot = TRUE, quant = c(0.95, 0.99),
 bw = c("reference", "mlcv", "lscv", "scv", "pi"), bdiag=TRUE,
 method = c("integral", "summation"), tol = 0.001, order.max = NULL,
 fit.method=c("yule-walker", "burg", "ols", "mle", "yw"),  smoothed = TRUE,
 nwork=detectCores(),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{univariate numeric time series object or numeric vectors (\code{y} is missing in the univariate case).}
  \item{lag.max}{maximum lag at which to calculate Srho; the default is 10.}
  \item{B}{number of surrogate time series.}
  \item{plot}{logical. If \code{TRUE} (the default) produces a plot of Srho together with confidence bands under the null hypothesis of linearity at 95\% and 99\%.}
  \item{quant}{quantiles to be specified for the computation of the significant lags and the plot of confidence bands. Up to 2 quantiles can be specified.
  Defaults are 95\% and 99\%.}
  \item{bw}{see \code{\link{Srho.ts}}.}
  \item{bdiag}{see \code{\link{Srho.ts}}.}
  \item{method}{see \code{\link{Srho.ts}}.}
  \item{tol}{see \code{\link{Srho.ts}}.}
  \item{order.max}{see \code{\link{surrogate.ARs}}.}
  \item{fit.method}{see \code{\link{surrogate.ARs}}.}
  \item{smoothed}{logical. If \code{TRUE} (the default) uses the smoothed sieve bootstrap in \code{\link{surrogate.ARs}} to generate surrogates. Otherwise uses the classic sieve by calling \code{\link{surrogate.AR}}.}
  \item{nwork}{number of workers/processes to be used in parallel environments.}
  \item{...}{further arguments, typically passed to \code{\link[cubature]{hcubature}}.}
}

\details{
For each lag from 1 to \code{lag.max} \code{Srho.test.AR} computes a test for nonlinearity for time series based
on \code{\link{Srho.ts}}. The distribution under the null hypothesis of linearity is obtained through the sieve bootstrap.
The routine requires the package parallel to spawn multiple workers.
}

\value{An object of class "Srho.test", which is a list with the following elements:
    \item{.Data}{ vector of \code{lag.max} elements containing Srho computed at each lag.}
    \item{\code{call}:}{Object of class \code{"call"}: contains the call to the routine.}
    \item{\code{call.h}:}{Object of class \code{"call"}: contains the call to the routine used for obtaining the surrogates or the bootstrap replicates under the null hypothesis}
    \item{quantiles}{Object of class \code{"matrix"}: contains the quantiles of the surrogate distribution under the null hypothesis.}
    \item{\code{test.type}}{Object of class \code{"character"}: contains a description of the type of test performed.}
    \item{significant.lags}{Object of class \code{"list"}: contains the lags at which Srho exceeds the confidence bands at \code{quant}\% under the null hypothesis.}
    \item{p.value}{Object of class \code{"numeric"}: contains the bootstrap p-value for each lag.}
    \item{lags}{integer vector that contains the lags at which Srho is computed.}
    \item{stationary}{Object of class \code{"logical"}:  \code{TRUE} if the stationary version is computed. Set to \code{FALSE} by default as only the non-stationary version is implemented.}
    \item{data.type}{Object of class \code{"character"}: contains the data type.}
    \item{notes}{Object of class \code{"character"}: additional notes.}
}
\references{
Giannerini S., Maasoumi E., Bee Dagum E., (2015), Entropy testing
for nonlinear serial dependence in time series, \emph{Biometrika}, \bold{102(3)}, 661--675
\doi{10.1093/biomet/asv007}.
}
 \author{Simone Giannerini<simone.giannerini@unibo.it>}
\seealso{See Also \code{\link{Srho.ts}}, \code{\link{surrogate.AR}}, \code{\link{surrogate.ARs}}, \code{\link{Srho.test.AR}}.}
\examples{
\dontrun{
## ************************************************************
## WARNING: computationally intensive, increase B with caution
## ************************************************************

# modify nwork to match the number of available cores
set.seed(13)
x      <- arima.sim(n=120, model = list(ar=0.8));
result <- Srho.test.AR.p(x, lag.max = 5,  B = 100, bw='reference', method='integral', nwork=2)

## ** Compare timings **
system.time(Srho.test.AR.p(x, lag.max = 5,  B = 100, bw='reference', method='integral', nwork=4))
system.time(Srho.test.AR(x, lag.max = 5,  B = 100, bw='reference', method='integral'))
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
