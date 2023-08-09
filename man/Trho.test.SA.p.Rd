\name{Trho.test.SA.p}
\alias{Trho.test.SA.p}
\alias{Trho.test.SA}

\title{Entropy Tests For Nonlinearity In Time Series - Parallel Version}

\description{Entropy test of nonlinearity for time series based
on \code{\link{Srho.ts}} and surrogate data obtained through
Simulated Annealing. The parallel version requires \code{parallel}.}
 \usage{
Trho.test.SA(x, y, lag.max = 10,  B = 100, plot = TRUE, quant = c(0.95, 0.99),
 bw = c("reference","mlcv", "lscv", "scv", "pi"), bdiag=TRUE,
 method =c("integral","summation"), tol=1e-03, nlag=trunc(length(x)/4),
 Te=0.0015, RT=0.9, eps.SA=0.05, nsuccmax=30, nmax=300, che=100000,...) 
 
Trho.test.SA.p(x, y, lag.max = 10,  B = 100, plot = TRUE, quant = c(0.95, 0.99),
 bw = c("reference","mlcv", "lscv", "scv", "pi"), bdiag=TRUE,
  method =c("integral","summation"), tol=1e-03, nlag=trunc(length(x)/4), Te=0.0015,
 RT=0.9, eps.SA=0.05, nsuccmax=30, nmax=300, che=100000, nwork=detectCores(),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{univariate numeric time series object or numeric vectors (\code{y} is missing in the univariate case).}
  \item{lag.max}{maximum lag at which to calculate Trho; the default is 10.}
  \item{B}{number of surrogate time series.}
  \item{plot}{logical. If \code{TRUE}(the default) produces a plot of Trho together with confidence bands under the null hypothesis of linearity at 95\% and 99\%.}
   \item{quant}{quantiles to be specified for the computation of the significant lags and the plot of confidence bands. Up to 2 quantiles can be specified.
  Defaults are 95\% and 99\%.}
  \item{bw}{see \code{\link{Srho.ts}}.}
  \item{bdiag}{see \code{\link{Srho.ts}}.}
  \item{method}{see \code{\link{Srho.ts}}.}
  \item{tol}{see \code{\link{Srho.ts}}.}
  \item{nlag}{see \code{\link{surrogate.SA}}.}
  \item{Te}{see \code{\link{surrogate.SA}}.}
  \item{RT}{see \code{\link{surrogate.SA}}.}
  \item{eps.SA}{see \code{\link{surrogate.SA}}.}
  \item{nsuccmax}{see \code{\link{surrogate.SA}}.}
  \item{nmax}{see \code{\link{surrogate.SA}}.}
  \item{che}{see \code{\link{surrogate.SA}}.}
  \item{nwork}{number of workers/processes to be used in parallel environments.}
  \item{...}{further arguments, typically passed to \code{\link[cubature]{hcubature}}.}
}
 \details{
For each lag from 1 to \code{lag.max} \code{Trho.test.SA} computes a test for nonlinearity for time series based
on \code{\link{Srho.ts}}. The distribution under the null hypothesis of a linear Gaussian process is obtained through a generalization of surrogate data
methods. Surrogate time series are obtained through Simulated Annealing (SA). Sensible (N-dependent) defaults are derived for the parameters of the SA
algorithm, there should not be the need to change them. The routine requires the package parallel to spawn multiple workers.
}
\value{An object of class "Srho.test", which is a list with the following elements:
    \item{.Data}{ vector of \code{lag.max} elements containing Trho computed at each lag.}
    \item{\code{call}:}{Object of class \code{"call"}: contains the call to the routine.}
    \item{\code{call.h}:}{Object of class \code{"call"}: contains the call to the routine used for obtaining the surrogates or the bootstrap replicates under the null hypothesis.}
    \item{quantiles}{Object of class \code{"matrix"}: contains the quantiles of the surrogate distribution under the null hypothesis.}
    \item{\code{test.type}}{Object of class \code{"character"}: contains a description of the type of test performed.}
    \item{significant.lags}{Object of class \code{"list"}: contains the lags at which Trho exceeds the confidence bands at \code{quant}\% under the null hypothesis.}
    \item{p.value}{Object of class \code{"numeric"}: contains the bootstrap p-value for each lag.}
    \item{lags}{integer vector that contains the lags at which Trho is computed.}
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
\seealso{See Also \code{\link{Srho.ts}}, \code{\link{surrogate.SA}}, \code{\link{Trho.test.SA}}.}
\examples{
\dontrun{
# modify nwork to match the number of available cores
set.seed(13)
x      <- arima.sim(n=120, model = list(ar=0.8));
result <- Trho.test.SA.p(x, lag.max = 5, B = 100, bw='reference', nwork=2)

## ** Compare timings **
system.time(Trho.test.SA.p(x, lag.max = 5, B = 100, bw='reference', nwork=4))
system.time(Trho.test.SA(x, lag.max = 5, B = 100, bw='reference'))
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
