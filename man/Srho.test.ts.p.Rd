\name{Srho.test.ts.p}
\alias{Srho.test.ts.p}
\title{Entropy Tests Of Serial And Cross Dependence For Time Series - Parallel Version}

\description{Entropy test of serial and cross dependence for numeric time series (continuous state space) based on
\code{\link{Srho.ts}}. The distribution under the null hypothesis of independence is obtained by means of bootstrap/permutations methods
(see \code{ci.type}). Parallel version -- requires parallel.}

 \usage{
Srho.test.ts.p(x, y, lag.max = 10,  B = 100, plot = TRUE, quant = c(0.95, 0.99),
 bw = c("reference","mlcv", "lscv", "scv", "pi"), method =c("integral","summation"),
 maxpts=0, tol=1e-03, ci.type = c("mbb","perm"), nslaves=detectCores())
}
%- maybe also 'usage' for other objects documented here.
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{univariate numeric time series object or numeric vectors (\code{y} is missing in the univariate case).}
  \item{lag.max}{maximum lag at which to calculate Srho; default is \code{trunc(N/4)} where N is the number of observations.}
  \item{B}{number of bootstrap/permutation replications.}
  \item{plot}{logical. If \code{TRUE}(the default) produces a plot of Srho together with confidence bands under the null hypothesis at
  levels set by \code{quant}.}
  \item{quant}{quantiles to be specified for the computation of the significant lags and the plot of confidence bands. Up to 2 quantiles can be specified.
  Defaults are 95\% and 99\%.}
  \item{bw}{see \code{\link{Srho.ts}}.}
  \item{method}{see \code{\link{Srho.ts}}.}
  \item{maxpts}{see \code{\link{Srho.ts}}.}
  \item{tol}{see \code{\link{Srho.ts}}.}
  \item{ci.type}{confidence interval type. determines how the distribution under the null hypothesis is obtained. \code{mbb} uses a moving block bootstrap
  with block length equal to \code{lag.max}. The option "perm" uses permutation methods (each resampled series is a random permutation of the original series).
  The option \code{mbb} makes sense only in the bivariate case for which is the default.}
  \item{nslaves}{number of slaves/processes to be used in parallel environments.}
}

 \details{
\describe{
\item{Univariate version: test for serial dependence}{\preformatted{Srho.test.ts(x, lag.max = 10,
B = 100, plot = TRUE, quant = c(0.95, 0.99), 
bw = c("reference", "mlcv", "lscv", "scv", "pi"), method =c("integral","summation"), 
maxpts=0, tol=1e-03, ci.type = c("perm"), nslaves=detectCores())}}

\item{Bivariate version: test for cross dependence}{\preformatted{Srho.test.ts(x, y, lag.max = 10,
B = 100, plot = TRUE, quant = c(0.95, 0.99), 
bw = c("reference", "mlcv", "lscv", "scv", "pi"), method =c("integral","summation"), 
maxpts=0, tol=1e-03, ci.type = c("mbb","perm"), nslaves=detectCores())}}
}
For each lag from 1 to \code{lag.max} (serial dependence) or from \code{-lag.max} to \code{lag.max} (cross dependence) \code{Srho.test.ts} computes
 a test for serial/cross dependence for time series based on \code{\link{Srho.ts}}. The distribution under the null hypothesis of independence
 is obtained through either permutation or bootstrap methods. If the option \code{mbb} is chosen (bivariate case only)  the resampled series
 keep the (serial) dependence structure of the original series so that the test will have better size than the permutation version.
}

\value{An object of class "Srho.test", which is a list with the following elements:
    \item{.Data}{ vector containing Srho computed at each lag.}
    \item{\code{call}:}{Object of class \code{"call"}: contains the call to the routine.}
    \item{\code{call.h}:}{Object of class \code{"call"}: contains the call to the routine used for obtaining the surrogates or the bootstrap replicates under the null hypothesis.}
    \item{quantiles}{Object of class \code{"matrix"}: contains the quantiles of the distribution under the null hypothesis.}
    \item{\code{test.type}}{Object of class \code{"character"}: contains a description of the type of test performed.}
    \item{significant.lags}{Object of class \code{"list"}: contains the lags at which Srho exceeds the confidence bands at \code{quant}\% under the null hypothesis.}
    \item{p.value}{Object of class \code{"numeric"}: contains the bootstrap p-value for each lag.}
    \item{lags}{integer vector that contains the lags at which Srho is computed.}
    \item{stationary}{Object of class \code{"logical"}:  \code{TRUE} if the stationary version is computed. Set to \code{FALSE} by default as only the non-stationary version is implemented.}
    \item{data.type}{Object of class \code{"character"}: contains the data type.}
    \item{notes}{Object of class \code{"character"}: additional notes.}
}
\references{
Granger C. W. J., Maasoumi E., Racine J., (2004) A dependence metric for possibly nonlinear processes.
\emph{Journal of Time Series Analysis}, \bold{25(5)}, 649--669.

Maasoumi E., (1993) A compendium to information theory in economics and econometrics.
\emph{Econometric Reviews}, \bold{12(2)}, 137--181.
}
\author{Simone Giannerini<simone.giannerini@unibo.it>}
\seealso{See Also \code{\link{Srho.test.ts}} and \code{\link{Srho.ts}}. 
The function \code{\link{Srho.test}} implements the same test for integer/categorical data. 
For a test for nonlinear serial dependence see \code{\link{Srho.test.AR}}, 
\code{\link{Trho.test.AR}}, \code{\link{Trho.test.SA}}, together with their parallel 
versions: \code{\link{Srho.test.AR.p}}, \code{\link{Trho.test.AR}}, \code{\link{Trho.test.SA}}.}

\examples{
\dontrun{
## ************************************************************
## WARNING: computationally intensive, increase B with caution
## ************************************************************
set.seed(13)
n      <- 120
w      <- rnorm(n)
x      <- arima.sim(n, model = list(ar=0.8));
y      <- arima.sim(n, model = list(ar=0.8));
z      <- lag(x,-1) + rnorm(n,sd=2) # dependence at lag 1
# UNIVARIATE VERSION
res1 <- Srho.test.ts.p(w, lag.max = 5,  B = 40, ci.type="perm") # independence
res2 <- Srho.test.ts.p(x, lag.max = 5,  B = 40, ci.type="perm") # dependence

# BIVARIATE VERSION
res3 <- Srho.test.ts.p(x, y, lag.max = 5,  B = 40, ci.type="mbb") # independence
res4 <- Srho.test.ts.p(x, z, lag.max = 5,  B = 40, ci.type="mbb") # dependence
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
