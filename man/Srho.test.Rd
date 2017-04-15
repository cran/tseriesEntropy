\name{Srho.test}
\alias{Srho.test}
\title{Entropy Test For Serial And Cross Dependence For Categorical Sequences}
\description{
Bootstrap/permutation tests of serial and cross dependence for integer or categorical sequences.
}
\usage{
Srho.test(x, y, lag.max, B = 1000, stationary = TRUE, plot = TRUE,
quant = c(0.95, 0.99), nor = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{integer or factor time series objects or vectors. (\code{y} is missing in the univariate case).}
  \item{lag.max}{maximum lag at which to calculate Srho; default is \code{trunc(N/4)} where N is the number of observations.}
  \item{B}{number of bootstrap/permutation replications.}
  \item{stationary}{logical. If \code{TRUE} assumes stationarity and computes marginal probabilities by using all the N observations. If \code{FALSE} uses N-k observations where k is the lag.}
  \item{plot}{logical. If \code{TRUE}(the default) produces a plot of Srho together with permutation confidence bands under the null hypothesis of independence.}
  \item{quant}{quantiles to be specified for the computation of the significant lags and the plot of confidence bands. Up to 2 quantiles can be specified.
    Defaults are 95\% and 99\%.}
  \item{nor}{logical. If \code{TRUE} normalizes Srho with respect to its attainable maximum. Defaults to \code{FALSE}.}
}

\details{
\describe{
\item{Univariate version: test for serial dependence}{\preformatted{Srho.test(x, lag.max, B = 1000,
stationary = TRUE, plot = TRUE, quant = c(0.95, 0.99), nor = FALSE)}}

\item{Bivariate version: test for cross dependence}{\preformatted{Srho.test(x, y, lag.max, B = 1000,
stationary = TRUE, plot = TRUE, quant = c(0.95, 0.99), nor = FALSE)}}
}
}

\value{An object of class "Srho.test", which is a list with the following elements:
    \item{.Data}{ vector of \code{lag.max} elements containing Srho computed at each lag.}
    \item{quantiles}{Object of class \code{"matrix"}: contains the quantiles of the bootstrap/permutation distribution under the null hypothesis.}
    \item{test.type}{Object of class \code{"character"}: contains a description of the type of test performed.}
    \item{significant.lags}{Object of class \code{"list"}: contains the lags at which Srho exceeds the confidence bands at \code{quant}\% under the null hypothesis.}
    \item{p.value}{Object of class \code{"numeric"}: contains the bootstrap p-value for each lag.}
    \item{lags}{integer vector that contains the lags at which Srho is computed.}
    \item{stationary}{Object of class \code{"logical"}:  \code{TRUE} if the stationary version is computed.}
    \item{data.type}{Object of class \code{"character"}: contains the data type.}
    \item{notes}{Object of class \code{"character"}: additional notes.}
} \section{Warning}{ Unlike \code{ccf} the lag k value returned
by \code{Srho.test(x,y)} estimates Srho between \code{x[t]} and
\code{y[t+k]}. The result is returned invisibly if plot is
TRUE. }
\references{
Granger C. W. J., Maasoumi E., Racine J., (2004) A dependence metric for possibly nonlinear processes.
\emph{Journal of Time Series Analysis}, \bold{25(5)}, 649--669.

Maasoumi E., (1993) A compendium to information theory in economics and econometrics.
\emph{Econometric Reviews}, \bold{12(2)}, 137--181.
}
\author{Simone Giannerini<simone.giannerini@unibo.it>}

\seealso{See also \code{\link{Srho}}, \code{\link{Srho.ts}}. The function \code{\link{Srho.test.ts}} implements the same test for numeric data.}
\examples{
set.seed(12)
x <- as.integer(rbinom(n=30,size=4,prob=0.5))
y <- as.integer(rbinom(n=30,size=4,prob=0.5))
z <- as.integer(c(4,abs(x[-30]*2-2))-rbinom(n=30,size=1,prob=1/2))

# no dependence
Srho.test(x,lag.max=4)   # univariate
Srho.test(x,y,lag.max=4) # bivariate

# lag 1 dependence
Srho.test(x,z,lag.max=4) # bivariate
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
