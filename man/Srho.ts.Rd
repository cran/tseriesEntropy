\name{Srho.ts}
\alias{Srho.ts}
\title{Entropy Measure Of Serial And Cross Dependence}
\description{
Entropy based measure of serial and cross dependence for continuous data. For integer/categorical data see \code{\link{Srho}}.
Implements a normalized version of the Hellinger/Matusita distance. As shown in the references the metric measure is a proper distance.}
\usage{
Srho.ts(x, y, lag.max = 10, bw = c("reference", "mlcv", "lscv", "scv", "pi"),
 method = c("integral", "summation"), plot = TRUE,
 maxpts = 0, tol = 0.001)
}
\arguments{
  \item{x, y}{univariate numeric time series object or numeric vectors (\code{y} is missing in the univariate case).}
  \item{lag.max}{maximum lag at which to calculate Srho; default is 10}
  \item{bw}{Object of class \code{"character"}: bandwidth selection method, can be \code{"reference"}, \code{"mlcv"}, \code{"lscv"}, \code{"scv"}, \code{"pi"}.}
  \item{method}{Object of class \code{"character"}: computation method, can be \code{"integral"} or \code{"summation"}.}
  \item{plot}{logical. If \code{TRUE} (the default) Srho is plotted.}
  \item{maxpts}{maximum number of function evaluations, passed to \code{\link[cubature]{adaptIntegrate}}, defaults to 0 (unlimited).}
  \item{tol}{max. tolerance, passed to \code{\link[cubature]{adaptIntegrate}}.}
}
\details{
\describe{
\item{Univariate version: serial entropy}{\preformatted{Srho.ts(x, lag.max = 10,
 bw = c("reference", "mlcv", "lscv", "scv", "pi"),
 method = c("integral", "summation"), plot = TRUE, maxpts = 0, tol = 0.001)}}

\item{Bivariate version: cross entropy}{\preformatted{Srho.ts(x, y, lag.max = 10,
 bw = c("reference", "mlcv", "lscv", "scv", "pi"),
 method = c("integral", "summation"), plot = TRUE, maxpts = 0, tol = 0.001)}}
}
The bandwidth selection methods are the following:
\describe{
  \item{\code{reference}:}{reference criterion.}
  \item{\code{mlcv}:}{maximum likelihood cross-validation.}
  \item{\code{lscv}:}{least-squares cross-validation, see \code{\link[ks]{Hlscv}}.}
  \item{\code{scv}:}{smoothed cross-validation, see \code{\link[ks]{Hscv}}}
  \item{\code{pi}:}{plugin, see \code{\link[ks]{Hpi}}}
}
}
\value{An object of class "Srho.ts", with the following slots:
    \item{.Data}{Object of class \code{"numeric"}: contains Srho computed on the data set.}
    \item{method}{Object of class \code{"character"}: computation method}
    \item{bandwidth}{Object of class \code{"character"}: bandwidth selection method.}
    \item{lags}{Object of class \code{"integer"}: contains the lags at which Srho is computed.}
    \item{stationary}{Object of class \code{"logical"}:  \code{TRUE} if the stationary version is computed.}
    \item{data.type}{Object of class \code{"character"}: contains the data type.}
    \item{notes}{Object of class \code{"character"}: additional notes.}
}
\references{
Granger C. W. J., Maasoumi E., Racine J., (2004) A dependence metric for possibly nonlinear processes.
\emph{Journal of Time Series Analysis}, \bold{25(5)}, 649--669.

Maasoumi E., (1993) A compendium to information theory in economics and econometrics.
\emph{Econometric Reviews}, \bold{12(2)}, 137--181.

Giannerini S., Maasoumi E., Bee Dagum E., (2015), Entropy testing
for nonlinear serial dependence in time series, \emph{Biometrika}, \bold{102(3)}, 661--675
\url{http://doi.org/10.1093/biomet/asv007}.
}
 \author{Simone Giannerini<simone.giannerini@unibo.it>}

\seealso{\code{\link{Srho.test.ts}}, \code{\link[cubature]{adaptIntegrate}}, \code{\link[ks]{ks}}.
The function \code{\link{Srho}} implements the same measure for integer/categorical data.}
\examples{
set.seed(11)
x <- arima.sim(list(order = c(1,0,0), ar = 0.8), n = 50)
S <- Srho.ts(x,lag.max=5,method="integral",bw="mlcv")

# creates a nonlinear dependence at lag 1
y <- c(runif(1),x[-50]^2*0.8-0.3)
S <- Srho.ts(x,y,lag.max=3,method="integral",bw="mlcv")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
