\name{surrogate.AR}
\alias{surrogate.AR}
\title{Surrogate Time Series Through AR Modeling (Sieve Bootstrap)}
\description{
Starting from a time series \code{x} given as input, the function generates surrogate series by means of the sieve bootstrap.
The surrogates can be used for testing for non linearity in time series.
}
\usage{
surrogate.AR(x, order.max = 10, fit.method = c("yule-walker",
 "burg", "ols", "mle", "yw"), nsurr)
}
\arguments{
  \item{x}{a univariate numeric time series object or a numeric vector.}
  \item{order.max}{maximum order of the AR model to fit.}
  \item{fit.method}{character string giving the method used to fit the AR model.
  Must be one of the strings in the default argument (the first few characters are sufficient). Defaults to "yule-walker".}
  \item{nsurr}{number of surrogates.}
}
\details{
\code{N} is the length of the series \code{x}. The best AR model is chosen by means of the AIC criterion. The residuals of the model are
resampled with replacement. Surrogate series are obtained by driving the fitted model with the resampled residuals.
}
\value{
A list with the following elements:
    \item{surr}{a matrix with \code{N} rows and \code{nsurr} columns, in each column is stored a surrogate.}
    \item{call}{contains the call to the routine.}
}
\references{
Giannerini S., Maasoumi E., Bee Dagum E., (2015), Entropy testing
for nonlinear serial dependence in time series, \emph{Biometrika}, \bold{102(3)}, 661--675
\url{http://doi.org/10.1093/biomet/asv007}.

Buhlmann, P., (1997). Sieve bootstrap for time series.
\emph{Bernoulli}, \bold{3}, 123--148.
}

\author{Simone Giannerini<simone.giannerini@unibo.it>}

\seealso{See also \code{\link{surrogate.AR}}, \code{\link{Trho.test.AR}}, \code{\link{surrogate.SA}}, \code{\link{Trho.test.SA}}.}
\examples{
set.seed(1345)
# Generates a AR(1) series
x      <- arima.sim(n=120, model = list(ar=0.8));
x.surr <- surrogate.AR(x, order.max=10, nsurr=3);
plot.ts(x.surr$surr,col=4);


## Check that the surrogates have the same ACF of x
corig <- acf(x,10,plot=FALSE)$acf[,,1];
csurr <- acf(x.surr$surr[,1],10,plot=FALSE)$acf[,,1];
round(cbind(corig,csurr,"abs(difference)"=abs(corig-csurr)),3)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
