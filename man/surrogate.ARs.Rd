\name{surrogate.ARs}
\alias{surrogate.ARs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Surrogate Time Series Through A Modeling (Smoothed Sieve Bootstrap)}
\description{
Starting from a time series \code{x} given as input, the function generates surrogate series by means of the smoothed sieve bootstrap.
The surrogates can be used for testing for non linearity in time series.
}
\usage{
surrogate.ARs(x, order.max = NULL,
 fit.method = c("yule-walker","burg", "ols", "mle", "yw"), nsurr)
}
\arguments{
  \item{x}{a univariate numeric time series object or a numeric vector.}
  \item{order.max}{maximum order of the AR model to fit. Passed to \code{\link{ar}}.}
  \item{fit.method}{character string giving the method used to fit the AR model. It is passed to \code{\link{ar}} and
  has to be one of the strings in the default argument (partial matching works). Defaults to "yule-walker".}
  \item{nsurr}{number of surrogates.}
}
\details{
Let \code{N} be the length of the series \code{x}. The best AR model is chosen by means of the AIC criterion. Surrogate series are obtained by driving the fitted model with
the smoothed resampled residuals. Smoothing is performed through Kernel density estimation with a Gaussian Kernel by using the
dafaults of \code{\link{density}}.}
\value{
A list with the following elements:
    \item{surr}{a matrix with \code{N} rows and \code{nsurr} columns, in each column is stored a surrogate.}
    \item{call}{contains the call to the routine.}
}
\references{
Giannerini S., Maasoumi E., Bee Dagum E., (2015), Entropy testing
for nonlinear serial dependence in time series, \emph{Biometrika}, \bold{102(3)}, 661--675
\doi{10.1093/biomet/asv007}.

Bickel, P., Buhlmann, P., (1999). A new mixing notion and
functional central limit theorems for a sieve bootstrap in time
series. \emph{Bernoulli} \bold{5}, 413--446. }

\author{Simone Giannerini<simone.giannerini@unibo.it>}

\seealso{See also \code{\link{surrogate.AR}}, \code{\link{Trho.test.AR}}, \code{\link{surrogate.SA}}, \code{\link{Trho.test.SA}}.}
\examples{
set.seed(1345)
# Generates a AR(1) series
x      <- arima.sim(n=120, model = list(ar=0.8));
x.surr <- surrogate.ARs(x, order.max=NULL, nsurr=3);
plot.ts(x.surr$surr,col=4);


## Check that the surrogates have the same ACF of x
corig <- acf(x,10,plot=FALSE)$acf[,,1];
csurr <- acf(x.surr$surr[,1],10,plot=FALSE)$acf[,,1];
round(cbind(corig,csurr,"abs(difference)"=abs(corig-csurr)),3)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
