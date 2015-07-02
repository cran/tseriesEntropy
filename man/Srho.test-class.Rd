\name{Srho.test-class}
\docType{class}
\alias{Srho.test-class}
\alias{plot,Srho.test,missing-method}
\alias{show,Srho.test-method}

\title{Class "Srho.test"}
\description{A class of tests for serial dependence and nonlinearity based upon \code{Srho}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Srho.test", ...)}.
}
\section{Slots}{
     \describe{
    \item{\code{.Data}:}{Object of class \code{"numeric"}: contains Srho computed on the data set.}
    \item{\code{call}:}{Object of class \code{"call"}: contains the call to the routine.}
    \item{\code{call.h}:}{Object of class \code{"call"}: contains the call to the routine used for obtaining the surrogates or the bootstrap replicates under the null hypothesis.}
    \item{\code{quantiles}:}{Object of class \code{"matrix"}: contains the quantiles of the bootstrap/permutation distribution under the null hypothesis.}
    \item{\code{test.type}:}{Object of class \code{"character"}: contains a description of the type of test performed.}
    \item{\code{significant.lags}:}{Object of class \code{"list"}: contains the lags at which Srho exceeds the confidence bands at \code{quant} under the null hypothesis.}
    \item{\code{p.value}:}{Object of class \code{"numeric"}: contains the bootstrap p-value for each lag.}
    \item{\code{lags}:}{Object of class \code{"integer"}: contains the lags at which Srho is computed.}
    \item{\code{stationary}:}{\code{TRUE} if the stationary version is computed.}
    \item{\code{data.type}:}{Object of class \code{"character"}: contains the data type.}
    \item{\code{notes}:}{Object of class \code{"character"}: additional notes.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Srho}"}, directly.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "Srho.test", y = "missing")}: ... }
    \item{show}{\code{signature(object = "Srho.test")}: ... }
     }
}
\author{ Simone Giannerini <simone.giannerini@unibo.it> }

\seealso{
See Also \code{\linkS4class{Srho}}
}
\examples{
showClass("Srho.test")
}
\keyword{classes}
