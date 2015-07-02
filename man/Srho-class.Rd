\name{Srho-class}
\docType{class}
\alias{Srho-class}
\alias{plot,Srho,missing-method}
\alias{show,Srho-method}

\title{Class "Srho"}
\description{A class for \code{Srho} and its extensions}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Srho", ...)}.}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"numeric"}: contains Srho computed on the data set.}
    \item{\code{lags}:}{Object of class \code{"integer"}: contains the lags at which Srho is computed.}
    \item{\code{stationary}:}{Object of class \code{"logical"}:  \code{TRUE} if the stationary version is computed.}
    \item{\code{data.type}:}{Object of class \code{"character"}: contains the data type.}
    \item{\code{notes}:}{Object of class \code{"character"}: additional notes.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "Srho", y = "missing")}: ... }
    \item{show}{\code{signature(object = "Srho")}: ... }
	 } }

\author{Simone Giannerini <simone.giannerini@unibo.it> }

\seealso{
See Also \code{\linkS4class{Srho.test}}
}
\examples{
showClass("Srho")
}
\keyword{classes}
