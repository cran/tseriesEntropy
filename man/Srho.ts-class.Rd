\name{Srho.ts-class}
\docType{class}
\alias{Srho.ts-class}
\alias{show,Srho.ts-method}

\title{Class "Srho.ts"}
\description{A class for \code{Srho} for continuous data/time series.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Srho.ts", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"numeric"}: contains Srho computed on the data set.}
    \item{\code{method}:}{Object of class \code{"character"}: computation method, can be "integral" or "summation".}
    \item{\code{bandwidth}:}{Object of class \code{"character"}: bandwidth selection method.}
    \item{\code{lags}:}{Object of class \code{"integer"}: contains the lags at which Srho is computed.}
    \item{\code{stationary}:}{Object of class \code{"logical"}:  \code{TRUE} if the stationary version is computed.}
    \item{\code{data.type}:}{Object of class \code{"character"}: contains the data type.}
    \item{\code{notes}:}{Object of class \code{"character"}: additional notes.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Srho}"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "Srho.ts")}: ... }
	 }
}
 \author{Simone Giannerini<simone.giannerini@unibo.it>}

\seealso{
See Also \code{\linkS4class{Srho.test}}, \code{\linkS4class{Srho}}
}
\examples{
showClass("Srho.ts")
}
\keyword{classes}
