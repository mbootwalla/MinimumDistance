\name{TrioSet-class}
\Rdversion{1.1}
\docType{class}
\alias{TrioSet-class}

\title{Class \code{"TrioSet"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TrioSet", assayData, phenoData, featureData, experimentData, annotation, protocolData, ...)}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
    \item{\code{phenoData2}:}{Object of class \code{"array"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
No methods defined with class "TrioSet" in the signature.
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
R <- array(NA, dim=c(10, 2, 3))
dimnames(R) <- list(letters[1:10], c("A", "B"), c("F", "M", "O"))
B <- array(runif(10*2*3, 0, 1), dim=dim(R), dimnames=dimnames(R))
x <- new("TrioSet", logRRatio=R, BAF=B)
x[1:2, ]
x[1:2,1]
x[1:2,1,1]
x <- x[, , 1] ## should get all fathers
x <- x[, , 2] ## should get all mothers
x <- x[, , 3] ## should get all offspring

R <- array(NA, dim=c(10, 2, 3))
dimnames(R) <- list(letters[1:10], c("A", "B"), c("F", "M", "O"))
B <- array(runif(10*2*3, 0, 1), dim=dim(R), dimnames=dimnames(R))
fD <- annotatedDataFrameFrom(R, byrow=T)

ldPath("~/")
R <- createFF("r", dim=c(10,2,3), vmode="double")
B <- createFF("r", dim=c(10,2,3), vmode="double")
x2 <- new("TrioSet", logRRatio=R, BAF=B,featureData=fD)



tmp <- list(x, x2)

showClass("TrioSet")
}
\keyword{classes}
