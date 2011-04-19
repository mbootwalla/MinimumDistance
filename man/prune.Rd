\name{prune}
\alias{prune}
\docType{data}
\title{

  Prune CBS ranges as a function of the absolute difference in adjacent
  segments (in units of standard deviation) and coverage.

}
\description{

  Ranges with low coverage must have a correspondingly higher difference
  in consecutive ranges.  Ranges with high coverage can have relatively
  smaller differences for adjacent segment means, but still be
  statistically significant.  An exponential decay function of the form
  lambda * exp(-lambda * coverage) is used to provide a rough rule of
  the relationship of coverage, distance, and statistical
  significance. Given two adjacent segments, we assess whether the
  vertical distance between their segment means exceeds lambda *
  exp(-lambda * coverage), where coverage refers to the coverage of the
  5-prime (left) segment.

}

\usage{
prune(genomdat, range.object, physical.pos, lambda = 0.05, MIN.CHANGE =
0.1, MAX.CHANGE = 3, MIN.COVERAGE = 3, weighted = FALSE, weights = NULL)
}

\arguments{

  \item{genomdat}{

    The log R ratios for the markers in the \code{range.object}.

}
\item{range.object}{

  The ranges for one individual and one chromosome.

}
\item{physical.pos}{

  Physical position of the markers in the range.object.  The physical
  position should correspond to the ordering of \code{genomedat}.

}
  \item{lambda}{

    The decay constant (> 0).

}
\item{MIN.CHANGE}{

  The minimum number of standard deviations for a shift to be
  biologically meaningful.  This determines how low the decay function
  will go.  See details.

}
\item{MAX.CHANGE}{

  The vertical height where the decay curve begins.  The vertical height
  has units of standard deviations (e.g, 1 = 1 standard deviation of the
  minimum distance).

}
\item{MIN.COVERAGE}{

  All segments with less than MIN.COVERAGE markers will be pruned
  automatically.

}
  \item{weighted}{

    whether to weight the mean of the log r ratios for a segment by
    weights.

  }

    \item{weights}{

      weights for the markers.  See the \code{segment}  function in the
      \code{DNAcopy} package.

  }

}
\details{


}
\value{

  \code{RangedData} object with ranges pruned.

}

\author{
R.Scharpf
}

\seealso{

  \code{\link{segment}}

}

\examples{
	## plot the number of sds required for coverage 3 - 100
	numberSds <- calculateChangeSd(3:100, lambda=0.05, MIN.CHANGE=0.2, MAX.CHANGE=3)
	par(mfrow=c(1,2))
	plot(3:100, numberSds, ylab="number of MADs", xlab="coverage")
	plot(3:100, 0.2 + f(0.05)/.015)
}
\keyword{manip}

