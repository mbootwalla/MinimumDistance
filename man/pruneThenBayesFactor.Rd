\name{pruneThenBayesFactor}
\alias{pruneThenBayesFactor}

\title{

  Prune the ranges and then compute Bayes factors on the pruned ranges.

}

\description{

  CBS applied to the minimum distance of the log R ratios produces a set
  of genomic intervals for the offspring.  This function first prunes
  the set of ranges as a function of the amplitude of the distance and
  the coverage (see the prune function for additional details).  For the
  set of pruned ranges, we compute a bayes factor comparing the most
  likely state for the trio to the normal state (F, M, O normal).

  For computing the Bayes factor, we assume that the BAFs are
  independent of the log R ratios conditional on the underlying hidden
  state.  Like PennCNV, we assume that the log R ratios for the F-M-O
  trio are conditionally independent and that the BAFs for the F-M-O
  trio are conditionally independent.

  As the boundaries of the intervals are assumed to be known from the
  segmentation of the minimum distance, the transition probabilities
  within an interval cancel in the numerator and denominator (assuming
  that the probability of staying in the same state is the same for
  different states).  Hence, the only thing we need to factor for the
  alternative is the probability of transitioning out of the normal
  state and the probability of transitioning back to the normal
  state. For the first segment in a chromosome, we replace the first
  transition probability with the an initial state probability.  Hence
  for each segment, we calculate a likelihood that is the sum of the log
  emission probabilities for each marker multiplied by 'tau1' and
  'tau2', where tau1 and tau2 for alternative states are the
  probabilities of leaving the normal state to the alternative state and
  tau2 is the probability of transitioning back to the normal state.
  For the denominator in the bayes factor (the normal state),
  tau1=tau2=Pr(normal state at marker t+1 given normal state at marker
  t).  For the first segment in a sample, we add the log of the intial
  state probabilities to the log emission probability.

}


\usage{
pruneThenBayesFactor(CHR, bsSet, minDistanceSet, ids, lambda = 0.05, MIN.CHANGE = 0.1, MAX.CHANGE = 3, MIN.COVERAGE = 10, states = 0:4, baf.sds = c(0.02, 0.03, 0.02), THR = -50, mu.logr = c(-2, -0.5, 0, 0.3, 0.75), log.pi, tau, normal.index = 61, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

    \item{CHR}{
      integer indicating the chromosome
}

  \item{bsSet}{
    object of class LogRatioSet
  }

  \item{minDistanceSet}{

    Object of class CopyNumberSet containing the minimum distance.  See
    details below.

  }

  \item{ids}{

    First five letters of identifier for offspring (this is sufficient
    to uniquely identify an offspring for the Cleft data).  We will only
    prune the genomic ranges for these individuals.  This is primarily
    for fast testing of the function on trios of interest.

  }

  \item{lambda}{
    passed to \code{prune}
}


  \item{MIN.CHANGE}{
    passed to \code{prune}
}
  \item{MAX.CHANGE}{
    passed to \code{prune}
}
  \item{MIN.COVERAGE}{
    passed to \code{prune}
  }


  \item{states}{
    Integer indicating copy number states.  For instance, 0:4.  Passed
    to \code{computeLoglik} and \code{trioStates}.

}
\item{baf.sds}{

Passed to \code{computeLoglik}.

  }

  \item{THR}{

    Passed to \code{computeLoglik}

}

  \item{mu.logr}{

Passed to \code{computeLoglik}.

}

\item{log.pi}{

  The initial state probabaility -- only relevant for the first segment
  in a chromosome.  The length of log.pi must the the same as the number
  of states and should sum to 1.


}

\item{tau}{

  A transition probability matrix.  tau[i,j] is the probability of
  transitioning from state i to state j.  Passed to \code{computeLoglik}.

}

\item{normal.index}{

  The row index of the 'normal state' in the matrix returned by
  \code{trioStates}. See the \code{trioStates} function for additional
  details.

}
  \item{verbose}{

    Whether to print messages during processing.

}
}
\details{

  Describe how the minimum distance is computed in the CopyNumberSet object.

}

\value{

An object of class RangedData with additional columns indicating the
posterior odds ('bayes.factor') and the argmax of the triostates. The
argmax corresponds to the row index in \code{trioStates()}.

}

\author{
  R. Scharpf
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{prune}}, \code{\link{trioStates}}, \code{\link{computeLoglik}}
}

\examples{
data(bsSet)
trio.states <- trioStates()
tau <- transitionProbability2(states=0:4, epsilon=0.5)
log.pi <- log(initialStateProbs(states=0:4, epsilon=0.5))
stopifnot(sum(exp(log.pi))==1)
ranges.bf <- pruneThenBayesFactor(CHR=22,
				  bsSet=bsSet,
				  minDistanceSet=minDistanceSet,
				  lambda=0.05,
				  ids="22048",
				  MIN.CHANGE=0.5,
				  MAX.CHANGE=3,
				  MIN.COVERAGE=10,
				  states=0:4,
				  baf.sds=c(0.02, 0.03, 0.02),
				  THR=-50,
				  mu.logr=c(-2, -0.5, 0, 0.3, 0.75),
				  log.pi=log.pi,
				  tau=tau,
				  a=0.0009,
				  normal.index=61)
}

\keyword{classif}

