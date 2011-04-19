\name{computeLoglik}
\alias{computeLoglik}

\title{

  Compute the emission probabilities for the the BAFs and log R ratios
  for each individual in a trio under several possible true copy number
  states.

}

\description{

    Compute the emission probabilities for the the BAFs and log R ratios
  for each individual in a trio under several possible true copy number
  states.

}

\usage{
computeLoglik(id, bsSet, CHR, mu.logr = c(-2, -0.5, 0, 0.3, 0.75), states = 0:4, baf.sds = c(0.02, 0.03, 0.02), THR = -50)
}

\arguments{
  \item{id}{

    The individual ids for members in one trio in the order of F-M-O.
    The ids should match to the sampleNames of the \code{LogRatioSet}
    object.

  }

  \item{bsSet}{

    Object of class \code{LogRatioSet} containing BAFs and logR Ratios
    for the Cleft study.

  }

  \item{CHR}{

    integer specifying which chromosome (this determines the rows of the
    \code{bsSet} object)

}

\item{mu.logr}{

    mu for normal distribution that is used to calculate the emission
    probabilities for log R ratios.  The mus correspond to the copy
    number states 0-4.

}

\item{states}{

    Integer indicating copy number states.

  }

\item{baf.sds}{

      sd of B allele frequency for (1) truncated normal with mu=0, (2)
    normal distribution with mu > 0 and mu < 1, and (3) for truncated
    normal with mu=1. For all normal.

}
  \item{THR}{

        Threshold for emission probability.  Used such that extreme values
	are not overly influential.

      }

    }

\details{

  An object of class 'LikSet' is instantiated that contains slots for
  'loglik' and 'range.index'.  The 'loglik' is a 2 x M x 3 x S array.

  For example,

  loglik(object)[i, j, k, l, m] extracs the emission probabiliities for
  measurement i (1 = log R ratios, 2 = B allele frequencies; i = 1 or
  2), markers j (j can be a numeric vector s.t. 1 <= j <= nrow(object)),
  trio member k (1=F, 2= M, 3=O), and copy number state S (1= null, 2=
  hemizygous deletion, 3=normal, 4=single copy duplication, 5=2 copy
  duplication).

  The minimum distance of the log r ratios between parents and offspring
  were segmented into intervals of constant copy number using circular
  binary segmentation. For instance, one chromosome for a trio may have
  been segmented into 5 intervals with different segment means.  The
  range.index identifier in the \code{featureData} slot provides a
  mapping of each marker to the genomic interval.

}

\value{

  object of class \code{LikSet}.

}

\author{
R. Scharpf
}

\seealso{

}
\examples{

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
