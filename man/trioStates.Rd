\name{trioStates}
\alias{trioStates}

\title{

  The possible copy number states for a F-M-O trio.

}

\description{

  The possible copy number states for a F-M-O trio.

}

\usage{
trioStates(states = 0:4)
}
\arguments{
  \item{states}{
    The possible copy number states for one member in a trio.
}
}

\value{
  R x 3 matrix, where R is the number of possible trio states and the
  columns correspond to F, M, and O, respectively.

  For instance, a row with 3-3-3 would correspond to
  normal-normal-normal if the states are specified as 0-4 (the third
  index in the 0-4 vector corresponds to copy number 2).

}

\author{
R. Scharpf
}


\examples{
trioStates(0:4)
}
\keyword{manip}
