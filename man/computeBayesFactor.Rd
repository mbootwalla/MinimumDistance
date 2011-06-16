\name{computeBayesFactor}
\alias{computeBayesFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute bayes factors for the minimum distance ranges
}
\description{
  Compute bayes factors for the minimum distance ranges
}
\usage{
computeBayesFactor(object, trioSet, id, states = 0:4, baf.sds = c(0.02, 0.03, 0.02), mu.logr = c(-2, -0.5, 0, 0.3, 0.75), log.pi = rep(1/length(states), length(states)), tau, normal.index = 61, a = 9e-04, prGtCorrect = 0.999, d0 = 10, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{trioSet}{
%%     ~~Describe \code{trioSet} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
}
  \item{states}{
%%     ~~Describe \code{states} here~~
}
  \item{baf.sds}{
%%     ~~Describe \code{baf.sds} here~~
}

  \item{mu.logr}{
%%     ~~Describe \code{mu.logr} here~~
}
  \item{log.pi}{
%%     ~~Describe \code{log.pi} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{normal.index}{
%%     ~~Describe \code{normal.index} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{prGtCorrect}{
%%     ~~Describe \code{prGtCorrect} here~~
}
  \item{d0}{
%%     ~~Describe \code{d0} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	library(Beaty)
	library(BeatyExperimentData)
	data(trioSets)
	data(prunedRanges)
	tau <- Beaty:::transitionProbability(states=0:4, epsilon=0.5)
	log.pi <- log(Beaty:::initialStateProbs(states=0:4, epsilon=0.5))
	BATCHSIZE <- 300
	NN <- splitIndicesByLength(seq(length=ncol(trioSets)), BATCHSIZE)
	idList <- lapply(NN, function(i, x) x[i], sampleNames(trioSets))##ssampleNames(minDistanceSet))
##	id <- idList[[10]][88]
##	rd <- computeBayesFactor(object=trioSets, ranges=prunedRanges,
##					 id=ids,
##					 tau=tau, log.pi=log.pi)
	\dontrun{
		## one trio (id = 21153), all chromosomes in trioSets
		Rprof(interval=2)
		system.time(rd <- computeBayesFactor(object=trioSets, ranges=prunedRanges,
						     id="21153",
						     tau=tau, log.pi=log.pi))
		Rprof(NULL)

		## one trio (id = 21153), chromosomes 8 and 22
		rd <- computeBayesFactor(object=trioSets[c(8, 22)], ranges=prunedRanges,
					 id="21153",
					 tau=tau, log.pi=log.pi)

		## two trios in trioSets, two chromosomes
		rd <- computeBayesFactor(object=trioSets[c(8,22)], ranges=prunedRanges,
					 id=sampleNames(trioSets)[1:2],
					 tau=tau, log.pi=log.pi)

		## all trios in trioSets, all chromosomes in trioSets
		rd <- computeBayesFactor(object=trioSets, ranges=prunedRanges,
					 tau=tau, log.pi=log.pi)


	}
	## one trio (id = 21153), chromosome 22
	Rprof(interval=0.5)
	computeBayesFactor(object=trioSets[[10]], ranges=prunedRanges,
			   id="21153",
			   tau=tau, log.pi=log.pi)
	Rprof(NULL)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}

