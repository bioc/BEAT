\encoding{latin1}
\name{epimutation_calls}
\alias{epimutation_calls}
\title{Returns epimutation rates and sites.}
\description{
Returns epimutation rates per genome and per genomic feature, as well as individual genomic sites at which epimutations were called.
}
\usage{
epimutation_calls(params, outputPath = getwd())
}
\arguments{
\item{params}{BEAT parameter object.}
\item{outputPath}{Path to which output files will be written, the default is the current working directory.}
}
\format{
Necessary function arguments are passed via a BEAT parameter object, which includes working path, sample names, reference sample name, model parameters and region sizes.
\describe{
	\item{params}{Parameter object created by calling \code{\link{makeParams}}.}
}
}
\value{
The function \code{epimutation_calls} returns :
\item{resultsList}{ A list is returned consisting of the two data.frames methSites and demethSites. methSites contains all regions at which methylating epimutations were observed, while demethSites contains all sites at which demethylating epimutations were observed. Each data.frame describes the genomic regions covered by the given sample and the reference sample using the columns: 'chr' (chromosome), 'start' (starting position), 'stop' (last position), 'meth' (methylated counts), 'unmeth' (unmethylatd counts) and 'epimutation_call_test' (epimutation call, 1 for methylating epimutation and $-1$ for demethylating epimutation). }
\item{methEpicalls}{ For each single-cell sample, a methEpicalls.RData object is saved in the working directory, which is a data.frame of all sites at which methylating epimutations were observed, consisting of the columns: chr, pos, endpos, meth, unmeth and methstate. }
\item{demethEpicalls}{ For each single-cell sample, a demethEpicalls.RData object is saved in the working directory, which is a data.frame of all sites at which demethylating epimutations were observed, consisting of the columns: chr, pos, endpos, meth, unmeth and methstate. }
}
\seealso{
See also \code{\link{makeParams}}.
}
\examples{
# Local working directory
localpath <- system.file('extdata', package = 'BEAT')
# Names of samples, expected filenames are e.g. reference.positions.csv
sampNames <- c("reference", "sample")
# Empirical BS-conversion rates, e.g. estimated from non-CpG methylation
convrates <- c(0.8,0.5)
# Vector denoting reference vs. single-cell status of given samples
is.reference <- c(TRUE,FALSE)
params <- makeParams(localpath, sampNames, convrates, is.reference, pminus = 0.2, regionSize = 10000, minCounts = 5, verbose = TRUE, computeRegions = TRUE, computeMatrices = TRUE, writeEpicallMatrix = TRUE)
# pool CG positions into regions
positions_to_regions(params)
# compute model statistics
generate_results(params)
# call epimutations
methDemethPlusMinus <- epimutation_calls(params)
}
\keyword{BEAT}
\author{Kemal Akman <akmank@mpipz.mpg.de>
}
