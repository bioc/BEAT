\encoding{latin1}
\name{generate_results}
\alias{generate_results}
\title{Computes model methylation states for genomic region counts of a list of samples.}
\description{
Computes model methylation states for genomic region counts of a list of samples.
}
\usage{
generate_results(params)
}
\arguments{
\item{params}{BEAT parameter object.}
}
\format{
Necessary function arguments are passed via a BEAT parameter object, which includes working path, sample names, reference sample name, model parameters and region sizes.
\describe{
	\item{params}{Parameter object created by calling \code{\link{makeParams}}.}
}
}
\value{
The function \code{generate_results} returns :
\item{cpgregionList}{ A list of data.frame objects, one for each sample. Each data.frame consists of one row per genomic region and 9 columns of chromosome, position, methylated counts, unmethylated counts and methylation state data computed by the model. }
\item{results}{ For each sample, the data.frame as described above is saved in a results.RData object in the working directory. }
}
\seealso{
See also \code{\link{makeParams}}.
}
\examples{
# Local working directory
localpath = system.file('example', package = 'BEAT')
# Names of samples, expected filenames are e.g. reference.positions.csv
sampNames <- c("reference", "singlecell")
# Empirical BS-conversion rates, e.g. estimated from non-CpG methylation
convrates <- c(0.8,0.5)
# Vector denoting reference vs. single-cell status of given samples
is.reference <- c(TRUE,FALSE)
params = makeParams(localpath, sampNames, convrates, is.reference, pminus = 0.2, regionSize = 10000, minCounts = 5, verbose = TRUE, computeRegions = TRUE, computeMatrices = TRUE, writeEpicallMatrix = TRUE, doPlots = TRUE)
# pool CG positions into regions
positions_to_regions(params)
# compute model statistics
generate_results(params)
}
\keyword{BEAT}
\author{Kemal Akman <akmank@mpipz.mpg.de>
}
