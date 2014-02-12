\encoding{latin1}
\name{positions_to_regions}
\alias{positions_to_regions}
\title{Converts methylation counts for single genomic positions to counts for genomic regions}
\description{
Converts methylation counts of a data.frame of single genomic positions into a data.frame of counts for genomic regions.
}
\usage{
positions_to_regions(params, outputPath = getwd())
}
\arguments{
\item{params}{BEAT parameter object.}
\item{outputPath}{Path to which output files will be written, the default is the current working directory.}
}
\format{
Necessary function arguments are passed via a BEAT parameter object, which includes working path, sample names, reference sample name, model parameters and region sizes.
\describe{
	\item{params}{Parameter object created by calling \code{\link{makeParams}}.}
	\item{sample.positions.csv}{For each sample referenced by the params argument sampleNames, a corresponding csv must be present in the working directory. The csv input must contain one row per genomic position for each CG site. Columns must be chr, pos, meth, unmeth (chromosome, position, methylated counts, unmethylated counts).}
}
}
%% \value{
%% The function \code{positions_to_regions} returns :
%% \item{regionsList}{A list of cpgregions objects, one for each given sample. These are data.frame objects of the counts per region, consisting of the following columns: chromosome, start position, end position, methylated counts, unmethylated counts.}
%% \item{regions.RData}{For each sample, the data.frame as described above is saved in a regions.RData object in the working directory.}
%% }
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
# Pool CG positions into regions
positions_to_regions(params)
}
\keyword{BEAT}
\author{Kemal Akman <akmank@mpipz.mpg.de>
}
