\encoding{latin1}
\name{makeParams}
\alias{makeParams}
\title{Creates a parameter object of arguments to be used with other BEAT functions.}
\description{
Creates a parameter object of arguments to be used with other BEAT functions.
}
\usage{
makeParams(localpath, sampNames, convrates, is.reference, pminus = 0.2, regionSize = 10000, minCounts = 5, verbose = TRUE, computeRegions = FALSE, computeMatrices = TRUE, writeEpicallMatrix = TRUE, doPlots = FALSE)
}
\arguments{
\item{localpath}{Full path to working directory from which files are read and where results are saved.}
\item{sampNames}{Vector of sample names to be analyzed.}
\item{convrates}{Vector of empirically determined bisulfite conversion rates per sample. Determines p+, the model parameter for incomplete conversion (false negative rates).}
\item{is.reference}{Vector of reference (TRUE) vs. single-cell (FALSE) status per sample.}
\item{pminus}{Model parameter for false conversion (false positive rate).}
\item{regionSize}{Region size in nucleotides into which genomic sites are grouped.}
\item{minCounts}{Minimum counts necessary for each region to be included in epimutation modeling and analysis.}
\item{verbose}{Shows more verbose console output during computation steps.}
\item{computeRegions}{If set to TRUE, regions will be recomputed from individual positions and saved as cpgregions.RData objects for each sample.}
\item{computeMatrices}{If set to TRUE, model parameters will be recomputed and saved as results.RData objects for each sample.}
\item{writeEpicallMatrix}{If set to TRUE, epimutation calls will be written as RData object.}
\item{doPlots}{If set to TRUE, plots of epimutation rates will be generated and saved as PDF in the working directory.}
}
\value{
The function \code{makeParams} returns :
\item{params}{Parameter object to be used in other BEAT functions.}
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
# Example usage of the params object
positions_to_regions(params)
}
\keyword{BEAT}
\author{Kemal Akman <akmank@mpipz.mpg.de>
}
