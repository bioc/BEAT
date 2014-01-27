library(BEAT)
localpath = system.file('example', package = 'BEAT')
positions = read.csv(paste(localpath, "/", "singlecell.positions.csv", sep=""))
head(positions)

sampNames <- c('reference','singlecell')
sampNames
is.reference <- c(TRUE, FALSE)

convrates <- c(0.8,0.5)

params = makeParams(localpath, sampNames, convrates, is.reference, pminus = 0.2, regionSize = 10000, minCounts = 5, verbose = TRUE, computeRegions = TRUE, computeMatrices = TRUE, writeEpicallMatrix = TRUE)
params

regions <- positions_to_regions(params)
head(regions[[1]])

results <- generate_results(params)
class(results)
head(results[[1]])

epiCalls <- epimutation_calls(params)
head(epiCalls$methSites$singlecell,3)
head(epiCalls$demethSites$singlecell,3)
