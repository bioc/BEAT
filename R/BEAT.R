#
# BEAT: BS-Seq Epimutation Analysis Toolkit
#
# Version 0.99.1
#
# Written 2012-2013 by Kemal Akman <akman@mpipz.mpg.de>
#
# This software is licensed according to version 3.0 of the GNU LESSER GENERAL PUBLIC LICENSE
# The license text is available under LICENSE or via http://www.gnu.org/licenses/lgpl-3.0.txt
#


###
### utility function section
###

positionsReadCSV <- function(csvFilename) {
	read.csv(csvFilename)
}

positionsWriteCSV <- function(positions, csvFilename) {
	write.csv(positions, file=csvFilename, row.names=FALSE)
}

minCovRegions <- function(cpgregions, minCov) {
	length(which(cpgregions$meth + cpgregions$unmeth >= minCov))
}

cpg2bed <- function(infile, outfile) {
	load(infile)
	positions <- addChrPosMeth(positions)
	positions$methLevel = round(positions$methLevel * 1000)
	write.table(positions[,c("chr","pos","pos","methLevel")],quote=FALSE,sep="\t",file=outfile,row.names=F,col.names=F)
}

cpgAndLogFc2TSV <- function(positions, outfile) {
	write.table(positions[,c("chr","pos","bj3","bj3bleo","logFc")],quote=FALSE,sep=", ",file=outfile,row.names=F,col.names=F)
}

methPerSiteBoth <- function(infile, outfilePrefix="lane1.cpg", skipSAMLines=1) {
	cpg <- methPerSite(infile,"z","Z", skipSAMLines)
	positions <- sumStrands(cpg)
	save(cpg, file=paste(outfilePrefix,".RData",sep=""))
	save(positions, file=paste(outfilePrefix,"both.RData",sep=""))
}

methPerSiteBothSplitData <- function(infileList, outfilePrefix="lane1.cpg", skipSAMLines=1) {
	splitCpg <- list()
	splitCpgboth <- list()
	i <<- 1

	for(infile in infileList) {
		cat(paste("Processing file: ", infile, "\n",sep=""))
		splitCpg[[i]] <- methPerSite(infile,"z","Z", skipSAMLines)
		splitCpgboth[[i]] <- sumStrands(splitCpg[[i]])
		
		i <<- i + 1
	}

	cpg <- splitCpg[[1]]
	positions <- splitCpgboth[[1]]

	cat("Appending...\n")
	for(j in 2:i) {
		cpg <- rbind(cpg, splitCpg[[j]])
		positions <- rbind(positions, splitCpgboth[[j]])
	}
	save(cpg, file=paste(outfilePrefix,".RData",sep=""))
	save(positions, file=paste(outfilePrefix,"both.RData",sep=""))
	cat("methPerSiteBothSplitData done\n")
}

methPerSite <- function(file, nonMeth="z", Meth="Z", skipSAMLines=1) {
	out <- data.frame(stringsAsFactors=F)

	data <- read.delim(file,skip=skipSAMLines, stringsAsFactors=F)
	colnames(data) <- c("spot", "strand", "chr", "pos", "meth")
	chrs <- sort(as.character(unique(data$chr)))
	for(a in chrs)
	{
		gc()
		cat("Summing sites for ", a, "\n", sep="")
		xc <- data[data$chr == a,]
		xcpos <- sort(unique(xc$pos))
		npos <- length(xcpos)

		### 1 line for crick (+), 1 line for watson (-)
		outc <- data.frame(strand=c(rep("+",npos),rep("-",npos)), chr=rep(a,npos*2), pos=c(xcpos,xcpos), meth=rep(0,npos*2), unmeth=rep(0,npos*2), stringsAsFactors=F)
		# apply() that runs once for every unique position on this chromosome
		outc2 <- apply(outc, 1, function(x) {
			results <- xc[as.integer(xc$pos) == as.integer(x[3]),]
			results <- results[as.character(results$strand) == as.character(x[1]),] ;# result := all reads in SAM file matching the current position
			x[4] <- as.character(dim(results[results$meth == Meth,])[1]) ;# add count of methylated C's from result to outc2
			x[5] <- as.character(dim(results[results$meth == nonMeth,])[1]) ;# add count of unmethylated C's from result to outc2
			x
		})
		outc2 <- data.frame(t(outc2), stringsAsFactors=F)
		stopifnot((sum(as.integer(outc2$meth) + as.integer(outc2$unmeth))) == dim(xc)[1])
		out <- rbind(out, outc2)
		rm(xc); rm(outc2)
	}

	colnames(out) <- c("strand", "chr", "pos", "meth", "unmeth")
	out$meth <- as.integer(out$meth)
	out$unmeth <- as.integer(out$unmeth)
	out
}

sumStrands <- function(mframe) {
	out <- data.frame(stringsAsFactors=F)

	chrs <- sort(as.character(unique(mframe$chr)))
	for(a in chrs)
	{
		cat(paste("Counting sites for ",a,"\n",sep=""))
		xc <- mframe[mframe$chr == a,]
		xcpos <- sort(unique(xc$pos))
		npos <- length(xcpos)

		# strand insensitive output
		outc <- data.frame(chr=rep(a,npos), pos=c(xcpos), meth=rep(0,npos), unmeth=rep(0,npos), stringsAsFactors=F)
		outc2 <- apply(outc, 1, function(x) {
			results <- xc[as.integer(xc$pos) == as.integer(x[2]),]
			x[3] = sum(as.integer(as.character(results$meth)))
			x[4] = sum(as.integer(as.character(results$unmeth)))
			x
		})
		outc2 <- data.frame(t(outc2), stringsAsFactors=F)
		stopifnot((sum(as.integer(outc2$meth) + as.integer(outc2$unmeth))) == sum(as.integer(xc$meth) + as.integer(xc$unmeth)))
		out <- rbind(out, outc2)
	}	

	colnames(out) <- c("chr", "pos", "meth", "unmeth")
	out
}

### CpG sites supported by at least n reads
supported <- function(cpg, atLeast=1) {
	dim(filterMinCoverage(cpg, atLeast))[1]
}

supportedExact <- function(cpg, whichCoverage=1) {
	nrow(cpg[which((as.integer(cpg$unmeth) + as.integer(cpg$meth)) == whichCoverage),])
}

filterMinCoverage <- function(cpg, atLeast=3) {
	cpg[which((as.integer(cpg$unmeth) + as.integer(cpg$meth)) >= atLeast),]
}

filterMaxCoverage <- function(cpg, atMost=200) {
	cpg[which((as.integer(cpg$unmeth) + as.integer(cpg$meth)) <= atMost),]
}

### add "chr.pos" field and methylation %
addChrPosMeth <- function(cpg) {
	cpg$pos <- as.integer(cpg$pos)
	cpg <- cbind(cpg, site=paste(cpg$chr, cpg$pos, sep="."), methLevel=(as.double(cpg$meth) / (as.double(cpg$meth)+as.double(cpg$unmeth))), stringsAsFactors=F)
	cpg
}

### add "chr.pos" field and methylation %
addChrPosMethRegions <- function(cpg) {
	cpg$start <- as.integer(cpg$start)
	cpg$stop <- as.integer(cpg$stop)
	cpg <- cbind(cpg, site=paste(cpg$chr, cpg$start, cpg$stop, sep="."), methLevel=(as.double(cpg$meth) / (as.double(cpg$meth)+as.double(cpg$unmeth))), stringsAsFactors=F)
	cpg
}

### filter to fully methylated (default definition > 99%) and fully unmethylated (default definition < 1%) sites
filterHighMethylation <- function(cpg, hiMeth=0.90) {
	cpg[which(cpg$methLevel >= hiMeth),]
}

filterLowMethylation <- function(cpg, loMeth=0.10) {
	cpg[which(cpg$methLevel <= loMeth),]
}

csv2cpgBoth <- function(siteCountFile, saveAs=NULL) {
	rawcpg <- read.delim2(siteCountFile,stringsAsFactors=F)
	meth <- rawcpg[,2]
	unmeth <- rawcpg[,3]

	chr <- unlist(lapply(strsplit(rawcpg[,1],split="\\."), function(x) {x[1]}))
	pos <- unlist(lapply(strsplit(rawcpg[,1],split="\\."), function(x) {x[2]}))
	positions <- data.frame(chr, pos, meth, unmeth, stringsAsFactors=F)
	if(!is.null(saveAs)) { save(positions, file=saveAs) }
}

csv2cpgStrandSpecific <- function(siteCountFile, saveAs=NULL) {
	rawcpg <- read.delim2(siteCountFile,stringsAsFactors=F)
	meth <- rawcpg[,2]
	unmeth <- rawcpg[,3]

	strand <- unlist(lapply(strsplit(rawcpg[,1],split="\\."), function(x) {x[1]}))
	chr <- unlist(lapply(strsplit(rawcpg[,1],split="\\."), function(x) {x[2]}))
	pos <- unlist(lapply(strsplit(rawcpg[,1],split="\\."), function(x) {x[3]}))
	cpgeach <- data.frame(chr, pos, meth, unmeth, strand, stringsAsFactors=F)
	if(!is.null(saveAs)) { save(cpgeach, file=saveAs) }
}

getCPGOverlap <- function(cpgFileOne, cpgFileTwo) {
	load(cpgFileOne)
	one <- addChrPosMeth(positions)
	rm(positions)
	load(cpgFileTwo)
	two <- addChrPosMeth(positions)

	common <- intersect(one$site, two$site)
	onec <- one[which(one$site %in% common),]
	onec <- onec[order(onec$site),]

	twoc <- two[which(two$site %in% common),]
	twoc <- twoc[order(twoc$site),]

	cpg <- data.frame(chr=onec$chr, pos=onec$pos, methA=as.numeric(onec$meth), totalA=as.numeric(onec$meth + onec$unmeth), methB=as.numeric(twoc$meth), totalB=as.numeric(twoc$meth + twoc$unmeth))
	cpg
}

###
### region segmentation section
###
library("GenomicRanges")

reduceToBinary <- function(cpg) {
	cpg <- cpg[which((cpg$unmeth == 0) | (cpg$meth == 0)),]
	cpg[which(cpg$unmeth == 0),]$meth = 1
	cpg[which(cpg$meth == 0),]$unmeth = 1
	cpg
}

cpg2regions <- function(cpg, windowSize, minCountsPerRegion) {
	cpg <- filterMaxCoverage(cpg, 100) ;# remove PCR artifact peaks
	chrs <- unique(cpg$chr)

	regionCountsAll <- data.frame()

	for(chr in chrs) {
		oneChrCpg <- cpg[cpg$chr == chr,]
		oneChrCpg$pos <- as.integer(oneChrCpg$pos)
		oneChrCpg <- oneChrCpg[order(oneChrCpg$pos),]

		### index regions
		cbin = oneChrCpg$pos %/% windowSize
		meth = as.integer(oneChrCpg[,"meth"])
		unmeth = as.integer(oneChrCpg[,"unmeth"])
		
		regionsMeth = tapply(meth,INDEX = as.factor(cbin),sum)
		regionsUnmeth = tapply(unmeth,INDEX = as.factor(cbin),sum)

		regioncounts <- data.frame(chr=rep(chr, length(regionsMeth)), start=as.integer(names(regionsMeth))*windowSize, stop=(as.integer(names(regionsMeth))*windowSize)+windowSize, meth=regionsMeth, unmeth=regionsUnmeth)

		### cutoff: include regions #CGs >= n
		regioncounts <- regioncounts[which(as.integer(regioncounts$meth) + as.integer(regioncounts$unmeth) >= minCountsPerRegion),]
		regionCountsAll <- rbind(regionCountsAll, regioncounts)
	}

	regionCountsAll
}

library("ShortRead")
library("Biostrings")
library("BSgenome")
library("IRanges")

####
#### ASSOCIATE BAM (RAW READ COUNTS) GENOMIC RANGES WITH GENES
####

bam2counts <- function(path, file, BAMtype="BAM", strandData=F, wildcardInput=F) {
	if(!wildcardInput) { file <- paste(file,"$",sep="") } ;# allow only exact matches of the file, no patterns

	aln <- readAligned(path, file, BAMtype)
	mapped <- !is.na(position(aln))
	aln <- aln[mapped]
	readLen <- width(sread(aln))
	endPos <- position(aln)+readLen ;# accounts for individual read len differences, e.g. due to clipping adaptor sequences

	### alntab <- data.frame(start=position(aln),end=position(aln)+readLen,chr=chromosome(aln))
	if(strandData) {
		#alntab <- table(paste(position(aln),endPos,chr=chromosome(aln),strands=GenomicRanges::strand(aln),sep="."))
		alntab <- table(paste(position(aln),endPos,chr=chromosome(aln),strands=strand(aln),sep="."))
	} else {
		alntab <- table(paste(position(aln),endPos,chr=chromosome(aln),sep="."))
	}

	getStartPos <- sapply(strsplit(names(alntab),split="\\."),function(x) { x[1] })
	getEndPos <- sapply(strsplit(names(alntab),split="\\."),function(x) { x[2] })
	getChrs <- sapply(strsplit(names(alntab),split="\\."),function(x) { x[3] })
	getStrands <- sapply(strsplit(names(alntab),split="\\."),function(x) { x[4] })

	if(strandData) {
		alnframe <- data.frame(counts=as.integer(alntab), start=as.integer(getStartPos), end=as.integer(getEndPos), chr=getChrs, strands=getStrands, stringsAsFactors=F)
	} else {
		alnframe <- data.frame(counts=as.integer(alntab), start=as.integer(getStartPos), end=as.integer(getEndPos), chr=getChrs, stringsAsFactors=F)
	}
	## length(which(alnframe$counts > 1))
	alnframe
}

bam2genomeInterval <- function(path, file, BAMtype="BAM") {
	aframe <- bam2counts(path, file, BAMtype)
	new("Genome_intervals", matrix(t(aframe[,2:3]),byrow=T,ncol=2),
			annotation = data.frame(seq_name = factor(aframe$chr), inter_base = rep(F,dim(aframe)[1]), counts=aframe$counts, stringsAsFactors=F))
}

###
### finding and analyzing clusters of imprinted genes
###
findClusteredGenes <- function(AVCandidates, maxGeneClusterDistance) {
	geneIsClustered <- rep(F, times=dim(AVCandidates[[1]])[1])

	## lists of vectors of genes, one vector per chromosome
	geneByChr <- list(); geneByChrIdx <- list(); j <- 1
	chrNames <- sort(unique(AVCandidates[[1]]$chr))

	# for each chromosome, find the positions and their global index in AVCandidates
	for(i in chrNames) {
		geneByChr <- append(geneByChr,0)
		geneByChrIdx[[j]] <- which(AVCandidates[[1]]$chr == i)
		geneByChr[[j]] <- AVCandidates[[1]][which(AVCandidates[[1]]$chr == i),]$chrpos
		j <- j+1
	}

	# for each chromosome, determine which genes are clustered = less than maxGeneClusterDistance apart
	for(i in 1:length(geneByChr)) {
		posInChr <- geneByChr[[i]]

		# distance of all vs. all gene positions on the same chromosome
		distmat <- sapply(posInChr,function(x) {abs(x-posInChr)})

		isClustered <- c()
		if(length(posInChr) > 1) {
			for(k in 1:length(posInChr)) { isClustered[k] <- length(which(distmat[,k] < maxGeneClusterDistance)) > 1 }
		} else if(length(posInChr) == 1) {
			isClustered[1] <- length(which(distmat < maxGeneClusterDistance)) > 1
		} else {
			isClustered <- F
		}

		clusteredIdx <- geneByChrIdx[[i]][which(isClustered)]
		geneIsClustered[clusteredIdx] <- T
	}

	lapply(AVCandidates, function(x) { cbind(x,isClustered=geneIsClustered) })
}

getCluster <- function(geneName, AVCandidates, maxGeneClusterDistance) {
	geneIsClustered <- rep(F, times=dim(AVCandidates[[1]])[1])

	queryPos <- AVCandidates[[1]][which(AVCandidates[[1]]$gene == geneName),]$chrpos
	queryChr <- AVCandidates[[1]][which(AVCandidates[[1]]$gene == geneName),]$chr

	sameChrIdx <- which(AVCandidates[[1]]$chr == queryChr)
	posInChr <- AVCandidates[[1]][sameChrIdx,]$chrpos

	isClustered <- c()

	# distance of the queried gene to all others on the same chromosome
	distmat <- sapply(posInChr,function(x) {abs(x-queryPos)})
	for(k in 1:length(posInChr)) { isClustered[k] <- (distmat[k] < maxGeneClusterDistance) }

	SNPNames <- AVCandidates[[1]][sameChrIdx[isClustered],]$gene
	geneNames <- unique(unlist(lapply(strsplit(SNPNames, "\\."), function(x) {x[1]})))
}

###
### Gene Ontology (GO) analysis
###
library(goseq)
library(GO.db)

### high level function for GO analysis
GOstats <- function(sampleID, geneNames, negativeGeneNames=NULL, GOSpecies = "hg18", doStatPlots=FALSE) {
	GOtable <- getGOtable(sampleID, geneNames, negativeGeneNames)
	pwf <- getGOpwf(GOtable, GOspecies = GOSpecies, GOgeneID = "ensGene", sampleID, plotName=paste(sampleID, ".GO.geneLengthBiasPWF.ps", sep=""))
	GO.wall=goseq(pwf, GOSpecies, "ensGene")

	if(doStatPlots) {
		GO.samp=goseq(pwf, GOSpecies, "ensGene", method="Sampling", repcnt=1000)
		postscript(paste(sampleID,".GO.WallVsSampling.ps",sep=""))
		plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.samp[,1],GO.wall[,1]),2]), xlab="log10(Wallenius p-values)", ylab="log10(Sampling p-values)", xlim=c(-3,0))
		abline(0,1,col=3,lty=2)
		dev.off()
	}

	cat(paste("GO categories overrepresented in sample ",sampleID,":\n",sep=""))
	cat(getGOOverRepresented(GO.wall, pvalthr=0.05))
	cat(paste("GO categories underrepresented in sample ",sampleID,":\n",sep=""))
	cat(getGOUnderRepresented(GO.wall, pvalthr=0.05))

	GO.wall
}

getGOtable <- function(sampID, geneNames, negativeGeneNames) {
	allCand <- sapply(as.character(geneNames), function(x) { if(is.na(x)) "NA" else strsplit(x,"\\.")[[1]][1] })
	posSet <- unique(allCand)
	negSet <- negativeGeneNames
	GOtable <- c(rep(1,length(posSet)), rep(0, length(negSet)))
	names(GOtable) <- c(as.character(posSet), as.character(negSet))
	GOtable <- GOtable[which(!is.na(names(GOtable)))]
	GOtable
}

getGOpwf <- function(GOtable, GOspecies = "hg19", GOgeneID = "ensGene", sampID, plotName) {
	pwf=nullp(GOtable,GOspecies,GOgeneID) ;# knownGene = UCSC gene names
	file.rename("Rplots.pdf",paste(sampID, ".GO.geneLengthBiasPWF.ps", sep=""))
	pwf
}

getGOOverRepresented <- function(GO.wall, pvalthr = 0.05) {
	overrep <- c()
	GO.overrep = GO.wall[which(GO.wall$upval < pvalthr),]$category
	for(go in GO.overrep){
		# print(GOTERM[[go]]) ; cat("--------------------------------------\n")
		overrep <- c(overrep, paste(Ontology(GOTERM[[go]]),":",Term(GOTERM[[go]]),sep=""))
	}
	overrep
}

getGOUnderRepresented <- function(GO.wall, pvalthr = 0.2) {
	underrep <- c()
	GO.underrep = GO.wall[which(GO.wall$dpval < pvalthr),]$category
	for(go in GO.underrep){
		# print(GOTERM[[go]]) ; cat("--------------------------------------\n")
		underrep <- c(underrep, paste(Ontology(GOTERM[[go]]),":",Term(GOTERM[[go]]),sep=""))
	}
	underrep
}

###
### epimutation calling section
###
#######################################################
### Algorithm for finding epimutations
#######################################################
# 1. Compare lane1, lane2, livy1, 99% conversion rate
# 2. Find overlapping CpG sites: intersect
# 3. Compare lane1 vs livy1 and lane2 vs livy:
# 4. List of common CpG sites > 90% methylated := fullyMethylated
# 5. List of common CpG sites < 10% methylated := fullyUnmethylated
# 6. setdiff(fullyMethylated, fullyUnmethylated in livy1) := methylatingEventSites
#    setdiff(a,b) present in a but not b
# 7. setdiff(fullyUnmethylated, fullyMethylated in livy1) := demethylatingEventSites

getDemethylatingEpimutationSites2 <- function(SingleCell, Control) {
	STATE_METH = 1
	STATE_UNMETH = (-1)

	commonPositions <- intersect(SingleCell$site, Control$site)
	commonInSingle <- SingleCell[which(SingleCell$site %in% commonPositions),]
	commonInControl <- Control[which(Control$site %in% commonPositions),] ;# same sites as above, but counts for refrence

	# Silvia: "we look for the ones that are fully methylated in the CTR (>90%)"
	Control_FullMeth <- commonInControl[which(commonInControl$methstate == STATE_METH),] ;# fully methylated sites in reference, common to lane1

	# "then we go back to the single cell and we ask whether those are fully methylated in the single cell"
	commonInSingle <- commonInSingle[which(commonInSingle$site %in% Control_FullMeth$site),]

	# "If not, we ask whether methylation is <10% in the single cell. If so, we call this epimutation."
	SingleCell_SitesFullyMethInControl_FullyUnmethylatedInSingleCell <- commonInSingle[which(commonInSingle$methstate == STATE_UNMETH),]

	DemethylatingEpimutations <- SingleCell_SitesFullyMethInControl_FullyUnmethylatedInSingleCell
	return(DemethylatingEpimutations)
}

# analog algorithm to demethylation search, see above
getMethylatingEpimutationSites2 <- function(SingleCell, Control) {
	STATE_METH = 1
	STATE_UNMETH = (-1)

	commonPositions <- intersect(SingleCell$site, Control$site)
	commonInSingle <- SingleCell[which(SingleCell$site %in% commonPositions),]
	commonInControl <- Control[which(Control$site %in% commonPositions),] ;# same sites as above, but counts for refrence

	Control_FullUnmeth <- commonInControl[which(commonInControl$methstate == STATE_UNMETH),]
	commonInSingle <- commonInSingle[which(commonInSingle$site %in% Control_FullUnmeth$site),]
	SingleCell_SitesFullyUnmethInControl_FullyMethylatedInSingleCell <- commonInSingle[which(commonInSingle$methstate == STATE_METH),]

	MethylatingEpimutations <- SingleCell_SitesFullyUnmethInControl_FullyMethylatedInSingleCell
	return(MethylatingEpimutations)
}

### function works for both methylation and demethylation rate
getEpimutationRate <- function(SingleCell, Control, eventSites, calcMethylatingEvents) {
	commonPositions <- intersect(SingleCell$site, Control$site)
	commonInControl <- Control[which(Control$site %in% commonPositions),]

	if(calcMethylatingEvents == TRUE) {
		## "number of methylating events in the single cell divided by total number of fully demethylated cytosines in the liver ctr (for those CpG sites covered in both samples)"
		#Control_FullUnmeth <- filterLowMethylation(commonInControl)
		#epimutRateOLD <- dim(eventSites)[1] / dim(Control_FullUnmeth)[1]
		epimutRate <- dim(eventSites)[1] / length(commonPositions) ;# methylating epimutation rate = # methylating events / # total shared CG sites
	} else {
		## "number of demethylating events in the single cell divided by total number of fully methylated cytosines in the liver ctr (for those CpG sites covered in both samples)"
		#Control_FullMeth <- filterHighMethylation(commonInControl) ;# fully methylated sites in reference, common to lane1
		#epimutRateOLD <- dim(eventSites)[1] / dim(Control_FullMeth)[1]		
		epimutRate <- dim(eventSites)[1] / length(commonPositions) ;# demethylating epimutation rate = # demethylating events / # total shared CG sites
	}

	return (round(epimutRate,5))
}

printEpimutStats <- function(name, SingleCell, ControlCells, methEvents, demethEvents, minReads, verboseChromosomeStats=FALSE) {
	STATE_METH = 1
	STATE_UNMETH = (-1)

	commonPositions <- intersect(SingleCell$site, ControlCells$site)
        commonInControl <- ControlCells[which(ControlCells$site %in% commonPositions),]
	commonInSingle <- SingleCell[which(SingleCell$site %in% commonPositions),]
	Control_FullMeth <- commonInControl[which(commonInControl$methstate == STATE_METH),]
	Control_FullUnmeth <- commonInControl[which(commonInControl$methstate == STATE_UNMETH),]

	methRate <- getEpimutationRate(SingleCell, ControlCells, eventSites=methEvents, calcMethylatingEvents=TRUE)
	demethRate <- getEpimutationRate(SingleCell, ControlCells, eventSites=demethEvents, calcMethylatingEvents=FALSE)
	totalRate <- sum(c(methRate,demethRate))

	cat(paste("\n\nStatistics for sample: ", name, " (min. coverage: ", minReads, " reads/site)\n============================================================\n",sep=""))
	cat(paste("Total shared CG sites(=regions!) between ", name, " and reference/CTRL: ", length(commonPositions), "\n",sep=""))
	cat(paste("Median CG sites of these total shared regions (REFERENCE): ", median(commonInControl$meth + commonInControl$unmeth), "\n",sep=""))
	cat(paste("Median CG sites of these total shared regions (SINGLE CELL): ", median(commonInSingle$meth + commonInSingle$unmeth), "\n",sep=""))
	cat(paste("Control has ", dim(Control_FullMeth)[1], " fully methylated and ", dim(Control_FullUnmeth)[1], " fully unmethylated sites\n"),sep="")
	cat(paste("Methylating Epimutation Rate: ", methRate, "\nDemethylating Epimutation Rate: ", demethRate, "\nTotal Epimutation Rate: ", totalRate, "\n",sep=""))

	### For Table 2 stats (methods paper):
	commonInSingle <- SingleCell[which(SingleCell$site %in% commonPositions),]
	cisMeth <- commonInSingle[which(commonInSingle$methstate == STATE_METH),]
	cisUnmeth <- commonInSingle[which(commonInSingle$methstate == STATE_UNMETH),]
	cat(paste("Fully meth in single: ", length(unique(cisMeth$site)), ", Fully unmeth in single: ", length(unique(cisUnmeth$site)), "\n",sep=""))

	if(verboseChromosomeStats) {
		cat(paste("Sites fully methylated in control are on chromosomes: ", paste(unique(Control_FullMeth$chr), collapse=","), "\n",sep=""))
		cat(paste("Sites fully unmethylated in control are on chromosomes: ", paste(unique(Control_FullUnmeth$chr), collapse=","), "\n",sep=""))
		cat(paste("Meth. epimutations are on chromosomes: ", paste(unique(methEvents$chr),collapse=","), "\n",sep=""))
		cat(paste("Demeth. epimutations are on chromosomes: ", paste(unique(demethEvents$chr),collapse=","), "\n",sep=""))
	}
	c(methRate, demethRate)
}

### assumes that verbose information (e.g. annotation) is in cpg1
### result = cpg1 sites matching ranges of cpg2
cpgIntersection <- function(cpg1, cpg2) {
	cpg1 <- GRanges(seqnames=cpg1$chrom, ranges=IRanges(start=cpg1$pos, end=cpg1$pos))
	cpg2 <- GRanges(seqnames=cpg2$chrom, ranges=IRanges(start=cpg2$pos, end=cpg2$pos))
}

###
### PARAMETER HANDLING SECTION
### 

makeParams <- function(localpath, sampNames, convrates, is.reference, pminus=0.2, regionSize=10000, minCounts=5, verbose=TRUE, computeRegions=FALSE, computeMatrices=TRUE, writeEpicallMatrix=TRUE, doPlots=FALSE) {
	if(length(pminus) == 1) { pminus = rep(pminus, length(sampNames)) }

	pplus = 1 - convrates ;# !!!

	stopifnot(length(sampNames) == length(pminus))
	stopifnot(length(pminus) == length(pplus))
	stopifnot(length(pplus) == length(is.reference))

	localpath = paste(localpath, "/", sep="")

	params <- list(localpath, sampNames, pplus, is.reference, pminus, regionSize, minCounts, verbose, computeRegions, computeMatrices, writeEpicallMatrix, doPlots)
	names(params) <- c("localpath", "sampNames", "pplus", "is.reference", "pminus", "regionSize", "minCounts", "verbose", "computeRegions", "computeMatrices", "writeEpicallMatrix", "doPlots")
	names(params$pplus) <- sampNames
	names(params$pminus) <- sampNames
	names(params$is.reference) <- sampNames

	params
}


###
### HIGH-LEVEL API SECTION
###

### PART 1: POSITIONS TO REGIONS
positions_to_regions <- function(params) {
	regionsList <- list()
	localpath = params[["localpath"]]
	sampNames = params[["sampNames"]]
	regionSize = params[["regionSize"]]
	minCounts = params[["minCounts"]]
	verbose = params[["verbose"]]

	lanesCpg <- paste(localpath, sampNames, ".positions.csv",sep="")
	lanesRegions <- paste(localpath, sampNames, ".regions.", regionSize, ".", minCounts, ".RData",sep="")

	if(params[["computeRegions"]]) {
		for(i in 1:length(lanesCpg)) {
			if(verbose) cat(paste("Sample: ", basename(lanesCpg[i]), " regionSize: ", regionSize, " minCounts: ", minCounts, "\n",sep=""))
			#load(lanesCpg[i])
			positions <- positionsReadCSV(lanesCpg[i])
			cpgregions <- cpg2regions(positions, windowSize=regionSize, minCountsPerRegion=minCounts)
			if(verbose) cat(paste("Processed ", basename(lanesCpg[i]), ", yielding ", nrow(cpgregions), " regions of ", regionSize, " nt\n",sep=""))
			rm(positions)
			regionsList[[lanesCpg[i]]] <- cpgregions
			save(cpgregions, file=lanesRegions[i])
		}
	} else {
		for(i in 1:length(lanesCpg)) {
			load(lanesRegions[i])
			regionsList[[lanesCpg[i]]] <- cpgregions
		}
	}
	regionsList
}

### PART 2: METHYLATION STATUS CORRECTION
generate_results <- function(params) {
	localpath = params[["localpath"]]
	sampNames = params[["sampNames"]]
	isReferenceVector = params[["is.reference"]]
	pplus = params[["pplus"]]
	pminus = params[["pminus"]]
	computeMatrices = params[["computeMatrices"]]
	verbose = params[["verbose"]]
	regionSize = params[["regionSize"]]
	minCounts = params[["minCounts"]]

	ret <- list()

	for(sampName in sampNames) {
		cpgregion <- MethCallGetCorrectedRegion(sampName, isReferenceVector, pplus[sampName], pminus[sampName], localpath, precompute=computeMatrices, regionSize, minCounts, verbose)
		# resulting region will also be saved to: paste(localpath,a,".results.", regionSize, ".", minCounts, ".RData",sep="")
		ret[[sampName]] <- cpgregion
	}
	ret
}

### PART 3: EPIMUTATION CALLING
epimutation_calls <- function(params) {
	localpath = params[["localpath"]]
	sampNames = params[["sampNames"]]
	isReferenceVector = params[["is.reference"]]

	pminus = params[["pminus"]]
	pplus = params[["pplus"]]
	computeMatrices = params[["computeMatrices"]]
	verbose = params[["verbose"]]
	doPlots = params[["doPlots"]]

	regionSize = params[["regionSize"]]
	minCounts = params[["minCounts"]]
	writeEpicallMatrix = params[["writeEpicallMatrix"]]

	referenceName = sampNames[isReferenceVector]
	stopifnot(length(sampNames[isReferenceVector]) <= 1) ;# analysis only works with a single reference

	### filename for ref_results
	referenceFile = paste(localpath, referenceName, ".results.", regionSize, ".", minCounts, ".RData",sep="")

	methRatesList <- list()
	demethRatesList <- list()
	methSites <- list()
	demethSites <- list()

	singleCellSamples <- sampNames[!isReferenceVector]

	for(sampName in singleCellSamples) {
		### filename for single_results
		corrFile = paste(localpath, sampName, ".results.", regionSize, ".", minCounts, ".RData",sep="")

		ret = epimutationRegionStats(singleCellFile=corrFile, referenceFile, singleName=sampName, minReads=minCounts, human=F, windowRange=regionSize)

		### Rates and ANNOTATION-SPEIFIC rates
		methRatesList[[sampName]] = ret[[3]]
		demethRatesList[[sampName]] = ret[[4]]
		methSites[[sampName]] = ret[[1]][,1:6]
		demethSites[[sampName]] = ret[[2]][,1:6]

		if(writeEpicallMatrix) {
			fileMeth = paste(localpath, sampName, ".methEpicalls.",regionSize, ".", minCounts,"p+=",pplus[sampName],"p-=",pminus[sampName],".RData", sep="")
			meth = methSites[[sampName]][,c("chr","pos","endpos","meth","unmeth","methstate"),]
			names(meth) = c("chr","start","stop","meth","unmeth","epimutation_call_test")
			meth <- addRegionID(meth)
			#meth <- merge(meth,vitMerge,by="id")
			save(meth,file=fileMeth)
			if(verbose) cat(paste("Written epi-methylation calls to matrix: \n", basename(fileMeth), "\n", sep=""))

			fileDemeth = paste(localpath, sampName, ".demethEpicalls.",regionSize, ".", minCounts,"p+=",pplus[sampName],"p-=",pminus[sampName],".RData", sep="")
			demeth = demethSites[[sampName]][,c("chr","pos","endpos","meth","unmeth","methstate"),]
			names(demeth) = c("chr","start","stop","meth","unmeth","epimutation_call_test")
			demeth <- addRegionID(demeth)
			#demeth <- merge(demeth,vitMerge,by="id")
			save(demeth,file=fileDemeth)
			if(verbose) cat(paste("Written epi-demethylation calls to matrix:\n", basename(fileDemeth), "\n", sep=""))
		}
	}

	### dev version
	### ret <- list(methRatesList, demethRatesList, methSites, demethSites, pminus, pplus)
	### names(ret) <- c("methRates","demethRates","methSites","demethSites","pminus","pplus")

	### stable version
	ret <- list(methSites, demethSites)
	names(ret) <- c("methSites","demethSites")

	ret
}

epimutationRegionStats <- function(singleCellFile, referenceFile, singleName, minReads, human=F, corr=NULL, corrRef=NULL, windowRange) {
	load(singleCellFile) ; single <- cpgregions ; rm(cpgregions) ; names(single) <- c("chr","pos","endpos","meth","unmeth","methstate") ; single <- addChrPosMeth(single)
	load(referenceFile) ; ref <- cpgregions ; rm(cpgregions) ; names(ref) <- c("chr", "pos","endpos","meth","unmeth","methstate") ; ref <- addChrPosMeth(ref)

	### Remove NA posterior values from methylation level estimation
	### These occur only in very high peaks (PCR artifacts) of abs(meth - unmeth) > 1000
	### On average, this removes about 0.000008 biased positions
	single <- single[which(!is.na(single$methLevel)),]
	ref <- ref[which(!is.na(ref$methLevel)),]

	if(minReads > 1) {
		single <- filterMinCoverage(single, minReads)
		ref <- filterMinCoverage(ref, minReads)
	}

	methEventSites <- getMethylatingEpimutationSites2(single, ref)
	demethEventSites <- getDemethylatingEpimutationSites2(single, ref)
	overlaplen <- length(intersect(single$site, ref$site))

	commonSites <- ref[which(ref$site %in% intersect(single$site, ref$site)),]

	methDemethRates <- printEpimutStats(singleName, single, ref, methEventSites, demethEventSites, minReads)
	totalMethRate <- methDemethRates[1]
	totalDemethRate <- methDemethRates[2]

	list(methEventSites,demethEventSites,totalMethRate,totalDemethRate)
}

compareHMM <- function(probs1File, probs2File) {
	load(probs1File) ; probs1 <- unlist(viterbi) ; rm(viterbi)
	load(probs2File) ; probs2 <- unlist(viterbi) ; rm(viterbi)
	rm(hmm)

	commonIDs <- intersect(names(probs1),names(probs2))
	probs1comm <- probs1[which(names(probs1) %in% commonIDs)]
	probs2comm <- probs2[which(names(probs2) %in% commonIDs)]

	stopifnot(length(probs1comm) == length(probs2comm))

	res = table(probs1comm,probs2comm)
	rownames(res) = c("lik_epimeth","lik_epidemeth","lik_ambig")
	colnames(res) = rownames(res)
	res
}

###
### MODEL SECTION
###

############################################################################################
###
### Methylation calling/estimation section
###
############################################################################################


### bsconv = convrates[sampleName]
### pminus = 0.1 by default
### precompute = TRUE means static matrix precomputation
MethCallGetCorrectedRegion <- function(sampleName, isReferenceVector, pplus, pminus, localpath, precompute=T, regionSize, minCounts, verbose=T, nmaxStatic=50) {
	a <- sampleName

	if(verbose) cat(paste("Sample: ",a,"\n",sep=""))

	if(isReferenceVector[a])
	{
		multiOrSingle <- "strict_call"
		multiOrSingleMat <- "strict_call_matrix"
	} else {
		multiOrSingle <- "permissive_call"
		multiOrSingleMat <- "permissive_call_matrix"
	}

	### Begin hardcoded parameters by Achim
	# r contains the points at which the posterior distribution will be calculated
	steps = 100
	r = seq(1/steps,1-1/steps,length=steps)
	# parameters of the beta mixture prior
	evi = 0.5 # evidence weight ("peakedness") for both beta distributions
	ameth = 9*evi; bmeth = 1*evi # shape of the methylated beta distribution
	aunmeth = 1*evi; bunmeth = 9*evi # shape of the unmethylated beta distribution
	lambda = 0.7 # mixing parameter, relative weight of the methylated beta distr.
	# parameter for the confidence level at which the methylation calls shall be made
	confidence_level = 0.75 ;# sep 2012 changed from 0.9
	### End hardcoded parameters

	if(precompute) {
		### depends only on p+ and p- (i.e. BS-conversion rate), not region size!
		if(verbose) cat("Generating conversion matrix...\n")
		convmat = generate.table(nmaxStatic,pminus,pplus,r,lambda,ameth,bmeth,aunmeth,bunmeth,confidence_level)
		save(convmat, file=paste(localpath,a,".convMat.",pminus,".",pplus,".RData",sep="")) ;# convMat is generic for all region sizes
	}

	load(paste(localpath, a, ".regions.", regionSize, ".", minCounts, ".RData",sep=""))
	load(paste(localpath, a, ".convMat.",pminus,".",pplus,".RData",sep=""))

	### dynamic calculation rates for new parameters (0.3/0.7/confidence 0.75)
	if(isReferenceVector[a]) {
		x <- (convmat[["strict_call_matrix"]][,nmaxStatic]) ;# FOR REFERENCE
	} else {
		x <- (convmat[["permissive_call_matrix"]][,nmaxStatic]) ;# FOR SINGLE-CELL
	}

	unmethCutoff <- round(max(which(x == -1))/length(x),4)
	methCutoff <- round(min(which(x == 1))/length(x),4)

	methUnmeth <- cpgregions[,4:5]
	methUnmeth[,1] <- as.integer(methUnmeth[,1])
	methUnmeth[,2] <- as.integer(methUnmeth[,2])

	results <- apply(methUnmeth, 1, function(x) {
		k <- x[1]
		n <- k+x[2]

		if((k <= nmaxStatic) && (n <= nmaxStatic)) {
			### Achim: zeile k+1 steht fuer Wert k (zeile 1 fuer k=0)
			stat <- convmat[[multiOrSingleMat]][k+1,n]
			return(stat)
		} else {
			#cat(".");
			#res=statistiken(k,n,pminus,pplus=1-bsconv,lambda,ameth,bmeth,aunmeth,bunmeth,confidence_level)
			#res[[multiOrSingle]]
			methLevel <- k/n
			if(methLevel <= unmethCutoff) {
				#cat(paste("DEBUG 1 k: ",k," n: ",n,"\n",sep=""))
				#fullMat[[multiOrSingleMat]][k+1,n] <- (-1)
				#cat(paste("DEBUG: ",fullMat[[multiOrSingleMat]][k+1,n],"\n",sep=""))
				stat <- -1
				return(stat)
			} else if(methLevel >= methCutoff) {
				#fullMat[[multiOrSingleMat]][k+1,n] <- 1
				stat <- 1
				return(stat)
			} else {
				#fullMat[[multiOrSingleMat]][k+1,n] <- 0
				stat <- 0
				return(stat)
			}
		}
		stopifnot(FALSE) ;# should return before
	})
	### fix numeric(0) for empty regions
	results <- lapply(results, function(x) { if(length(x) == 0) 0 else x })

	resultsMethEst <- apply(methUnmeth, 1, function(x) {
		k <- x[1]
		n <- k+x[2]
		if((k <= nmaxStatic) && (n <= nmaxStatic)) {
			### Achim: zeile k+1 steht fuer Wert k (zeile 1 fuer k=0)
			methestimate <- convmat[["methest_matrix"]][k+1,n]
			return(methestimate)
		} else {
			### ??? Achim vorher: methylation estimate = round(k/nmaxStatic)+1 ???
			methestimate = ((k/n) - pplus) / (1 - pminus - pplus)
			#fullMat[["methest_matrix"]][k+1,n] <- methestimate
			return(methestimate)
		}
		stopifnot(FALSE)
	})
	### fix numeric(0) for empty regions
	resultsMethEst <- lapply(resultsMethEst, function(x) { if(length(x) == 0) 0 else x })

	if(verbose) cat("Adding methylation states and 9 matrices...\n")
	stopifnot(length(results) == nrow(cpgregions))
	stopifnot(length(unlist(results)) == nrow(cpgregions))

	### for only non-null output
	#non_null = which((cpgregions$meth + cpgregions$unmeth) > 0)
	#cpgregions <- cpgregions[non_null,]
	#cpgregions1 <- cbind(cpgregions,methstate=unlist(results[non_null]))

	cpgregions1 <- cbind(cpgregions,methstate=unlist(results))
	cpgregions1 <- cpgregions_to_results(cpgregions1,convmat)
	cpgregions <- cpgregions1

	#stopifnot(length(resultsMethEst) == nrow(cpgregions))
	#stopifnot(length(unlist(resultsMethEst)) == nrow(cpgregions))
	#cpgregions <- cbind(cpgregions,methest=resultsMethEst)
	#cpgregions <- cpgregions_to_results(cpgregions, convmat)

	save(cpgregions,file=paste(localpath,a,".results.", regionSize, ".", minCounts, ".RData",sep=""))
	## first additional column: meth = 1, unmeth = -1 // further columns: 9 matrices
	cpgregions
}

## calculation of one term in the sum p(k|n,r,p+,p-) = sum_j sum_m { koeff(j,m,k,n,r,p+,p-) }

koeff = function(j,m,k,n,r,pminus,pplus){
	if ((m > j) | (k-m > n-j) | (k-m < 0)) stop("Illegal (j,m,k,n) indices.") 
	res = dbinom(x=m,size=j,prob=1-pminus) * 
		dbinom(x=k-m,size=n-j,prob=pplus) * 
		dbinom(x=j,size=n,prob=r)
	return(res)
}



## calculation of the likelihood p(r|k,n,p-,p+) for r a given grid of methylation rates

likr = function(r,k,n,pminus,pplus){
	j = as.vector(outer(0:n,0:n,function(x,y){x}))
	m = as.vector(outer(0:n,0:n,function(x,y){y}))
	sel = which(!((m > j) | (k-m > n-j) | (k-m<0)))
	j = j[sel]
	m = m[sel]
	res = apply(cbind(j,m),1,function(x){koeff(x[1],x[2],k,n,r,pminus,pplus)})
	res = rowSums(res)
	return(res)
}



## calculation of the prior p(r|...)= 
##  = lambda*Beta(r;ameth,bmeth) + (1-lambda)*Beta(r;aunmeth,bunmeth)
## for r a given grid of methylation rates

priorr = function(r,lambda,ameth,bmeth,aunmeth,bunmeth){
	res = lambda * dbeta(r,ameth,bmeth) + (1-lambda) * dbeta(r,aunmeth,bunmeth)
	return(res)
}



## calculates the posterior distribution of the methylation rates
## the vector r contains the x-grid points at which the 
## posterior distribution p(r|k,n,p+,p-) is to be evaluated
## additional parameters are: pminus (the rate of false non-methylation calls),
## pplus (the rate of false methylation calls), and the parameters speciying the
## prior, which is a 2-Beta mixture: lambda*Beta(ameth,bmeth) + (1-lambda)*Beta(aunmeth,bunmeth)
## ameth, bmeth, aunmeth, bunmeth
## the functiion outputs the density values of the posterior at the given grid points r
## (later, these density values are commonly denoted by pr)

postr = function(r,k,n,pminus,pplus,lambda,ameth,bmeth,aunmeth,bunmeth){
	res = likr(r,k,n,pminus,pplus) * priorr(r,lambda,ameth,bmeth,aunmeth,bunmeth)
	res = res/sum(res)
	return(res)
}



## auswertung expects a vector r containing the x-grid values at which the density
## of the posterior p(r|k,n,p+,p-) is evaluated, and pr contains the corresponding 
## values of p(r|k,n,p+,p-). confidence_level contains the probability for the 
## methylation levels being above/below the methl/nonmeth level thresholds,
## sum(pr[r>methlevel]) resp. sum(pr[r<nonmethlevel]))
 
auswertung = function(r,pr,confidence_level=0.9,methlevel=0.7,nonmethlevel=0.3){

	# methestimate is the expectation value of the posterior distribution pr
	methestimate = sum(r*pr) 

	# right and left are the left and right boundaries of the <confidendce_level>
	# credible intervals around the expectation value (methestimate)
	cumdist = cumsum(pr)
	rightprob = min(methestimate+confidence_level/2,1) + (-min(methestimate-confidence_level/2,0))
	leftprob  = max(methestimate-confidence_level/2,0) - (max(methestimate+confidence_level/2,1) - 1)
	right = r[which.min(abs(cumdist-rightprob))]
	left  = r[which.min(abs(cumdist-leftprob))]

	# p_strict_meth and p_strict_unmeth are the probabilities supporting the highly meth/unmeth calls
	p_strict_unmeth = sum(pr[r<nonmethlevel]) ### sep 2012 changed from 0.2 / 0.8
	p_strict_meth = sum(pr[r>methlevel])
	# the strict call is made if pmeth resp. punmeth are above the <confidence_level>
	# with values strict_call =1 for highly meth, =0 for no strict call, =-1 for sparsely meth
	if (p_strict_meth > confidence_level) {
		strict_call = 1
	} else if (p_strict_unmeth > confidence_level) {
		strict_call = -1
	} else {
		strict_call = 0
	}

	# p_perm_meth and p_perm_unmeth are the probabilities supporting the increased meth/unmeth calls
	p_perm_unmeth = sum(pr[r<0.5])
	p_perm_meth = sum(pr[r>0.5])
	# the permissive call is made if epimeth resp. epiunmeth are above the <confidence_level>
	# with values permissive_call =1 for increased meth, =0 for no call, =-1 for decreased meth
	if (p_perm_meth>confidence_level) {
		permissive_call = 1
	} else if (p_perm_unmeth>confidence_level) {
		permissive_call = -1
	} else {
		permissive_call = 0
	}

	return(list(methestimate=methestimate,left=left,right=right,strict_call=strict_call,permissive_call=permissive_call,
			p_strict_meth=p_strict_meth,p_strict_unmeth=p_strict_unmeth,
			p_perm_unmeth=p_perm_unmeth,p_perm_meth=p_perm_meth))
}



## statistiken generates all interesting statistics for a given region with k methylation calls
## among n calls. 
## the input parameters (pminus,pplus,lambda,ameth,bmeth,aunmeth,bunmeth) specify the likelihood and the prior
## the additional parameters (confidence_level,methlevel,nonmethlevel define the methylation calls
## and r is a grid of values between 0 and 1 at which the posterior will be evaluated (should contain >=100 points)
statistiken = function(k,n,pminus,pplus,r,
				lambda=0.5,ameth=0.7*0.5,bmeth=(1-0.7)*0.5,aunmeth=0.3*0.5,bunmeth=(1-0.3)*0.5,
				confidence_level=0.9,methlevel=0.7,nonmethlevel=0.3){
	### 3 lines added by Achim 3/10/13
	if (n==0) return(list(methestimate=methestimate,left=left,right=right,strict_call=strict_call,permissive_call=permissive_call,
		p_strict_meth=p_strict_meth, p_strict_unmeth=p_strict_unmeth,
		p_perm_unmeth=p_perm_unmeth, p_perm_meth=p_perm_meth) )
	reslikr = likr(r,k,n,pminus,pplus)
	respriorr = priorr(r,lambda,ameth,bmeth,aunmeth,bunmeth)
	posterior = reslikr * respriorr
	pr = posterior/sum(posterior)
	res = auswertung(r,pr,confidence_level)
	return(res)
}



## generate.table
generate.table = function(nmax,pminus,pplus,r,
					lambda,ameth,bmeth,aunmeth,bunmeth,
					confidence_level,methlevel,nonmethlevel){
	methest_matrix 		= matrix(NA,ncol=nmax,nrow=nmax+1)
	lower_bound_matrix 	= matrix(NA,ncol=nmax,nrow=nmax+1)
	upper_bound_matrix 	= matrix(NA,ncol=nmax,nrow=nmax+1)
	p_strict_meth_matrix	= matrix(NA,ncol=nmax,nrow=nmax+1)
	p_strict_unmeth_matrix	= matrix(NA,ncol=nmax,nrow=nmax+1)
	p_perm_unmeth_matrix	= matrix(NA,ncol=nmax,nrow=nmax+1)
	p_perm_meth_matrix	= matrix(NA,ncol=nmax,nrow=nmax+1)
	strict_call_matrix 	= matrix(NA,ncol=nmax,nrow=nmax+1)
	permissive_call_matrix 	= matrix(NA,ncol=nmax,nrow=nmax+1)

	for (n in 1:nmax){
		for (k in 0:n){
			res = statistiken(k,n,pminus,pplus,r,
				lambda,ameth,bmeth,aunmeth,bunmeth,
				confidence_level,methlevel,nonmethlevel)
			methest_matrix[k+1,n] 		= res$methestimate
			lower_bound_matrix[k+1,n] 	= res$left
			upper_bound_matrix[k+1,n] 	= res$right
			p_strict_meth_matrix[k+1,n]	= res$p_strict_meth
			p_strict_unmeth_matrix[k+1,n]	= res$p_strict_unmeth
			p_perm_meth_matrix[k+1,n]	= res$p_perm_meth
			p_perm_unmeth_matrix[k+1,n]	= res$p_perm_unmeth
			strict_call_matrix[k+1,n] 	= res$strict_call
			permissive_call_matrix[k+1,n] = res$permissive_call
		}
	}

	return(list(methest_matrix=methest_matrix, 
			lower_bound_matrix=lower_bound_matrix, upper_bound_matrix=upper_bound_matrix,
			p_strict_meth_matrix=p_strict_meth_matrix, p_strict_unmeth_matrix=p_strict_unmeth_matrix,
			p_perm_meth_matrix=p_perm_meth_matrix, p_perm_unmeth_matrix=p_perm_unmeth_matrix,
			strict_call_matrix=strict_call_matrix,
			permissive_call_matrix=permissive_call_matrix))
}

## very useful function (why does this not exist in R?)
## function by Achim, version from 3/11/2013
matrixreadout = Vectorize(function(x,y,mat){
    if (y==0) return(NA)
    if (y <= ncol(mat)) return(mat[x,y])
    return(mat[round((x-1)/y*ncol(mat))+1,ncol(mat)])
    },vectorize.args = c("x","y"))

addRegionID <- function(cpgregions) {
	id = paste(cpgregions$chr, cpgregions$start, cpgregions$stop, sep=".")
	cbind(cpgregions, id=id, stringsAsFactors=FALSE)
}

### eliminate <= 0.01% of empty chromosome-end regions
commonRegions <- function(ref_results, single_results) {
	single_results <- addRegionID(single_results)
	ref_results <- addRegionID(ref_results)
	
	commonRegions <- intersect(single_results$id, ref_results$id)
	single_results <- single_results[which(single_results$id %in% commonRegions),]
	ref_results <- ref_results[which(ref_results$id %in% commonRegions),]

	stopifnot(single_results$id == ref_results$id)

	list(ref_results, single_results)
}

## expands the cpgregions matrix by the entries extracted from convmat
cpgregions_to_results = function(cpgregions,convmat) {
	kplusone = cpgregions[,"meth"] +1
	n = cpgregions[,"unmeth"] + cpgregions[,"meth"]

	methest = matrixreadout(kplusone,n,convmat$methest_matrix)
	lower_bound = matrixreadout(kplusone,n,convmat$lower_bound_matrix)
	upper_bound = matrixreadout(kplusone,n,convmat$upper_bound_matrix)
	p_strict_meth = matrixreadout(kplusone,n,convmat$p_strict_meth_matrix)
	p_strict_unmeth = matrixreadout(kplusone,n,convmat$p_strict_unmeth_matrix)
	p_perm_meth = matrixreadout(kplusone,n,convmat$p_perm_meth_matrix)
	p_perm_unmeth = matrixreadout(kplusone,n,convmat$p_perm_unmeth_matrix)
	strict_call = matrixreadout(kplusone,n,convmat$strict_call_matrix)
	permissive_call = matrixreadout(kplusone,n,convmat$permissive_call_matrix)

	cpgregions = cbind(cpgregions,methest,lower_bound,upper_bound,
		p_strict_meth,p_strict_unmeth,p_perm_meth,p_perm_unmeth,
		strict_call,permissive_call)

	return(cpgregions)
}

### 2 following functions added by Achim 10/3/13
## This function returns a 3-column matrix, each row corresponds to one region,
## the columns represent the probabilities for 
## an epimethylation event, no event, an apidemethylation event
calculate_emission_probabilities = function(ref_results,single_results) {
	rs <- commonRegions(ref_results, single_results)
	ref_results <- rs[[1]]
	single_results <- rs[[2]]

	# calculate P(event|k_ref,n_ref,k_single,n_single) probabilities
	p_epimeth = ref_results[,"p_strict_unmeth"] * single_results[,"p_perm_meth"]
	p_epidemeth = ref_results[,"p_strict_meth"] * single_results[,"p_perm_unmeth"]
	p_ambiguous = 1 - p_epimeth - p_epidemeth

	# calculate P(k_ref,n_ref,k_single,n_single|event) "likelihoods"
	sum_epimeth = sum(p_epimeth,na.rm=T)
	sum_epidemeth = sum(p_epidemeth,na.rm=T)
	sum_ambiguous = sum(p_ambiguous,na.rm=T)
	sumsum = sum_epimeth + sum_epidemeth + sum_ambiguous

	k = ref_results[,"meth"]

	### n is TOTAL counts here (meth+unmeth) -- 3/11/2013
	n = k + ref_results[,"unmeth"]

	kn_factor = factor(paste(k,"_",n,sep=""))
	kn_counts = table(kn_factor)
	kn_factor = kn_counts[kn_factor] 
	# kn_factor contains the number of regions with (k meth, n total) counts
	# likelihoods are all proportionally rescaled by sumsum
	# to make them live on a nicer scale
	lik_epimeth = p_epimeth*kn_factor/sum_epimeth * sumsum
	lik_epidemeth = p_epidemeth*kn_factor/sum_epidemeth * sumsum
	lik_ambiguous = p_ambiguous*kn_factor/sum_ambiguous * sumsum

	stopifnot(ref_results$id == single_results$id)
	probs = cbind(lik_epimeth,lik_epidemeth,lik_ambiguous,single_results$id)
	nonNA <- which(!(is.na(probs[,1]) & is.na(probs[,2]) & is.na(probs[,3])))
	probs = probs[nonNA,]
	colnames(probs) = c("lik_epimeth","lik_epidemeth","lik_ambiguous","id")
	return(probs)
}

