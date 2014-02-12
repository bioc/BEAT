#!/usr/bin/python

#
# Script to determine BS-conversion rates from
# Bismark >= 0.7 methylation_extractor output
#
# Example preparation using bismark:
# bismark_v0.10.1/bismark_methylation_extractor -s --merge_non_CpG SAMPLENAME.sam
# cat Non_CpG*.SAMPLENAME.txt | grep -v Bismark > NonCpG.SAMPLENAME.txt
#
# Example usage:
# python bsConversionRate.py NonCpG.SAMPLENAME.txt
#
# 2014 by Kemal Akman <akman@mpipz.mpg.de>
#

import pysam
import csv
import os
import sys

inFile = sys.argv[1];

print >> sys.stdout, "Input from bismark_metyhlation_extractor:\t", inFile, "";

methr = csv.reader(open(inFile, 'rb'), delimiter='\t',quotechar='"')

meth = 0;
unmeth = 0;

### fetch all reads from samfile and look at CIGAR "XM" = methylation call string
for line in methr:
	posID = line[3];
	methChar = line[4];

	if(methChar == "X"):
		meth = meth + 1;
	elif(methChar == "x"):
		unmeth = unmeth + 1;
	elif(methChar == "H"):
		meth = meth + 1;
	elif(methChar == "h"):
		unmeth = unmeth + 1;

print >> sys.stdout, "methylated counts: ", meth, " unmethylated counts: ", unmeth;
print >> sys.stdout, "non-CpG methylation rate: ", float(meth) / float(meth+unmeth);
print >> sys.stdout, "bs-conversion rate: ", float(1) - (float(meth) / float(meth+unmeth));

