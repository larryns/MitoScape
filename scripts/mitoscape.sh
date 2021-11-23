#!/bin/bash

# Helper script to run the MitoScape classification.

[ $# -eq 2 ]] || { echo "Usage: $0 <prefix> <prob>"; exit -1; }

# Takes a set of input bam files with prefix: $1, and generates a bam file
# with just mtDNA reads.

# Helper variables. 

# Helper variable so you don't have to write long variable names.
ROOTDIR=""

# Working/temp directory
WORKDIR=${ROOTDIR}/tmp

# Where the input bams are stored
BAMDIR=${ROOTDIR}/Input

# the MitoScape classifier directory.
MTCLDIR=""

# Where the Model (.RF) file is stored.
MODELDIR=""

# Number of threads to use
THREADS=5

# The prediction probability. Around 0.5 should do. 
PROB=$2

# Shouldn't have to change anything below.

# Sample name
sample=$1

# Perform the classification
java -Xmx20G -jar ${MTCLDIR}/target/scala-2.12/MTClassify.jar \
	--threads ${THREADS} \
	--prob ${PROB} \
	--ld ${MTCLDIR}/src/universal/mitomap.ld \
	--numt ${MTCLDIR}/src/universal/NUMTs_hg38.txt \
	--classifier ${MODELDIR}/MTClassifierModel.RF \
	--prefix ${sample} \
	--out ${sample}_MTDNA_${PROB}.bam

exit 0
