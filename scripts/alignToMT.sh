#!/bin/bash

# Align a fastq file to just the mitochondrial chromosome (rCRS).

[ $# -eq 3 ] || { echo "Usage: $0 <fastq prefix> <number of threads> <sorting memory>"; exit -1 }

# Args:
#
# 1 - Fastq prefix (including path), assuming the suffix is .1.fq and .2.fq
# 2 - Number of threads
# 3 - Sort Memory

# Note:
#
# We assume that we have a rCRS.fa file in DATADIR.
#
# We also assume that samtools and gsnap are in the path.
#
# rCRS_gsnap is the name of the indexed rCRS (chrM) mitochondrial chromosome
# indexed and prepared for gsnap. Something like:
# gmap_build -d rCRS_gsnap -c chrM, rCRS.fa 
# should work. 
# 
# We assume that the format of the fastqs are <prefix>.R1.fq, but feel free
# to change as necessary.

# Modify these variables as needed

# Root directory. A helper variable so you don't have to type long names.
ROOTDIR=""

# Working/temp directory.
WORKDIR=""

# Where you store the data.
DATADIR=${ROOTDIR}/Data

# Output directory
OUTDIR=${WORKDIR}/Bams

# Number of threads to use for gsnap
THREADS=$2

# How much memory to use for sorting
SORTMEM=$3

# Input directory for fastq files.
INDIR=""

# You shouldn't need to change anything below this line.

oFile=${OUTDIR}/${1}_MT.bam
fastqPath=${INDIR}/$1

# Align the fastq files to the rCRS.
gsnap \
	-D $DATADIR/rCRS_gsnap \
	-d rCRS_gsnap \
	--input-buffer-size=25000 \
	--gunzip \
	--nthreads=${THREADS} \
	--npaths=1 \
	-A sam \
	--read-group-id=$1 \
	--read-group-name=MT_$1 \
	--read-group-library=bar \
	--read-group-platform=illumina \
	$1.R1.fq $1.R2.fq | \
	samtools view -f 0x3 -u -| \
	samtools sort -n -o $oFile -T ${oFile}.tmp -m ${SORTMEM} -O BAM - 
	
# Get the MD file
samtools \
	-e -Q --output-fmt BAM \
	$oFile ${DATADIR}/rCRS.fa > ${OUTDIR}/${1}_MT_MD.bam

exit 0
