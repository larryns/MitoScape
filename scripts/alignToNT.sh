#!/bin/bash

# Align fastq files to just the nuclear genome (remove the mitochondrial 
# chromosome). Use hg38 and newer, hg19 has the wrong mitochondrial chromosome.

[ $# -eq 3 ] || { echo "Usage: $0 <input bam file> <number of threads> <sorting memory>"; exit -1 }

# Args:
#
# 1 - Input bam file
# 2 - Number of threads
# 3 - Sort Memory

# Note:
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

sample="$1"
r1=${OUTDIR}/${sample}.R1.fq
r2=${OUTDIR}/${sample}.R2.fq
oFile=${OUTDIR}/${sample}_NT.bam

# First convert the bam files to fastq. Note that gsnap doesn't handle 
# interleaved fastq files, so piping isn't an option.
samtools fastq \
	-0 /dev/null -1 ${r1} -2 ${r2} \
	-f 0x3 \
	-N \
	${INDIR}/${sample}_MT.bam 

# Align the fastq files to the rCRS.
gsnap \
    -D ${DATADIR}/hg38_nuconly_gsnap \
    -d hg38_nuconly_gsnap \
    --input-buffer-size=10000 \
    --gunzip \
    --nthreads=${THREADS} \
    --npaths=1 \
    --print-snps \
    -A sam \
    --read-group-id=$1 \
    --read-group-name=MT_$1 \
    --read-group-library=bar \
    --read-group-platform=illumina ${r1} ${r2}| \
    /opt/apps/samtools/current/bin/samtools view -f 0x3 -u -| \
    /opt/apps/samtools/current/bin/samtools sort -o $oFile \
        -m ${SORTMEM} -O BAM - \

# Cleanup
rm ${r1} ${r2}
