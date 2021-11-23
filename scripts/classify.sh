#!/bin/bash

# Run the classifier with GNU parallel. Only use this for a small number of
# samples and a single machine.

ROOTDIR="/data/singhln/Projects/Schizophrenia/mtDNA"
WORKDIR=${ROOTDIR}/tmp
BAMDIR=${ROOTDIR}/Unclassified
BINDIR=${ROOTDIR}/MitoScape
MTCLDIR="/data/singhln/Projects/MTClassify"
MODELDIR="/data/singhln/Projects/MTClassifier"
THREADS=3
PROB=0.4

ls ${BAMDIR}/*_sort_MT.bam| cut -f1-2 -d'_'|parallel "${BINDIR}/mtclassifier.sh {} ${PROB}"  
exit 0
