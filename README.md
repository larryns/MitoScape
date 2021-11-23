# MitoScape
A big-data, machine-learning workflow for obtaining mtDNA sequence from NGS data.

## Cavatica Approach
Workflow and additional documenatation available at: : https://cavatica.sbgenomics.com/public/apps#d3b-bixu/app-publisher/mitoscape-wf/, including full instructions on how to run MitoScape using the Cavatica framework. If you require help on running MitoScape on Cavatica, please look at the Cavatica documentation and reach out to the Cavatica team.

## Local cluster/server Approach
If you wish to run MitoScape independently of Cavatica, the following are the steps. You can use the fat jar including here or clone the source and create the fat jar with 
```
sbt assembly
```

The basic workflow is as follows:

### Alignment

MitoScape was tested on gsnap, which is able to handle the circular mtDNA chromosome, but in theory any aligner should do.

1. Align your fastq files to the mitochondrial revised Cambridge reference sequence (rCRS) **only**. Note:
- The output bam filename **must** end in"\_MT.bam". 
- The output bam file **must** be sorted by coordinate (```samtools sort```).
- I suggest keeping only reads that are mapped properly (```samtools view -f 0x3```).
2. Using the file generated from step 1, create a bam file with MD tags using samtools. The file from this step **must** end in "\_MT\_MD.bam".
```
samtools calmd -e -Q --output-fmt BAM input_MT.bam > input_MT_MD.bam
```
3. Convert the files generated in step 1 ("\_MT.bam") to fastq. 
3. Align these fastq files in step 3 to the full human genome minus the rCRS. That is your reference fasta file must be the entire genome with chrM/MT removed. 
- The filename of created ***must*** end in "\_NT.bam".
- The output file **must** be sorted by **chromsome and location**.
- I suggest again that you keep only reads that are mapped properly (```samtools view -f 0x3```).

### MitoScape

4. Run mitoscape
```
java -Xmx16G -jar MitoScape-1.0.jar \
  --threads <number of threads> \
  --prob <classification probability> \
  --ld NUMTs_hg38.txt \
  --classifier MTClassifierModel.RF \
  --prefix <prefix, see below for explanation> \
  --out <output prefix for bam file>
```
The prefix matches the prefix of the bam filenames generated in the steps above. MitoScape will look for <prefix>\_MT\_MD.bam and <prefix>\_NT.bam. Be sure to "untar" the MTClassifierModel.RF.tar file.

### Variant Calling

5. In theory any variant caller will do, but MitoScape was tested on Mutect2 using mitochondrial mode.


Please feel free to contact me should you have questions. Please see the scripts directory for sample scripts.
