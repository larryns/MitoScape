# MitoScape
A big-data, machine-learning workflow for obtaining mtDNA sequence from NGS data. For more details see this [video](https://www.youtube.com/watch?v=lB9kBeq8m8c).

## Cavatica Approach
Workflow and additional documenatation available at: : https://cavatica.sbgenomics.com/public/apps/d3b-bixu/app-publisher/mitoscape-wf/, including full instructions on how to run MitoScape using the Cavatica framework. If you require help on running MitoScape on Cavatica, please look at the Cavatica documentation and reach out to the Cavatica team. Note that you **must** create a cavatica account to access the URL above.

**Note:** The Cavatica approach uses an older Docker version of the code. I have since made several bug fixes and upgraded the Spark engine and libraries. The jar version described below is the latest and most well-tested version. As far as I know the Docker version in the Cavatica group has not been updated to use the current jar file, and I am no longer maintaining this version. Use at your own peril.

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
  --ld mitomap.ld \
  --numt NUMTs_hg38.txt \
  --classifier MTClassifierModel.RF \
  --prefix <prefix, see below for explanation> \
  --out <output prefix for bam file>
```
The prefix matches the prefix of the bam filenames generated in the steps above. MitoScape will look for <prefix>\_MT\_MD.bam and <prefix>\_NT.bam. Be sure to "untar" the MTClassifierModel.RF.tar file.

### Variant Calling

5. In theory any variant caller will do, but MitoScape was tested on Mutect2 using mitochondrial mode.


Please feel free to contact me should you have questions. Please see the scripts directory for sample scripts.
  
## Building MitoScape

The fat jar is somewhat onerous to maintain, and I will not be maintaining that going forward. The latest versions of sbt and java will not work to
create the fat jar in the current state. You would have to either use the older version of sbt and sbt-assembly or modify the build.sbt to get the sbt assembly to work. At some point, these older versions may no longer be avialablle. Note that you will also need to use Java 8 to get the fat far to build.
Given all these issues, I will be switching back to Docker in the near future, and only supporting Docker builds.

### Further Documentation
  
If you use MitoScape, please cite: 

```"MitoScape: A big-data, machine-learning platform for obtaining mitochondrial DNA from next-generation sequencing data". PLoS Comput Biol. 2021 Nov 11;17(11):e1009594. doi: 10.1371/journal.pcbi.1009594. eCollection 2021 Nov.```

More information on the algorithm and design are also included in this paper.
