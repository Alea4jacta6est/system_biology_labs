# RNAseq assembly

Download fastq files from your dataset. Some datasets are bigger than others, in case you have more than 2 libraries, select two of your choice (preferably, coming from the same experimental condition).

Check read quality with Fastqc and multiqc

If you detect any problems with your read data, use tools to solve them: remove adapters with Trimmomatic, cutadapt or fastp. Remove possible overrepresented sequences/contaminants with sortmerna or with alignment tools like bbmap, bowtie, STAR (by aligning reads to contaminant reference)

Perform QC after read trimming/filtering. Repeat step 3 if you still got some problems with the dataset

Run Trinity and rnaspades to assemble de novo transcriptomes from your filtered read data. Important: If you don’t have enough computational resources use one library instead. If that is not enough, take half the reads in the library etc. until you can assemble RNAs.

Run rnaquast and BUSCO on resultant fasta files.

Report:

- Read QC before and after filtering/trimming

- rnaquast results for both rnaspades and Trinity assemblies

- BUSCO results for both assemblies (on the same figure)


Possible options (not mandatory, do only if you have free time and not much else to do):

Compare transcriptome assemblies on downsamples of libraries (How much worse is the transcriptome assembled from half the reads?)

Use Transrate for assembly QC

Some datasets were created from DE analysis. You could try calculating DE using assembled reference

Combine points 1 and 3. Is DE analysis suffering from poor assembly?
