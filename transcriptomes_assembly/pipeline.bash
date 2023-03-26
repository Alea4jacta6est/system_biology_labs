#!/bin/bash

##### Step 1: preparation
# starting with only one sample
# downloading data from SRA https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA453118&o=acc_s%3Aa
# SRR7058642 - control rep 3 for CDDP; SRR7058643 - cisplatinum rep 3
CURRENT_DIR=`pwd`
DATADIR="$CURRENT_DIR/data"; mkdir -p "$DATADIR"

for file in $DATADIR/*.fastq.gz; do
TAG=$(basename $file .fastq.gz);
echo $TAG
# TAG=SRR7058642

date
echo "STARTED DOWNLOADING THE DATA"
fastq-dump -O "$DATADIR" --gzip $TAG --split-3
date

##### Step 2: fastqc

ANALYSISDIR="$CURRENT_DIR/analysis"; mkdir -p "$ANALYSISDIR"
OUTDIR="${ANALYSISDIR}/fastqc/${TAG}"; mkdir -p "$OUTDIR"
fastqc -o "$OUTDIR" "${DATADIR}/${TAG}.fastq.gz" |& tee "$OUTDIR/${TAG}.fastqc.log"


#### Step 3: multiqc for everything

multiqc -x .Rproj.user -f "$CURRENT_DIR/analysis" -o "$CURRENT_DIR/analysis" #; done
 
#### Step 4: remove possible overrepresented sequences/contaminants

fastp -i "${DATADIR}/${TAG}.fastq.gz" -o "${DATADIR}/${TAG}_postprocessed.fastq.gz" --overrepresentation_analysis -P 20; done

Step 5: QC after read trimming/filtering

POSTOUTDIR="${ANALYSISDIR}/fastqc_postprocessed/${TAG}"; mkdir -p "$POSTOUTDIR"
fastqc -o "$POSTOUTDIR" "${DATADIR}/${TAG}_postprocessed.fastq.gz" |& tee "$OUTDIR/${TAG}_postprocessed.fastqc.log"

TRIMDIR="$ANALYSISDIR/trimming"
mkdir -p "$TRIMDIR"
TRIMLOG="$TRIMDIR/${TAG}.trimming.log"
TRIMMEDFASTQ="$TRIMDIR/${TAG}.trimmed.fastq.gz"
trimmomatic SE -threads 4 -phred33 "${DATADIR}/${TAG}_postprocessed.fastq.gz" "$TRIMMEDFASTQ" ILLUMINACLIP:adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 |& tee "$TRIMLOG"
fastqc -o "$TRIMDIR" "$TRIMMEDFASTQ" |& tee "$TRIMDIR/${TAG}.fastqc.log"; done

multiqc -x .Rproj.user -f "$CURRENT_DIR/analysis" -o "$CURRENT_DIR/analysis" # ; done

# Step 6: Run Trinity and SPAdes
TRINITYDIR="$ANALYSISDIR/trinity"
mkdir -p "$TRINITYDIR"
TRINITYLOG="$TRINITYDIR/${TAG}.trinity.log"
FASTQ_1="${DATADIR}/${TAG}_1.fastq"
FASTQ_2="${DATADIR}/${TAG}_2.fastq"
Trinity --seqType fq --left $FASTQ_1 --right $FASTQ_2 --max_memory 12G --min_kmer_cov 3 --output $TRINITYDIR |& tee $TRINITYLOG

SPADESDIR="$ANALYSISDIR/spades"
mkdir -p "$SPADESDIR"
SPADESLOG="$SPADESDIR/${TAG}.rnaspades.log"
python rnaspades.py -1 $FASTQ_1 -2 $FASTQ_2 -o "$SPADESDIR" -t 4 |& tee $SPADESLOG; done
# python rnaspades.py -1 data/SRR7058643_1.fastq.gz -2 data/SRR7058643_2.fastq.gz -o analysis/spades -t 4

python rnaQUAST.py -r reference.fasta -c $ANALYSISDIR/spades/transcripts.fasta -o $ANALYSISDIR/rnaquast_results

# Step 7: Run BUSCO
busco --in $ANALYSISDIR/spades/SRR7058642/transcripts.fasta -o $ANALYSISDIR/busco_output_SRR7058642 -l arthropoda_odb10 -m tran -f
generate_plot.py -wd analysis/busco/busco_SRR7058643