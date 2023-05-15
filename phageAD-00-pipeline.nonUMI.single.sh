#!/usr/bin/env bash
set -e

# Setting up parameters, using environmental variables in bash
## Defining the path to scripts for this data analysis project
### CUSTOMIZE HERE (2 lines) ###
batchAccn='202304'
dataType='nonUMI'
### CUSTOMIZE ENDS ###

projectDir='/home/diaorch/projects/crisprAdaptation/crisprAdaptation-'"$batchAccn"'/'
source "$projectDir"'crispr'"$batchAccn"'_adapt-00-envSetUp.sh'

# Setting up / activating conda environment
conda activate "$CONDA_OVERALL_NAME"

# Step 00 to Step 01 are the same as UMI sample analysis

# Step 02:
## trimming for spacer
echo '=== Step 02: trimming ==='
# LEADS_R26='cttatgaaataaggatttcccgtcgaagtattgtagcactgcgaaatgagaaagggagctacaac' # LeadS + 5 b gap + Repeat 26
R26='attgtagcactgcgaaatgagaaagggagctacaac' # Repeat 26

mkdir 02-spacer/
echo '=====             Spacer               ====='
echo '=== removing Repeat 26 at 5 prime end  ==='
"$CUTADAPT_PATH" -g "$R26" --discard-untrimmed -o example.cutSpacer.fastq.gz 00-seq/example_R1.fastq.gz > example.cutSpacer.log 2> example.cutSpacer.err
echo '=== removing Repeat 26 at 3 prime end  ==='
"$CUTADAPT_PATH" -a "$R26" --discard-untrimmed --minimum-length 28 --maximum-length 32 --max-n 0 --quality-cutoff   30,30 -o example.spacer.fastq.gz example.cutSpacer.fastq.gz > example.spacer.log 2> example.spacer.err

# extracting UMI - not applicable for non-UMI samples

# Step 03:
# deduplicating spacer-umi combination
# not applicable for non-UMI samples

# Step 04:
# building alignment reference and aligning to reference genome
# building alignment reference is the same as in the UMI sample analysis

mkdir 04-align/
printf '' > align.all.log
echo '=== Aligning to custom genome ==='
"$BOWTIE2_PATH" --very-sensitive -x 04-align_ref/phageAD_2023-ref_Nme_MDA -U 03-dedup/example.spacer.extract.fastq.gz -S example.sam >> example.align.log 2>&1
cat example.align.log >> align.all.log
echo "=== Converting SAM to sorted BAM ==="
"$SAMTOOLS_PATH" view -Sbh example.sam > example.bam
"$SAMTOOLS_PATH" sort -o example.sorted.bam example.bam
rm example.sam
rm example.bam
"$SAMTOOLS_PATH" index example.sorted.bam
echo "=== Filter DEDUP bam file for uniquely mapped reads ==="
echo "===== forward ====="
"$SAMTOOLS_PATH" view -h -F 20 example.sorted.dedup.bam | grep -v 'XS:i' > example.uniq.forward.sam
if grep -v -q '^@' example.uniq.forward.sam;
then
    "$SAMTOOLS_PATH" view -Sbh example.uniq.forward.sam -o example.uniq.forward.bam
else
    echo 'ERROR: no uniquely aligned forward reads for sample: example'
    exit 1
fi
rm example.uniq.forward.sam
"$SAMTOOLS_PATH" index example.uniq.forward.bam
echo "===== reverse ====="
"$SAMTOOLS_PATH" view -h -f 16 -F 4 example.sorted.dedup.bam | grep -v 'XS:i' > example.uniq.reverse.sam
if grep -v -q '^@' example.uniq.reverse.sam;
then
    "$SAMTOOLS_PATH" view -Sbh example.uniq.reverse.sam -o example.uniq.reverse.bam
else
    echo 'ERROR: no uniquely aligned reverse reads for sample: example'
    exit 1
fi
rm example.uniq.reverse.sam
"$SAMTOOLS_PATH" index example.uniq.reverse.bam

# Step 05:
# getting padded sequences
mkdir 05-pad/
echo '=== Step 05: padding ==='

echo '=== Padding forward alignments ==='
"$BEDTOOLS_COVERAGEBED" -bg -5 -ibam 04-align/example.uniq.forward.bam > example.forward.bedgraph
python2 "$padBedFilesByLength" example.forward.bedgraph 25 50 forward > example.forward.padded.bed
# reverse
echo "=== Padding reverse alignments ==="
"$BEDTOOLS_COVERAGEBED" -bg -5 -ibam 04-align/example.uniq.reverse.bam >example.reverse.bedgraph
python2 "$padBedFilesByLength" example.reverse.bedgraph 25 50 reverse >example.reverse.padded.bed
# merge
echo "=== Merging BED files of padded locations ==="
cat example.forward.padded.bed example.reverse.padded.bed >example.padded.bed
# NOTE: beginning of edits
# check spacer% at each site
echo "=== Apply max spacer% at each uniq spacer site ==="

conda activate bio_r-3.6.2
Rscript "$applyJackpotCap" example.padded.bed 0.05
conda deactivate

echo "=== Retrieving genomic sequences by BED into FASTA ==="
conda activate "$CONDA_OVERALL_NAME"
python3 "$getFastaFromBed" example.padded.bed 04-align_ref/phageAD_2023-ref_Nme_MDA.fa > example.padded.fasta

# Step 06
# visualizaing PAMs
# Part 1 and 2
# reverse complementing extracted sequences, calculate PWM
mkdir 06-rc/
echo '=== Step 06: extraction of genomic contexts ==='
python3 "$getReverseComplement" 05-pad/example.padded.fa 06-rc/example.rc.fasta

# Step 07
# separating files by Nme or MDA-phage sources
source phageAD-06-analysisByChr.sh

echo 'ANALYSIS WAS SUCCESSFUL'
