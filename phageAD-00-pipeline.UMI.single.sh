#!/usr/bin/env bash
set -e
# shopt -s expand_aliases
# source ~/.bashrc

# Setting up parameters, using environmental variables in bash
## Defining the path to scripts for this data analysis project
### CUSTOMIZE HERE (2 lines) ###
batchAccn='202304'
dataType='UMI_FatI'
### CUSTOMIZE ENDS ###

# Setting up / activating conda environment
conda activate "$CONDA_OVERALL_NAME"

mkdir "$analysisDataDir"

# Step 00:
# Creating data directory
# manually create analysis data directory so that the pipeline.log and .err can be stored for running this script
# mkdir "$analysisDataDir"

# Step 00:
# Linking and organizing raw sequencing data files
# create symbolic links for raw data files (.fastq.gz) in working directory (dataDir/00-seq)
# mkdir 00-seq/
# ln -s "$filename" 00-seq/"$renaming""$fileSuffix"

# Step 01:
# Performing quality control for raw data
echo '=== Step 01: quality control ==='
# echo 'Step 01: Performing quality control on raw reads'
## run FastQC, version 0.11.8
mkdir -p 01-qc/fastqc/
fastqc 00-seq/*.fastq.gz -o 01-qc/fastqc/ --noextract
## run MultiQC, versiion 0.9
mkdir 01-qc/multiqc/
multiqc 01-qc/fastqc/* -o 01-qc/multiqc/
# Step 02:
## trimming for spacer
echo '=== Step 02: trimming ==='
R26='attgtagcactgcgaaatgagaaagggagctacaac' # Repeat 26

mkdir 02-spacer/
echo '=====             Spacer               ====='
echo '=== removing Repeat 26 at 5 prime end  ==='
"$CUTADAPT_PATH" -g "$R26" --discard-untrimmed -o example.cutSpacer.fastq.gz 00-seq/example_R1.fastq.gz > example.cutSpacer.log 2> example.cutSpacer.err
echo '=== removing Repeat 26 at 3 prime end  ==='
"$CUTADAPT_PATH" -a "$R26" --discard-untrimmed --minimum-length 28 --maximum-length 32 --max-n 0 --quality-cutoff   30,30 -o example.spacer.fastq.gz example.cutSpacer.fastq.gz > example.spacer.log 2> example.spacer.err

## trimming for UMIs
pSp25_KNOWN11='gacgtagcgtccatgcgcggcgcattacctttac' # reverse strand
L1='gacgtagcgtcatgactttaacgcacgttcgcttatcgcaacggctg' # reverse strand
OYZ_ADAPTOR='actacgcacgcgacga' # NOTE: add sequence, as in reverse strand
UMI_LENGTH=8

if [ "$dataType" = 'UMI_FatI' ]
then
    constructSeq="$pSp25_KNOWN11"
elif [ "$dataType" = 'UMI_BsmAI' ]
then
    constructSeq="$L1"
else
    echo "Unexpected sample type"
    exit 1
fi

mkdir 02-umi/
echo '=====              UMI                 ====='
echo '=== removing partial SP25 or L1 at 3 prime end  ==='
"$CUTADAPT_PATH" -a "$constructSeq" --discard-untrimmed --minimum-length 8 -o example.cutUmi.fastq.gz example_R2.fastq.gz > example.cutUmi.log 2> example.cutUmi.err
echo '=== removing known adaptor sequence at 5 prime  ==='
"$CUTADAPT_PATH" -g "$OYZ_ADAPTOR"  --max-n 0 --quality-cutoff 30,30 --maximum-length 8 --minimum-length 8 -o example.umi.fastq.gz example.cutUmi.fastq.gz > example.umi.log 2> example.umi.err

echo 'Successfully extracted spacer and UMI sequences'

# Step 03:
# deduplicating spacer-umi combination
echo '=== Step 03: deduplication ==='
mkdir 03-dedup/

echo "=== REMOVE UNPAIRED SPACER-UMI SEQUENCES ==="
gunzip -c 02-spacer/example.spacer.fastq.gz > example.umi.seq.fastq
gunzip -c 02-umi/example.umi.fastq.gz > example.spacer.seq.fastq

"$FASTQPAIR_PATH" example.umi.seq.fastq example.spacer.seq.fastq >>example.fastq_pair.log 2>&1
rm example.umi.seq.fastq example.spacer.seq.fastq
gzip example.umi.seq.fastq.paired.fq
gzip example.spacer.seq.fastq.paired.fq

echo "=== DEDUPLICATION extract ==="
"$UMITOOLS_PATH" extract -I example.umi.seq.fastq.paired.fq.gz --bc-pattern=NNNNNNNN --read2-in=example.spacer.seq.fastq.paired.fq.gz --stdout=example.umi.extract.fastq.gz --read2-out=example.spacer.extract.fastq.gz --log=example.umitools.extract.log

# Step 04:
# aligning to N. meningitidis and phage genome
echo '=== Step 04: alignment ==='
echo '==== Aligning ===='

mkdir 04-align_ref/
"$BOWTIE2BUILD_PATH" 04-align_ref/Neisseria_meningitidis/NCBI/8013/NC_017501.1/GCF_000026965.1_ASM2696v1_genomic.fna,04-align_ref/download_ref_seq/MDA_kan.fa 04-align_ref/phageAD_2023-ref_Nme_MDA > 04-align_ref/phageAD_2023-ref_Nme_MDA.bt2build.log

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

echo '=== FINISH DEDUPLICATION ==='
"$UMITOOLS_PATH" dedup -I example.sorted.bam --output-stats=03-dedup/example.dedup.stats -S example.sorted.dedup.bam > example.dedup.log
"$SAMTOOLS_PATH" index example.sorted.dedup.bam

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
cat 04-align_ref/Neisseria_meningitidis/NCBI/8013/NC_017501.1/GCF_000026965.1_ASM2696v1_genomic.fna 04-align_ref/download_ref_seq/MDA_kan.fa > 04-align_ref/phageAD_2023-ref_Nme_MDA.fa

mkdir 05-pad/
echo '=== Step 05: padding ==='

echo '=== Padding forward alignments ==='
"$BEDTOOLS_COVERAGEBED" -bg -5 -ibam 04-align/example.uniq.forward.bam >example.forward.bedgraph
python2 "$padBedFilesByLength" example.forward.bedgraph 25 50 forward >example.forward.padded.bed
# reverse
echo "=== Padding reverse alignments ==="
"$BEDTOOLS_COVERAGEBED" -bg -5 -ibam 04-align/example.uniq.reverse.bam >example.reverse.bedgraph
python2 "$padBedFilesByLength" example.reverse.bedgraph 25 50 reverse >example.reverse.padded.bed
# merge
echo "=== Merging BED files of padded locations ==="
cat example.forward.padded.bed example.reverse.padded.bed >example.padded.bed
echo "=== Retrieving genomic sequences by BED into FASTA ==="
python3 "$getFastaFromBed" example.padded.bed 04-align_ref/phageAD_2023-ref_Nme_MDA.fa > example.padded.fasta

# Step 06
# correct to representative strand
mkdir 06-rc/
echo '=== Step 06: extraction of genomic contexts ==='
python3 "$getReverseComplement" example.padded.fasta 06-rc/example.padded.rc.fasta
grep -A 2 '>NC' 06-rc/example.padded.rc.fasta | sed '/^--$/d' > 06-rc/example.Nme.fa
grep -A 2 '>MDA' 06-rc/example.padded.rc.fasta | sed '/^--$/d' > 06-rc/example.MDA.fa

# Final:
# print success message
echo 'ANALYSIS WAS SUCCESSFUL'
