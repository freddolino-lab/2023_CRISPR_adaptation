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
# mkdir "$rawSubDir"
# ln -s "$filename" "$rawSubDir""$renaming""$fileSuffix"

# Step 01:
# Performing quality control for raw data
echo '=== Step 01: quality control ==='
# echo 'Step 01: Performing quality control on raw reads'
mkdir "$qcSubDir"
## run FastQC, version 0.11.8
mkdir "$qcSubDir"'fastqc/'
fastqc "$rawSubDir""$dataPrefix"*"$rawSuffix" -o "$qcSubDir"'fastqc/' --noextract
## run MultiQC, versiion 0.9
mkdir "$qcSubDir"'multiqc/'
multiqc "$qcSubDir"'fastqc/'* -o "$qcSubDir"'multiqc/'
# Step 02:
## trimming for spacer
echo '=== Step 02: trimming ==='
R26='attgtagcactgcgaaatgagaaagggagctacaac' # Repeat 26

mkdir "$newSpacerSubDir"
r1Suffix='_R1'
for filename in "$rawSubDir""$dataPrefix"*"$r1Suffix""$rawSuffix";
do
    echo $filename
    filenameBase=${filename#$rawSubDir$dataPrefix}
    i=${filenameBase%$r1Suffix$rawSuffix}

    naming="$newSpacerSubDir""$dataPrefix""$i"

    echo "$i"
    echo 'Sample '"$i"
    echo '=====             Spacer               ====='
    echo '=== removing Repeat 26 at 5 prime end  ==='
    "$CUTADAPT_PATH" -g "$R26" --discard-untrimmed -o "$naming""$cutSpacerSeqFqSuffix" "$filename" > "$naming"'.cutSpacer.log' 2> "$naming"'.cutSpacer.err'
    # echo '===     getting new Spacer sequences     ==='
    echo '=== removing Repeat 26 at 3 prime end  ==='
    "$CUTADAPT_PATH" -a "$R26" --discard-untrimmed --minimum-length 28 --maximum-length 32 --max-n 0 --quality-cutoff 30,30 -o "$naming""$spacerSeqFqSuffix" "$naming""$cutSpacerSeqFqSuffix" > "$naming"'.spacer.log' 2> "$naming"'.spacer.err'
done

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

mkdir "$umiSubDir"
r2Suffix='_R2'
for filename in "$rawSubDir""$dataPrefix"*"$r2Suffix""$rawSuffix";
do
    echo $filename
    filenameBase=${filename#$rawSubDir$dataPrefix}
    i=${filenameBase%$r2Suffix$rawSuffix}
    echo 'Sample '"$i"

    naming="$umiSubDir""$dataPrefix""$i"

    echo '=====              UMI                 ====='
    echo '=== removing partial SP25 or L1 at 3 prime end  ==='
    "$CUTADAPT_PATH" -a "$constructSeq" --discard-untrimmed --minimum-length 8 -o "$naming""$cutUmiSeqFqSuffix" "$filename" > "$naming"'.cutUmi.log' 2> "$naming"'.cutUmi.err'
    echo '=== removing known adaptor sequence at 5 prime  ==='
    "$CUTADAPT_PATH" -g "$OYZ_ADAPTOR"  --max-n 0 --quality-cutoff 30,30 --maximum-length 8 --minimum-length 8 -o "$naming""$umiSeqFqSuffix" "$naming""$cutUmiSeqFqSuffix" > "$naming"'.umi.log' 2> "$naming"'.umi.err'
done

echo 'Successfully extracted spacer and UMI sequences'

# Step 03:
# deduplicating spacer-umi combination
echo '=== Step 03: deduplication ==='
mkdir "$dedupSubDir"
for filename in "$umiSubDir"*"$umiSeqFqSuffix";
do
    echo "$filename"
    i=$(basename $filename | awk -F '.' '{print $1}')
    echo 'Sample '"$i"

    spacerSeqFilename="$newSpacerSubDir""$i""$spacerSeqFqSuffix"
    umiSeqFilename="$umiSubDir""$i""$umiSeqFqSuffix"

    naming="$dedupSubDir""$i"
    if test -f "$naming"'.spacer.extract.fastq.gz'; then
        continue
    fi


    echo "=== REMOVE UNPAIRED SPACER-UMI SEQUENCES ==="
    gunzip -c "$umiSeqFilename" > "$naming".umi.seq.fastq
    gunzip -c "$spacerSeqFilename" > "$naming".spacer.seq.fastq
    echo "$naming"
    "$FASTQPAIR_PATH" "$naming".umi.seq.fastq "$naming".spacer.seq.fastq >>"$naming".fastq_pair.log 2>&1
    rm "$naming".umi.seq.fastq "$naming".spacer.seq.fastq
    gzip "$naming".umi.seq.fastq.paired.fq
    gzip "$naming".spacer.seq.fastq.paired.fq

    echo "=== DEDUPLICATION extract ==="
    "$UMITOOLS_PATH" extract -I "$naming".umi.seq.fastq.paired.fq.gz --bc-pattern=NNNNNNNN --read2-in="$naming".spacer.seq.fastq.paired.fq.gz --stdout="$naming".umi.extract.fastq.gz --read2-out="$naming".spacer.extract.fastq.gz --log="$naming".umitools.extract.log
done

# Step 04:
# aligning to N. meningitidis and phage genome
echo '=== Step 04: alignment ==='
echo '==== Aligning ===='
readPath="$dedupSubDir"
readSuffix='.spacer.extract.fastq.gz'
logFullFilename="$alignLog"

##### DATA 202211: phage genome needed for alignment
mkdir "$alignRefSubDir"
"$BOWTIE2BUILD_PATH" "$NmeGenomeFa","$phageFa" "$bowtie2index" > "$bowtie2index"'.bt2build.log'
sampleBt2Idx="$bowtie2index"

mkdir "$alignSubDir"
printf '' > "$logFullFilename"

for filename in "$readPath"*"$readSuffix"
do
    echo "$filename"
    filenameBase=${filename#$readPath}
    i=$(echo "$filenameBase" | awk -F '-' '{print $2}')
    # echo "$i"
    echo 'Sample '"$i"

    sampleName=${filenameBase%$readSuffix}
    echo 'Sample name: '"$sampleName" >> "$logFullFilename"
    refNameBase=$(basename "$sampleBt2Idx")
    echo 'Reference genome: '"$refNameBase" >> "$logFullFilename"

    naming="$alignSubDir""$sampleName"

    if test -f "$naming"'.uniq.reverse.bam'; then
        continue
    fi

    echo '=== Aligning to custom genome ==='
    "$BOWTIE2_PATH" --very-sensitive -x "$sampleBt2Idx" -U "$filename" -S "$naming".sam >> "$naming"'.align.log' 2>&1
    cat "$naming"'.align.log' >> "$logFullFilename"

    echo "=== Converting SAM to sorted BAM ==="
    "$SAMTOOLS_PATH" view -Sbh "$naming"'.sam' > "$naming"'.bam'
    "$SAMTOOLS_PATH" sort -o "$naming"'.sorted.bam' "$naming"'.bam'
    rm "$naming"'.sam'
    rm "$naming"'.bam'
    "$SAMTOOLS_PATH" index "$naming"'.sorted.bam'

    echo '=== FINISH DEDUPLICATION ==='
    "$UMITOOLS_PATH" dedup -I "$naming"'.sorted.bam' --output-stats="$dedupSubDir$sampleName"'.dedup.stats' -S "$naming"'.sorted.dedup.bam' > "$naming"'.dedup.log'
    "$SAMTOOLS_PATH" index "$naming"'.sorted.dedup.bam'


    echo "=== Filter DEDUP bam file for uniquely mapped reads ==="
    echo "===== forward ====="
    "$SAMTOOLS_PATH" view -h -F 20 "$naming"'.sorted.dedup.bam' | grep -v 'XS:i' > "$naming"'.uniq.forward.sam'
    if grep -v -q '^@' "$naming"'.uniq.forward.sam';
    then
        "$SAMTOOLS_PATH" view -Sbh "$naming"'.uniq.forward.sam' -o "$naming"'.uniq.forward.bam'
    else
        echo 'ERROR: no uniquely aligned forward reads for sample '"$sampleName"
        exit 1
    fi
    rm "$naming"'.uniq.forward.sam'
    "$SAMTOOLS_PATH" index "$naming"'.uniq.forward.bam'
    echo "===== reverse ====="
    "$SAMTOOLS_PATH" view -h -f 16 -F 4 "$naming"'.sorted.dedup.bam' | grep -v 'XS:i' > "$naming"'.uniq.reverse.sam'
    if grep -v -q '^@' "$naming"'.uniq.reverse.sam';
    then
        "$SAMTOOLS_PATH" view -Sbh "$naming"'.uniq.reverse.sam' -o "$naming"'.uniq.reverse.bam'
    else
        echo 'ERROR: no uniquely aligned reverse reads for sample '"$sampleName"
        exit 1
    fi
    rm "$naming"'.uniq.reverse.sam'
    "$SAMTOOLS_PATH" index "$naming"'.uniq.reverse.bam'
done


# Step 05:
# getting padded sequences
cat "$NmeGenomeFa" "$phageFa" > "$ALL_CHR_FASTA"

mkdir "$padSubDir"
echo '=== Step 05: padding ==='

for filename in "$alignSubDir"*'.uniq.forward.bam'
do
    echo "$filename"
    # filenameBase=${filename##"$alignSubDir"}
    filenameBase=$(basename $filename | awk -F '.' '{print $1}')
    echo 'Sample name (base name): '"$filenameBase"

    naming="$padSubDir""$filenameBase"

    echo '=== Padding forward alignments ==='
    "$BEDTOOLS_COVERAGEBED" -bg -5 -ibam "$alignSubDir""$filenameBase"'.uniq.forward.bam' >"$naming"'.forward.bedgraph'
    python2 "$padBedFilesByLength" "$naming"'.forward.bedgraph' 25 50 forward >"$naming"'.forward.padded.bed'
    # reverse
    echo "=== Padding reverse alignments ==="
    "$BEDTOOLS_COVERAGEBED" -bg -5 -ibam "$alignSubDir""$filenameBase"'.uniq.reverse.bam' >"$naming"'.reverse.bedgraph'
    python2 "$padBedFilesByLength" "$naming"'.reverse.bedgraph' 25 50 reverse >"$naming"'.reverse.padded.bed'
    # merge
    echo "=== Merging BED files of padded locations ==="
    cat "$naming"'.forward.padded.bed' "$naming"'.reverse.padded.bed' >"$naming"'.padded.bed'
    echo "=== Retrieving genomic sequences by BED into FASTA ==="
    python3 "$getFastaFromBed" "$naming"'.padded.bed' "$ALL_CHR_FASTA" >"$naming""$paddedFaSuffix"
done

# Step 06
mkdir "$rcSubDir"
mkdir "$pwmSubDir"
echo '=== Step 06: extraction of genomic contexts ==='
for filename in "$padSubDir"*"$paddedFaSuffix"
do
    echo "$filename"
    tmp=${filename%%"$paddedFaSuffix"}
    filenameBase=${tmp#"$padSubDir"}
    echo 'Sample name (base name): '"$filenameBase"
    python3 "$getReverseComplement" "$filename" "$rcSubDir""$filenameBase""$rcPaddedFaSuffix"
done

# Step 07
# separating files by Nme or MDA-phage sources
source "$projectDir"'crispr'"$batchAccn"'_adapt-06-analysisByChr.sh'

# Final:
# print success message
echo 'ANALYSIS WAS SUCCESSFUL'

RCLONE_PATH="$HOME"'/packages/rclone-v1.50.2-linux-amd64/rclone'
"$RCLONE_PATH" copy "$rcSubDir" 'gdrive:cli_upload_rclone_tmp/data/crisprAdaptation-'"$batchAccn"'/'"$dataType"'/061-rc'
echo 'FASTA FILE UPLOADED'
