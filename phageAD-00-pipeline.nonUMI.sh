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

# Step 00 to Step 01 are the same as non-UMI sample analysis


# Step 02:
## trimming for spacer
echo '=== Step 02: trimming ==='
# LEADS_R26='cttatgaaataaggatttcccgtcgaagtattgtagcactgcgaaatgagaaagggagctacaac' # LeadS + 5 b gap + Repeat 26
R26='attgtagcactgcgaaatgagaaagggagctacaac' # Repeat 26

mkdir "$newSpacerSubDir"
r1Suffix='_R1' # NOTE update raw sequence file name
for filename in "$rawSubDir""$dataPrefix"*"$r1Suffix""$rawSuffix";
do
    echo $filename
    filenameBase=${filename#$rawSubDir$dataPrefix}
    i=${filenameBase%$r1Suffix$rawSuffix}
    echo 'Sample '"$i"
    naming="$newSpacerSubDir""$dataPrefix""$i"
    echo '=====             Spacer               ====='
    echo '=== removing Repeat 26 at 5 prime end  ==='
    "$CUTADAPT_PATH" -g "$R26" --discard-untrimmed -o "$naming""$cutSpacerSeqFqSuffix" "$filename" > "$naming"'.cutSpacer.log' 2> "$naming"'.cutSpacer.err'
    # echo '===     getting new Spacer sequences     ==='
    echo '=== removing Repeat 26 at 3 prime end  ==='
    "$CUTADAPT_PATH" -a "$R26" --discard-untrimmed --minimum-length 28 --maximum-length 32 --max-n 0 --quality-cutoff 30,30 -o "$naming""$spacerSeqFqSuffix" "$naming""$cutSpacerSeqFqSuffix" > "$naming"'.spacer.log' 2> "$naming"'.spacer.err'
done


# Step 03:
# deduplicating spacer-umi combination
# not applicable for non-UMI samples

# Step 04:
echo '==== Aligning ===='
readPath="$newSpacerSubDir"
readSuffix="$spacerSeqFqSuffix"
logFullFilename="$alignLog"

mkdir "$alignSubDir"
printf '' >> "$logFullFilename"

for filename in "$readPath"*"$readSuffix"
do
    # echo "$filename"
    tmp=${filename%"$readSuffix"}
    sampleName=${tmp#"$readPath"}

    echo 'Sample name (base name): '"$sampleName"
    echo 'Sample name (base name): '"$sampleName" >> "$logFullFilename"
    echo '=== Aligning to custom genome ==='

    naming="$alignSubDir""$sampleName"

    if test -f "$naming"'.uniq.reverse.bam'; then
        continue
    fi

    "$BOWTIE2_PATH" --very-sensitive -x "$bowtie2index" -U "$filename" -S "$naming".sam >> "$naming"'.align.log' 2>&1
    cat "$naming"'.align.log' >> "$logFullFilename"

    echo "=== Converting SAM to sorted BAM ==="
    "$SAMTOOLS_PATH" view -Sbh "$naming"'.sam' -o "$naming"'.bam'
    "$SAMTOOLS_PATH" sort -o "$naming"'.sorted.bam' "$naming"'.bam'
    rm "$naming"'.sam'
    rm "$naming"'.bam'
    "$SAMTOOLS_PATH" index "$naming"'.sorted.bam'
    echo "=== Filter sorted bam file for uniquely mapped reads ==="
    echo "===== forward ====="
    "$SAMTOOLS_PATH" view -h -F 20 "$naming".'sorted.bam' | grep -v 'XS:i' > "$naming"'.uniq.forward.sam'
    if grep -v -q '^@' "$naming"'.uniq.forward.sam';
    then
        "$SAMTOOLS_PATH" view -Sbh "$naming"'.uniq.forward.sam' -o "$naming"'.uniq.forward.bam'
    else
        echo 'ERROR: no uniquely aligned forward reads for sample '"$filenameBase"
        exit 1
    fi
    rm "$naming"'.uniq.forward.sam'
    "$SAMTOOLS_PATH" index "$naming"'.uniq.forward.bam'
    echo "===== reverse ====="
    "$SAMTOOLS_PATH" view -h -f 16 -F 4 "$naming"'.sorted.bam' | grep -v 'XS:i' > "$naming"'.uniq.reverse.sam'
    if grep -v -q '^@' "$naming"'.uniq.reverse.sam';
    then
        "$SAMTOOLS_PATH" view -Sbh "$naming"'.uniq.reverse.sam' -o "$naming"'.uniq.reverse.bam'
    else
        echo 'ERROR: no uniquely aligned reverse reads for sample '"$filenameBase"
        exit 1
    fi
    rm "$naming"'.uniq.reverse.sam'
    "$SAMTOOLS_PATH" index "$naming"'.uniq.reverse.bam'
done

# Step 05:
# getting padded sequences
mkdir "$padSubDir"
echo '=== Step 05: padding ==='

for filename in "$alignSubDir"*'.uniq.forward.bam'
do
    echo "$filename"
    filenameBase=$(basename $filename | awk -F '.' '{print $1}')
    echo 'Sample name (base name): '"$filenameBase"
    naming="$padSubDir""$filenameBase"
    echo '=== Padding forward alignments ==='
    "$BEDTOOLS_COVERAGEBED" -bg -5 -ibam "$alignSubDir""$filenameBase"'.uniq.forward.bam' >"$naming"'.forward.bedgraph'
    python2 "$padBedFilesByLength" "$naming"'.forward.bedgraph' 25 50 forward >"$naming"'.forward.padded.bed'
    # reverse
    echo "=== Padding reverse alignments ==="
    "$BEDTOOLS_COVERAGEBED" -bg -5 -ibam "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.reverse.bam' >"$naming"'.reverse.bedgraph'
    python2 "$padBedFilesByLength" "$naming"'.reverse.bedgraph' 25 50 reverse >"$naming"'.reverse.padded.bed'
    # merge
    echo "=== Merging BED files of padded locations ==="
    cat "$naming"'.forward.padded.bed' "$naming"'.reverse.padded.bed' >"$naming"'.padded.bed'
    # NOTE: beginning of edits
    # check spacer% at each site
    echo "=== Apply max spacer% at each uniq spacer site ==="

    conda activate bio_r-3.6.2
    Rscript "$applyJackpotCap" "$naming"'.padded.bed' 0.05
    conda deactivate

    echo "=== Retrieving genomic sequences by BED into FASTA ==="
    conda activate "$CONDA_OVERALL_NAME"
    python3 "$getFastaFromBed" "$naming"'.pad.cap.bed' "$ALL_CHR_FASTA" >"$naming""$paddedFaSuffix"
    # NOTE: ending of edits
done

# Step 06
# visualizaing PAMs
# Part 1 and 2
# reverse complementing extracted sequences, calculate PWM
mkdir "$rcSubDir"
mkdir "$pwmSubDir"
echo '=== Step 06: visualization ==='
for filename in "$padSubDir"*"$paddedFaSuffix"
do
    echo "$filename"
    tmp=${filename%%"$paddedFaSuffix"}
    filenameBase=${tmp#"$padSubDir"}
    echo 'Sample name (base name): '"$filenameBase"
    # Part 1: reverse complement
    python3 "$getReverseComplement" "$filename" "$rcSubDir""$filenameBase""$rcPaddedFaSuffix"
done

# Step 07
# separating files by Nme or MDA-phage sources
source "$projectDir"'crispr'"$batchAccn"'_adapt-06-analysisByChr.sh'

echo 'ANALYSIS WAS SUCCESSFUL'

RCLONE_PATH="$HOME"'/packages/rclone-v1.50.2-linux-amd64/rclone'
"$RCLONE_PATH" copy "$rcSubDir" 'gdrive:cli_upload_rclone_tmp/data/crisprAdaptation-'"$batchAccn"'/'"$dataType"'/061-rc'
echo 'FASTA FILE UPLOADED'
