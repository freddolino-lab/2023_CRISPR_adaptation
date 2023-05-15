#!/usr/bin/env bash
set -e

projectDir='/home/diaorch/projects/crisprAdaptation/crisprAdaptation-'"$batchAccn"'/'
source "$projectDir"'crispr'"$batchAccn"'_adapt-00-envSetUp.sh'

# Setting up / activating conda environment
conda activate "$CONDA_OVERALL_NAME"

# Step 06
echo '=== Step 06: visualization ==='
for filename in "$padSubDir"*"$paddedFaSuffix"
do
    echo "$filename"
    tmp=${filename%%"$paddedFaSuffix"}
    filenameBase=${tmp#"$padSubDir"}
    echo 'Sample name (base name): '"$filenameBase"
    # Part 1: reverse complement
    # python3 "$getReverseComplement" "$filename" "$rcSubDir""$filenameBase""$rcPaddedFaSuffix"
    # Part 1.5: separate FASTA sequences by "chromosome"
    grep -A 2 '>NC' "$rcSubDir""$filenameBase""$rcPaddedFaSuffix" | sed '/^--$/d' > "$rcSubDir""$filenameBase"'.Nme.fa'
    grep -A 2 '>MDA' "$rcSubDir""$filenameBase""$rcPaddedFaSuffix" | sed '/^--$/d' > "$rcSubDir""$filenameBase"'.MDA.fa'
done
