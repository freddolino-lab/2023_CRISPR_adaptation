# 2022 CRISPR apoCas9

For detailed descriptions of analysis, see the Materials and Methods section.
Data is available at NCBI GEO.

## Analysis steps

The numbering of analysis steps are:

00. retrieve raw sequencing reads (as FASTQ files)
01. quality control of raw sequencing reads
02. extract spacer sequences using adapter trimming of sequencing constructs, extract UMIs when applicable
03. when applicable, deduplicate spacer-UMI pairs
04. align spacer sequences to reference containing *Neisseria meningitidis* genome and MDA phage genome, and filter to keep only unique alignments using `bowtie2` flags
05. calculate the ranges of genomic positions after padding the aligned positions for spacers for upstream 25 bp and downstream 50 bp (BED files) and extract the sequences specific on the aligned strands (FASTA files) as genomic contexts
06. reverse complement the genomic context sequences for non-target sequence consistent for PAM identification conventions

The numbers in script names indicate the steps for which each scripts are written for, with scripts `phageAD-00-*` indicating the overall pipeline scripts or environment variables.

## Environments

+ `env/` contains two YML files of `conda` environments used in analysis
  + `adapt_py-3.7.yml`: Python 3.7, used for all Python scripts and other programs called from `bash` unless otherwise specified
  + `adapt_r-3.6.2.yml`: R 3.6.2 and corresponding packages, used for all R scripts
  + `adapt_umitools-3.7.yml`: environment where `umitools` is installed on top of the existing `adapt_py-3.7.yml`
+ Others:
  + `Cutadapt` 2.6
  + `ParDRe` 2.2.5
  + `bowtie2` 2.4.1
  + `samtools` 1.9 (using `htslib` 1.9)
  + `bedtools`: 2.19.1
  + `FastQC` 0.11.8
  + `MultiQC` 0.9
+ Genome: Neisseria meningitidis strain 8013 (RefSeq accession number: NC_017501.1, GenBank assembly accession: GCF_000026965.1, v1); MDA phage genomic sequence, supplied by collaborators

## Data directory structure

The sub-directories under the overall project data directories are organized by the analysis steps:

```
data directory
├ 00-seq/ # raw FASTQ files or symbolic links
├ 01-qc/ # data quality control
├ 02-spacer/ # extracted spacer sequences
├ (02-umi/) # extracted UMI sequences, UMI samples only
├ (03-dedupped/) # deduplication, spacer and UMI sequences, runtime logs, UMI samples only
├ 04-align/ # alignment to Nme genome and BAM file processing
├ 05-pad/ # extracted strand-specific genome context positions and sequences
└ 06-rc/ # non-target strand sequences for motif analysis
```
