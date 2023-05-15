#!/usr/bin/env bash
# 202304
# CRISPR adaptation project with Yan Zhang Lab

## ===== program paths =====
# bcl2fastq='/home/petefred/bin/bcl2fastq'
CONDA_OVERALL_NAME='adapt_py-3.7'
CUTADAPT_PATH="$CONDA_PATH"'cutadapt'
BOWTIE2_PATH='/usr/bin/bowtie2' # ver 2.4.1
BOWTIE2BUILD_PATH='/usr/bin/bowtie2-build'
SAMTOOLS_PATH='/usr/bin/samtools'
BEDTOOLS_PATH="$HOME"'/src/bedtools2/bin/bedtools'
BEDTOOLS_COVERAGEBED="$HOME"'/src/bedtools2/bin/genomeCoverageBed'
FASTQPAIR_PATH="$HOME"'/packages/github/fastq-pair/build/fastq_pair'
UMITOOLS_PATH="$HOME"'/local/miniconda3/envs/adapt_umitools-3.7/bin/umi_tools'

## ===== homemade script paths =====
padBedFilesByLength='apoCas9-05-getPaddedBed.py'
getFastaFromBed='apoCas9-05-getFastaFromBed.py'
getReverseComplement='apoCas9-06-revCompFasta.py'
