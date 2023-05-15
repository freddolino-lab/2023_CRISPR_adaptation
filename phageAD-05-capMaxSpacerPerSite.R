#!/usr/bin/env Rscript

# Project: CRISPR Adaptation apoCas9
# Goal: apply maximum spacer-per-site (or unique spacer) count caps to handle jackpot sites
# Author: diaorch
# Start date: 20230501

capThreshold <- 0.05

args <- commandArgs(trailingOnly = T)
# fn <- '~/data/crisprAdaptation-202303/nonUMI/05-pad/crispr_202303-22-1012.padded.bed'
fn <- args[1]
capThreshold <- as.numeric(args[2])


bedIn <- read.table(fn, sep = '\t', header = F, stringsAsFactors = F)
# print(head(bedIn))
# print(str(bedIn))

bedOut <- bedIn[NULL, ]

for (chrom in c('NC_017501.1', 'MDA_kan')){
  chromBed <- bedIn[bedIn[, 1] == chrom, ]
  chromScoreSum <- sum(chromBed[, 5])
  scoreCap <- round(chromScoreSum * capThreshold)
  chromBed[chromBed[, 5] > scoreCap, 5] <- scoreCap
  # print(chromBed[chromBed[, 5] == scoreCap, ])
  bedOut <- rbind(bedOut, chromBed)
}

write.table(x = bedOut, file = gsub(x = fn, pattern = 'padded', replacement = 'pad.cap'), sep = '\t', row.names = F, col.names = F, quote = F)
