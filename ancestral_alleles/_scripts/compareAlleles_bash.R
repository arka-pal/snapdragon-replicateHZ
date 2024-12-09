#### compareAlleles: Compute chromosome-wide allele comparisons between Amolle and Amajus
#### author: Arka Pal
#### written: 2024-Feb-1
#### last update: 2024-Jun-5

#### usage: Rscript compareAlleles_bash.R baseDIR chrom stitchRun outFile

## Read arguments from the bash
args <- commandArgs(trailingOnly = TRUE)
baseDIR <- as.character(args[1])
chrom <- as.character(args[2])
posFile <- as.character(args[3])
outFile <- as.character(args[4])

##
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(tictoc)

source("~/snap_hap_repHZ/ancestral_alleles/_scripts/functions_polarisation.R")

alleleSequence =  c('A','T','C','G','N','del')

alleleDat_chrom = compile_AlleleDat_chromSegments(alleleSequence, baseDIR, chrom, posFile)
fwrite(alleleDat_chrom, 
       file = outFile,
       sep = ',', 
       quote = F, 
       row.names = FALSE, 
       col.names = TRUE)