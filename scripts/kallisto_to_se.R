#!/usr/bin/env Rscript

library(tximport)
library(SummarizedExperiment)
library(Matrix)
library(tidyverse)
library(magrittr)
library(S4Vectors)

################################################################################
                                        # Load Argument Values
################################################################################
if (interactive()) {
    args = c("data/kallisto/utrome_r2",
             "/data/mayrc/data/mca/gff/utrome.e30.t5.gc25.pas3.f0.9999.w500.m200.tsv",
             "/scratch/fanslerm/hspcs_bulk.txs.rds")
} else{
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 3) {
        stop("Incorrect number of arguments!\nUsage:\n> kallisto_to_se.R <inputPath> <txToGeneFile> <txsOutFile>\n")
    }
}
arg.inputPath <- args[1]
arg.txToGeneFile <- args[2]
arg.txsOutFile <- args[3]

################################################################################
                                        # Load Data
################################################################################
## Load transcript annotations
tx2gene <- read_tsv(arg.txToGeneFile, skip=1, col_names=c("transcript_id", "transcript_id_merged", "gene_id"))

## get filenames
sample.files <- list.files(path=arg.inputPath, pattern="abundance.h5", full.names=TRUE, recursive=TRUE) %>%
    set_names(str_extract(., '(HSC|MPP)[1-4]?[a-d]'))

txsOut <- tximport(files=sample.files, type="kallisto", countsFromAbundance="no", txOut=TRUE)

se <- SummarizedExperiment(assays=list(counts=Matrix(txsOut$counts, sparse=TRUE)))

colData(se) <- colnames(se) %>%
    { DataFrame(sample=.,
                cell_type=str_extract(., '(HSC|MPP)[1-4]?'),
                batch=str_extract(., '[a-d]$'),
                row.names=.) }
rowData(se)['transcript_id'] <- rownames(se)

rowData(se) %<>% merge(tx2gene, by="transcript_id", sort=FALSE)

saveRDS(se, arg.txsOutFile)
