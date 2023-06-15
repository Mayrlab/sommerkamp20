#!/usr/bin/env Rscript

library(SummarizedExperiment)
library(Matrix)
library(readr)
library(dplyr)
library(stringr)

if (interactive()) {
    args <- c("data/se/hspcs_bulk.txs.raw.rds",
              "/data/mayrc/data/atlas-mm/utrs/utrome_txs_annotation.Rds",
              "/scratch/fanslerm/test.merged.txs.rds")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 3) {
        stop("Incorrect number of arguments!\nUsage:\n> collapse_se_txs.R <sceRDS> <annots> <outFile>\n")
    }
}
arg.seRDS     <- args[1]
arg.annots    <- args[2]
arg.outFile   <- args[3]

## Load SCE
se <- readRDS(arg.seRDS)

## Load transcript annots
annots.txs <- readRDS(arg.annots)

## Collapse using merging info
mat.merged <- fac2sparse(rowData(se)$transcript_id_merged) %*% assay(se)

se.merged <- SummarizedExperiment(assays=list(counts=mat.merged),
                                  rowData=rowData(se)[rownames(mat.merged),1, drop=FALSE],
                                  colData=colData(se))

rowData(se.merged) <- merge(rowData(se.merged), annots.txs, by="transcript_id", sort=FALSE)

saveRDS(se.merged, arg.outFile)
