#!/usr/bin/env Rscript

library(SummarizedExperiment)
library(Matrix)
library(readr)
library(dplyr)
library(stringr)

if (interactive()) {
    args <- c("data/se/hspcs_bulk.txs.raw.rds",
              "/data/mayrc/data/mca/gff/adult.utrome.e3.t200.f0.999.w500.merge.tsv",
              "/data/mayrc/data/atlas-mm/utrs/utr-metadata.tsv",
              "/scratch/fanslerm/test.merged.txs.rds")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 4) {
        stop("Incorrect number of arguments!\nUsage:\n> collapse_se_txs.R <sceRDS> <txMapFile> <annots> <outFile>\n")
    }
}
arg.seRDS     <- args[1]
arg.txMapFile <- args[2]
arg.annots    <- args[3]
arg.outFile   <- args[4]

## Load SCE
se <- readRDS(arg.seRDS)

## Load transcript map
tx2tx <- read_tsv(arg.txMapFile, skip=1, col_types='cc-',
                  col_names=c("transcript_name", "merged_transcript_name"))

## Load transcript annots
annots.txs <- read_tsv(arg.annots, col_types='c---------l--d-d-dd-d-d-c')


## Attach mapping data to matrix
rowData(se) <- merge(rowData(se), tx2tx, by="transcript_name", sort=FALSE)

## keep note of gene symbol mapping
tx2gene <- rowData(se) %>%
    as.data.frame() %>%
    select(transcript_name, merged_transcript_name, gene_symbol) %>%
    unique()

## WARNING! There are some transcripts that are crossing over into other genes!
## We will simply discard these for now...
bad.txs <- tx2gene %>%
    mutate(gene=str_match(merged_transcript_name, "^(.*)\\.\\d+$")[,2]) %>%
    filter(gene_symbol != gene)

print("[Warning] Removing the following transcripts for mismatching gene names:")
print(bad.txs)

bad.idx <- which(rownames(se) %in% bad.txs$transcript_name)

## Collapse using merging info
mat.merged <- fac2sparse(rowData(se[-bad.idx,])$merged_transcript_name) %*% assay(se[-bad.idx,])

se.merged <- SummarizedExperiment(assays=list(counts=mat.merged),
                                   rowData=rowData(se)[rownames(mat.merged),1:4],
                                  colData=colData(se))

rowData(se.merged) <- merge(rowData(se.merged), annots.txs, by="transcript_name", sort=FALSE)

saveRDS(se.merged, arg.outFile)
