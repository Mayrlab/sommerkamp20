#!/usr/bin/env Rscript

library(SummarizedExperiment)
library(Matrix)
library(readr)
library(dplyr)

if (interactive()) {
    args <- c("data/se/hspcs_bulk.txs.rds",
              "/data/mayrc/data/atlas-mm/utrs/genes-utr-metadata-lengths.tsv",
              "/scratch/fanslerm/test.genes.rds")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 3) {
        stop("Incorrect number of arguments!\nUsage:\n> se_txs_to_genes.R <seRDS> <annots> <outFile>\n")
    }
}
arg.seRDS  <- args[1]
arg.annots  <- args[2]
arg.outFile <- args[3]

################################################################################
                                        # Load Data
################################################################################
se <- readRDS(arg.seRDS)

df.genes <- rowData(se) %>%
    as.data.frame() %>%
    group_by(gene_id, gene_symbol, chromosome) %>%
    summarise(transcripts=paste(transcript_name, collapse=';')) %>%
    as.data.frame()
rownames(df.genes) <- df.genes$gene_id

## Load gene annots
annots.genes <- read_tsv(arg.annots, col_types='-ccddcdd')

################################################################################
                                        # Gene Counts
################################################################################
mat.genes <- fac2sparse(rowData(se)$gene_id) %*% assay(se)

se.genes <- SummarizedExperiment(assays=list(counts=mat.genes),
                                  rowData=df.genes[rownames(mat.genes),],
                                  colData=colData(se))

rowData(se.genes) <- merge(rowData(se.genes), annots.genes, by="gene_id", sort=FALSE, all.x=TRUE)

################################################################################
                                        # Export Data
################################################################################
saveRDS(se.genes, arg.outFile)
