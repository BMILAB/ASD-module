rm(list=ls())
options(stringsAsFactors = FALSE)
getwd()
setwd()

library(Matrix)
library(readr)
library(reticulate)
library(Seurat)

#============================================Load data============================================

# ASD data download source: autism.cells.ucsc.edu
# Raw data includes three files: matrix.mtx, barcodes.tsv, genes.tsv

# Read cell metadata
sample <- read.table("data/meta.txt",sep="\t",header=TRUE)
# Read count matrix
counts <- readMM("data/matrix.mtx")

# Read genes.tsv
genes <- read_tsv("data/genes.tsv",col_names=FALSE)
gene_ids <- genes$X2
genes <- rownames(rud)
# Read barcodes.tsv
cells <- read_tsv("data/barcodes.tsv",col_names = FALSE)
cells_ids <- cells$X1

# Assign gene IDs as row names and cell IDs as column names to the count matrix
rownames(counts) <- gene_ids
colnames(counts) <- cells_ids
identical(cells_ids,sample$cell) 
dim(counts)

# View cell types and the number of cells in each type
table(sample$cluster)