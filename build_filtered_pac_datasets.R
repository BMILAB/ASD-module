####*************************Construct and Annotate Per-Sample PACdataset Objects with Endogenous PAS Filtering***********************************#########
library(Sierra); library(movAPA); library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Rsamtools)

# Setup
GSM_list <- paste0("GSM", 9262917:9262957)
cell_prefixes <- paste0(GSM_list, "_")
known_pas <- readRDS("REFERENCE_DIR/Human_KnownPASs_fourDBs.rds")
gff_annotation <- parseGenomeAnnotation("REFERENCE_DIR/gencode.v44.annotation.gtf")

pacount_list <- list()

for (i in seq_along(GSM_list)) {
  cat("Processing", i, "/", length(GSM_list), GSM_list[i], "\n")
  
  count_dir <- file.path("COUNT_INPUT_DIR", paste0(GSM_list[i], "_trap.covert.sierra.merged_peak.count"))
  annotation_file <- "ANNOTATION_INPUT_DIR/trap.peak_covert_sierra.merged.peak_annotations.txt"
  
  count_matrix <- as.data.frame(ReadPeakCounts(
    data.dir = count_dir,
    mm.file = "matrix.mtx.gz",
    barcodes.file = "barcodes.tsv.gz"
  ))
  peak_anno <- read.delim(annotation_file)
  
  # Align rows
  common_peaks <- intersect(rownames(count_matrix), rownames(peak_anno))
  count_matrix <- count_matrix[common_peaks, , drop = FALSE]
  peak_anno <- peak_anno[common_peaks, , drop = FALSE]
  
  # Build coordinates and IDs
  peak_anno$coord <- ifelse(peak_anno$strand == "+", peak_anno$end, peak_anno$start)
  peak_anno$sierra.PA_id <- rownames(peak_anno)
  
  anno_df <- peak_anno[, c("seqnames", "start", "end", "strand", "coord", "sierra.PA_id")]
  colnames(anno_df) <- c("chr", "start", "end", "strand", "coord", "sierra.PA_id")
  anno_df <- unique(anno_df)
  count_matrix <- count_matrix[rownames(anno_df), , drop = FALSE]
  
  # Create PACdataset
  scPACds <- readPACds(pacFile = anno_df, noIntergenic = FALSE, PAname = "PA")
  coldata <- data.frame(barcode = colnames(count_matrix), row.names = colnames(count_matrix))
  scPACds@colData <- coldata
  
  # Remove internal priming artifacts
  scPACds@anno$chr <- paste0("chr", scPACds@anno$chr)
  scPACds <- removePACdsIP(scPACds, BSgenome.Hsapiens.UCSC.hg38, returnBoth = TRUE, up = -10, dn = 10, conA = 6, sepA = 7, chrCheck = TRUE)
  
  IP_ovp <- findOvpPACds(qryPACds = scPACds$ip, sbjPACds = known_pas, d = 50, qryMode = 'region', sbjMode = 'point', returnNonOvp = TRUE)
  scPACds <- rbind(scPACds$real, IP_ovp$ovp)
  
  # Annotate and extend 3'UTR
  scPACds <- annotatePAC(scPACds, gff_annotation)
  scPACds <- ext3UTRPACds(scPACds, ext3UTRlen = 2000, extFtr = '3UTR')
  
  # Update counts and metadata
  rownames(scPACds@anno) <- scPACds@anno$sierra.PA_id
  scPACds@counts <- as.matrix(count_matrix[scPACds@anno$sierra.PA_id, , drop = FALSE])
  rownames(scPACds@counts) <- scPACds@anno$sierra.PA_id
  
  colnames(scPACds@counts) <- paste0(cell_prefixes[i], colnames(scPACds@counts))
  rownames(scPACds@colData) <- colnames(scPACds@counts)
  
  pacount_list[[GSM_list[i]]] <- scPACds
  saveRDS(scPACds, file = paste0(GSM_list[i], "_sierra.merged_trap.peak.PAC.rds"))
}

saveRDS(pacount_list, file = "sierra.merge.pacount.list.rds")