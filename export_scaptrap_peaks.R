####*************************Build PACdataset from scAPAtrap Results and Export Peak Files***********************************#########
setwd("PROJECT_DIR")
library(movAPA)
library(TxDb.Hsapiens.UCSC.hg38.knownGene, quietly = TRUE, verbose = FALSE)

# Define sample ID range
sample_ids <- 9262917:9262957

for (sid in sample_ids) {
  file_path <- file.path("INPUT_DIR", paste0(sid, "/scAPAtrapData.rda"))
  load(file_path)
  peaks <- scAPAtrapData
  
  # Create and filter PAC dataset
  PACds <- createPACdataset(counts = peaks$peaks.count, anno = peaks$peaks.meta)
  rm(peaks)
  PACds <- subsetPACds(PACds, totPACtag = 10, minExprConds = 10, verbose = TRUE)
  
  # Annotate with genomic features
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  PACds <- annotatePAC(PACds, txdb)
  PACds <- ext3UTRPACds(PACds, 1000)
  PACds <- PACds[!PACds@anno$ftr == "intron"]
  
  # Format peak table
  peak <- PACds@anno[, c("gene", "chr", "strand", "start", "end")]
  peak$strand <- ifelse(peak$strand == "-", "-1", "1")
  peak$polyA_ID <- paste0(peak$gene, ":", peak$chr, ":", peak$start, "-", peak$end, ":", peak$strand)
  colnames(peak) <- c("Gene", "Chr", "Strand", "Fit.start", "Fit.end", "polyA_ID")
  
  # Save peak file
  GSM_id <- paste0("GSM", sid)
  output_file <- file.path("OUTPUT_DIR", paste0(GSM_id, "_trap.peak.txt"))
  write.table(peak, file = output_file, quote = FALSE, row.names = FALSE, sep = "\t")
}