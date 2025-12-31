####*************************Count UMI Support for Merged Peaks Across Samples***********************************#########
source("UTILS_DIR/count_polyA.R")

sample_ids <- 9262917:9262957
samples <- paste0("SRR", sample_ids)
gtf_file <- "REFERENCE_DIR/gencode.v44.annotation_nochr.gtf"
merged_peak_file <- "MERGED_PEAK_DIR/nochr_trap.peak_covert_sierra.merged.peaks.txt"

# Per-sample counting
for (sample in samples) {
  whitelist <- file.path("CELLRANGER_DIR", paste0("cellranger_", sample, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"))
  bamfile <- file.path("CELLRANGER_DIR", paste0("cellranger_", sample, "/outs/possorted_genome_bam.bam"))
  out_dir <- file.path("COUNT_OUTPUT_DIR", paste0("GSM", substring(sample, 4), "_trap.covert.sierra.merged_peak.count"))
  
  CountPeaks(
    peak.sites.file = merged_peak_file,
    gtf.file = gtf_file,
    bamfile = bamfile,
    whitelist.file = whitelist,
    output.dir = out_dir,
    countUMI = TRUE,
    ncores = 4
  )
}

# Optional: global counting (uncomment if needed)
# CountPeaks(
#   peak.sites.file = merged_peak_file,
#   gtf.file = gtf_file,
#   bamfile = bamfile,
#   whitelist.file = whitelist,
#   output.dir = "COUNT_OUTPUT_DIR",
#   countUMI = TRUE,
#   ncores = 4
# )