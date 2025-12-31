####*************************Annotate Merged Peaks Using GTF***********************************#########
source("UTILS_DIR/Sierra.Annotate.R")

gtf_file <- "REFERENCE_DIR/gencode.v44.annotation.gtf"  # must include "chr" prefix
merged_peak_file <- "MERGED_PEAK_DIR/nochr_trap.peak_covert_sierra.merged.peaks.txt"

AnnotatePeaksFromGTF2(
  peak.sites.file = merged_peak_file,
  gtf.file = gtf_file,
  output.file = file.path("ANNOTATION_OUTPUT_DIR", "annote_merged_peak")
)