####*************************Merge Peaks Across Samples Using Sierra***********************************#########
library(Sierra)
source("UTILS_DIR/sierra_dataset_merging.R")

sample_ids <- 9262917:9262957
peak_files <- file.path("PEAK_INPUT_DIR", sprintf("GSM%d_trap.peak.txt", sample_ids))
identifiers <- sprintf("GSM%d", sample_ids)

peak_dataset_table <- data.frame(
  Peak_file = peak_files,
  Identifier = identifiers,
  stringsAsFactors = FALSE
)

MergePeakCoordinates2(
  peak_dataset_table,
  output.file = file.path("OUTPUT_DIR", "trap.peak_covert_sierra.merged.peaks.txt"),
  ncores = 4
)