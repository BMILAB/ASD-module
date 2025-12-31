####*************************Seurat Clustering on PAS Counts (All vs. 3'UTR-Only)***********************************#########
library(Seurat); library(harmony); library(Matrix); library(dplyr)

# Load metadata
metadata <- read.table("METADATA_DIR/meta.txt", header = TRUE, sep = "\t")
rownames(metadata) <- metadata$cell
GSM_list <- paste0("GSM", 9262917:9262957)

# —————— All PAS clustering ——————
all_pas_list <- lapply(GSM_list, function(gsm) {
  ds <- readRDS(paste0(gsm, "_sierra.merged_trap.peak.PAC.rds"))
  CreateSeuratObject(counts = ds@counts, min.cells = 0, min.features = 0, project = "All_PAS")
})
names(all_pas_list) <- GSM_list

all_integrated <- merge(all_pas_list, y = all_pas_list[-1])
all_integrated <- all_integrated[, metadata$cell]

# —————— 3'UTR-only PAS clustering ——————
utr_pas_list <- lapply(GSM_list, function(gsm) {
  ds <- readRDS(paste0(gsm, "_sierra.merged_trap.peak.PAC.rds"))
  ds_3utr <- ds[ds@anno$ftr == "3UTR"]
  CreateSeuratObject(counts = ds_3utr@counts, min.cells = 0, min.features = 0, project = "3UTR_PAS")
})
names(utr_pas_list) <- GSM_list

utr_integrated <- merge(utr_pas_list, y = utr_pas_list[-1])
utr_integrated <- utr_integrated[, metadata$cell]

# —————— Add metadata and run clustering (using 3'UTR object as example) ——————
seurat_obj <- utr_integrated
idx <- match(colnames(seurat_obj), metadata$cell)
seurat_obj@meta.data$cluster    <- metadata$cluster[idx]
seurat_obj@meta.data$sample     <- metadata$sample[idx]
seurat_obj@meta.data$diagnosis  <- metadata$diagnosis[idx]

# Standard Seurat workflow
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 50)

# Harmony batch correction (by sample)
seurat_obj <- harmony::RunHarmony(seurat_obj, group.by.vars = "sample")
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)

# Assign cell types
celltype_colors <- c(
  "L2/3"="#96C3D8", "L5/6-CC"="#5D9BBE", "L4"="#F5B375", "L5/6"="#C0937E",
  "IN-VIP"="#67A59B", "IN-PV"="#A4D38E", "IN-SST"="#4A9D47", "IN-SV2C"="#F19294",
  "OPC"="#BDA7CB", "AST-PP"="#684797", "Oligodendrocytes"="#9983B7",
  "Neu-NRGN-I"="#CD9A99", "Neu-mat"="#DD4B52", "Endothelial"="#DA8F6F",
  "Neu-NRGN-II"="#F58135", "AST-FB"="#E45A5F", "Microglia"="#3477A9"
)

seurat_obj@meta.data$celltype <- metadata[colnames(seurat_obj), ]$cluster
save(seurat_obj, file = "CLUSTER_OUTPUT_DIR/utr_pA_cluster.RData")