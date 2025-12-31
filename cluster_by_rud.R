####*************************RUD-Based Cell Clustering***********************************#########
library(movAPA)

# Load annotated PAC dataset
scPACds <- readRDS("PAC_RUD_INPUT_DIR/sierra.merge_trap.peak.PAC_RUD.rds")

# Extract 3'UTR APA site pairs
apa3UTR <- movAPA::get3UTRAPAds(scPACds[scPACds@anno$ftr == "3UTR"])
apa_pairs <- get3UTRAPApd(apa3UTR, minDist = 20, maxDist = 10000, minRatio = 0.05, fixDistal = FALSE, addCols = 'pd')

# Compute smartRUD index
RUD <- movAPAindex(apa_pairs, method = "smartRUD", sRUD.oweight = FALSE, clearPAT = 1)
RUD[is.nan(RUD)] <- 0
saveRDS(RUD, file = "RUD_OUTPUT_DIR/RUD.RData")

# Build Seurat object
metadata <- read.table("METADATA_DIR/meta.txt", header = TRUE, sep = "\t")
rownames(metadata) <- metadata$cell

rud_seurat <- CreateSeuratObject(counts = RUD, project = "RUD")
rud_seurat <- rud_seurat[, metadata$cell]

idx <- match(colnames(rud_seurat), metadata$cell)
rud_seurat@meta.data$sample     <- metadata$sample[idx]
rud_seurat@meta.data$diagnosis  <- metadata$diagnosis[idx]
rud_seurat@meta.data$cluster    <- metadata$cluster[idx]

# Seurat pipeline
DefaultAssay(rud_seurat) <- "RNA"
rud_seurat <- NormalizeData(rud_seurat)
rud_seurat <- FindVariableFeatures(rud_seurat, nfeatures = 2000)
rud_seurat <- ScaleData(rud_seurat)
rud_seurat <- RunPCA(rud_seurat, npcs = 10)
rud_seurat <- harmony::RunHarmony(rud_seurat, group.by.vars = "sample")
rud_seurat <- RunUMAP(rud_seurat, reduction = "harmony", dims = 1:10)
rud_seurat <- FindNeighbors(rud_seurat, reduction = "harmony", dims = 1:10)
rud_seurat <- FindClusters(rud_seurat, resolution = 0.3)

rud_seurat@meta.data$celltype <- metadata[colnames(rud_seurat), ]$cluster
save(rud_seurat, file = "RUD_OUTPUT_DIR/RUD_cluster.RData")