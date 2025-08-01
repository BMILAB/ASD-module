library(scran)
library(sva)

#============================================Data Preprocessing============================================

# Remove nuclear and mitochondrial genes listed in MitoCarta2.0 human (1242 genes)
# MitoCarta2.0 human download: https://www.broadinstitute.org/files/shared/metabolism/mitocarta/human.mitocarta2.0.html
library(readxl)
mitocarta <- read_excel("Human.MitoCarta2.0.xls",sheet = 2)
index1 <- (which(genes$X1 %in% mitocarta$EnsemblGeneID))

grep("^MT",genes$X2[index1])
counts <- counts[-index1,]
genes <- genes[-index1,]
dim(counts);dim(genes)

# Select specific cell types for analysis due to computational constraints
# Analysis will be performed per cell type, so subset selection does not affect downstream steps
unique(sample$cluster)
index <- which(sample$cluster %in% c("AST-FB","Oligodendrocytes","AST-PP","OPC"))
length(index)
sample <- sample[index,]
counts <- counts[,index]
dim(counts)
unique(sample$cluster)
identical(colnames(counts),sample$cell) 

# Save processed data
saveRDS(counts,".Rdata")
write.csv(counts,"RUD_L23.csv")
saveRDS(sample,"RUD_L23_sample.Rdata")
write.csv(sample,"RUD_L23_sample.csv")


# Perform quality control using the scater package
# Reference: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scater/scater_01_qc.html

# Construct a SingleCellExperiment object
sce <- SingleCellExperiment(list(counts=counts),
                            colData=DataFrame(sample),
                            rowData=DataFrame(label=genes$X2),
                            metadata=list(study=sample))

sce <- SingleCellExperiment(list(counts=counts),
                            colData=DataFrame(sample),
                            rowData=DataFrame(label=genes),
                            metadata=list(study=sample))

head(colData(sce))
rm(counts)

# Remove lowly expressed genes and cells
dim(sce)
selected_c <- colnames(sce)[colSums(counts(sce)) > 1000]
selected_f <- rownames(sce)[rowSums(counts(sce)) > 20]
sce<- sce[selected_f, selected_c]
dim(sce)

# Remove spike-in genes (ERCC)
# Reference: https://www.singlecellcourse.org/basic-quality-control-qc-and-exploration-of-scrna-seq-datasets.html#dealing-with-confounders

# Assign ERCC spike-ins to altExp
altExp(sce,"ERCC") <- sce[grep("^ERCC",rownames(sce)), ]
sce <- sce[grep("^ERCC",rownames(sce),invert = T), ]
# Identify mitochondrial protein-coding genes
is.mito <-rownames(sce) %in% grep("^MT-",rownames(sce),value=T)
table(is.mito)

# Compute per-cell and per-feature QC metrics
# Includes total counts (library size), number of detected genes, mitochondrial counts, and % mitochondrial reads
sce_cell <- perCellQCMetrics(sce,subsets=list(Mito=is.mito))
sce_feature <- perFeatureQCMetrics(sce)

head(sce_cell)
head(sce_feature)
# - sum: total UMI counts (library size)
# - detected: number of expressed genes
# - subsets_Mito_percent: percentage of mitochondrial reads
# - altexps_ERCC_percent: percentage of ERCC spike-in reads

# Add QC metrics to the SingleCellExperiment object
sce <- addPerCellQC(sce, subsets=list(Mito=is.mito))
sce <- addPerFeatureQC(sce)

# Use quickPerCellQC to automatically filter out low-quality cells
reasons <- quickPerCellQC(sce_cell, sub.fields=c("subsets_Mito_percent", 
                                                 "altexps_ERCC_percent"))

# Add a metadata column indicating whether each cell is discarded
colSums(as.matrix(reasons))
table(reasons$discard)
sce$discard <- reasons$discard

sce_filter <- sce[,!colData(sce)$discard]
dim(sce_filter)

rm(sce)

# Data normalization and batch effect correction

# Perform initial clustering to avoid combining highly divergent cell types
# Compute size factors using scran to normalize for cell-specific biases
cluster <- quickCluster(sce_filter)
sce_filter <- computeSumFactors(sce_filter,clusters=cluster)
summary(sizeFactors(sce_filter))
sce_filter  <- logNormCounts(sce_filter )

# Merge multiple batches and correct batch effects using ComBat (from sva)
modcombat <- model.matrix(~as.factor(diagnosis) +as.factor(cluster), 
                          data=as.data.frame(colData(sce_filter )))
colData(sce_filter)$age_5 <- cut(sce_filter$age,5)
colData(sce_filter)$PMI_5 <- cut(sce_filter$post.mortem.interval..hours.,5)
colData(sce_filter)$RIN_3 <- cut(sce_filter $RNA.Integrity.Number,3)
batch <- paste(sce_filter$Seqbatch,
               sce_filter$sex,
               sce_filter$age_5,
               sce_filter$PMI_5,
               sce_filter$RIN_3,sep="-")
table(batch)
length(unique(batch))

# Convert log-normalized counts to matrix for ComBat
norm <- as.matrix(logcounts(sce_filter))
# removes 0 variance
data_filted <- norm[rowVars(norm)>0,]
dim(norm);dim(data_filted)

combat_edata <- ComBat(dat= data_filted,
                       batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

head(combat_edata)
dim(combat_edata)
class(combat_edata)

# Check for negative values introduced by ComBat
length(which(combat_edata[,1]<0))
# Replace negative values with zero
rm(norm,data_filted)
combat_edata <- ifelse(combat_edata<0,0,combat_edata)
summary(combat_edata[,1])
head(sample$cell)

# Extract metadata for cells in the ComBat-processed dataset
index <- match(colnames(combat_edata),sample$cell)
sample_combat <- sample[index,]
saveRDS(sample_combat,"sample_combat.Rdata")

# Extract data for L2/3 excitatory neurons
identical(colnames(combat_edata),sample_combat$sample_ID)
table(sample_combat$cluster)
L23.index <- which(sample_combat$cluster=="L2/3")
combat_L23 <- combat_edata[,L23.index]
dim(combat_L23)
saveRDS(combat_L23, "combat_L23.Rdata")