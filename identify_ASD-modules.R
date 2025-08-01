library(RcppCNPy)

#==============================Identify ASD-Related Gene Modules=======================================

#### t-test analysis####

t_test_result_list <- list()
sp_result_list <- list()

for (i in 1:5) {
  Wmap <- npyLoad(paste0("L23_W_array", i, ".npy"))
  train_sample <- read.csv(paste0("L23_train_sample_", i, ".csv"))
  
  diagnosis<-table(train_sample$diagnosis)
  
  index.control <- which(train_sample$diagnosis == "Control")
  
  W.case <- Wmap[,-index.control]
  W.control <- Wmap[,index.control]
  dim(W.case)
  dim(W.control)
  
  T.TEST.p.value.data <- data.frame(module.ID=NULL, p.value=NULL)
  
  for(j in seq_len(nrow(Wmap))) {
    pvalue <- t.test(W.case[j,], W.control[j,])$p.value
    
    temp <- data.frame(module.ID=j, p.value=pvalue)
    
    T.TEST.p.value.data <- rbind(T.TEST.p.value.data, temp)
  }
  
  # Adjust p-values using Benjamini-Hochberg (BH) method
  T.TEST.p.value.data$p.adjust <- p.adjust(T.TEST.p.value.data$p.value, "BH", n=nrow(T.TEST.p.value.data))

  T.summary.p.value <- summary(T.TEST.p.value.data$p.value)
  T.summary.p.adjust <- summary(T.TEST.p.value.data$p.adjust)
  
  # Identify significant modules using raw p-value threshold of 0.05
  T.index <- which(T.TEST.p.value.data$p.value <= 0.05 )
  
  rownames(Wmap) <- paste0("M",1:nrow(Wmap))
  colnames(Wmap) <- paste0("S",1:ncol(Wmap))
  
  t_test_result_list[[i]] <- list(
    diagnosis = diagnosis,
    T.summary.p.value = T.summary.p.value,
    T.summary.p.adjust = T.summary.p.adjust,
    T.TEST.p.value.data = T.TEST.p.value.data,
    T.significant_modules = rownames(Wmap)[T.index]
  )
  
  
####Spearman correlation test####

  # Create a vector Y: 1 for ASD (case), -1 for Control
  Y <- rep(1,time=ncol(Wmap))
  Y[which(train_sample$diagnosis=="Control")] <- -1
  table(Y)
  
  sp.data <- data.frame(NULL)
  
  for(k in 1:nrow(Wmap)){
    temp <- cor.test(Y,Wmap[k,],method="spearman")
    sp.data.temp <- data.frame(rho=temp$estimate,
                               p.value=temp$p.value,
                               module.ID=rownames(Wmap)[k])
    
    sp.data <- rbind(sp.data ,sp.data.temp)
  }
  
  # Adjust p-values using Benjamini-Hochberg (BH) method
  sp.data$p.adjust <- p.adjust(sp.data$p.value,"BH")  
  
  sp.summary.p.value <- summary(sp.data$p.value)
  sp.summary.p.adjust <- summary(sp.data$p.adjust)
  sp.summary.rho <- summary(abs(sp.data$rho))
  
  
  # Identify significant modules using raw p-value threshold of 0.05
  sp.index <- which(sp.data$p.value <= 0.05 )
  
  sp_result_list[[i]] <- list(
    sp.summary.p.value = sp.summary.p.value,
    sp.summary.p.adjust = sp.summary.p.adjust,
    sp.summary.rho = sp.summary.rho,
    sp.data = sp.data,
    sp.significant_modules = rownames(Wmap)[sp.index]
  )
}

####Find intersection of Wilcoxon and Spearman significant modules to define ASD-associated gene modules####
for (i in 1:10) {
  sig_module <- intersect(wilcox_test_result_list[[i]][["Wilcox.significant_modules"]],sp_result_list[[i]][["sp.significant_modules"]])
  assign(paste0("sig_module_", i), sig_module)
}

#L23_train_1 <- readRDS("G:/sierra_rud/RUD.Rdata")
#Load SMAF module matrix data (Umap)

for (i in 1:10) {
  filename <- paste0("SMAF_Umap_train_", i, ".npy")
  
  np_array <- import("numpy")$load(filename)
  
  umap_df <- as.data.frame(np_array)
  
  rownames(umap_df) <- rownames(L23_train_1)
  colnames(umap_df) <- paste0(1:ncol(umap_df))
  
  assign(paste0("Umap_", i), umap_df)
}

for (i in 1:10) {
  umap_var <- get(paste0("Umap_", i))
  col_select <- get(paste0("sig_module_", i))
  
  assign(paste0("Umap_", i), umap_var[, col_select,drop = FALSE])
}

for (i in 1:10) {
  umap_df <- get(paste0("Umap_", i))
  
  colnames(umap_df) <- paste0(i, "_", colnames(umap_df))
  
  assign(paste0("Umap_", i), umap_df) 
}


Neu_NRGN_II <- cbind(Umap_1,Umap_2,Umap_3,Umap_4,Umap_5,Umap_6,Umap_7,Umap_8,Umap_9,Umap_10)

colnames(Neu_NRGN_II) <- paste0("Neu_NRGN_II_", colnames(Neu_NRGN_II))

ASDsig_module <- cbind(AST_FB,AST_PP,Endothelial,IN_PV,IN_SST,IN_SV2C,IN_VIP,L4,L23,L56,L56_CC,Microglia,Neu_mat,Neu_NRGN_I,Neu_NRGN_II
                       ,Oligodendrocytes,OPC)
gene_module <- ASDsig_module


####Compute Z-scores for Genes in Each Module####
z_scores <- apply(gene_module, 2, function(x) {
  (x - mean(x))/sd(x)
})

threshold <- 1

selected_genes <- vector("list", ncol(gene_module))

# Iteratively adjust threshold
repeat {
  too_many_genes <- FALSE 
  
  for (i in 1:ncol(gene_module)) {
    module_genes <- rownames(gene_module)[which(gene_module[, i] != 0)] 
    zscores <- z_scores[module_genes, i] 
    
    selected <- module_genes[zscores > threshold]
    
    # If too many genes selected, increase threshold
    if (length(selected) > 1000) {
      too_many_genes <- TRUE
      break
    }
    
    #  Store selected genes, sorted by Z-score (descending)
    selected_genes[[i]] <- selected[order(zscores[selected], decreasing = TRUE)]
  }
  
  # Exit loop if all modules meet size criteria
  if (!too_many_genes) {
    break
  }
  
  # Increase threshold to reduce gene count
  threshold <- threshold + 0.05
}

# Output final threshold
cat("Final Z-score threshold:", threshold, "\n")

names(selected_genes) <- colnames(wil_sp_module)


####Calculate Module Recurrence Rate###

AST_FB_modules <- selected_genes[grep("^AST_FB_", names(selected_genes))]
AST_PP_modules <- selected_genes[grep("^AST_PP_", names(selected_genes))]
Endothelial_modules <- selected_genes[grep("^Endothelial_", names(selected_genes))]
IN_PV_modules <- selected_genes[grep("^IN_PV_", names(selected_genes))]

module_variables <- c(
  "AST_FB_modules", "AST_PP_modules", "Endothelial_modules",
  "IN_PV_modules", "IN_SST_modules", "IN_SV2C_modules", 
  "IN_VIP_modules", "L23_modules", "L4_modules", 
  "L56_modules", "L56_CC_modules", "Microglia_modules", 
  "Neu_mat_modules", "Neu_NRGN_I_modules", "Neu_NRGN_II_modules", 
  "Oligodendrocytes_modules"
)

result_list <- list()

for (module_var in module_variables) {
  current_module <- get(module_var)
  
  if (length(current_module) < 2) {
    cat("Skipping:", module_var, "- not enough modules to compare.\n")
    next
  }
  
  # Initialize data frame to store pairwise overlap rates
  module_overlap <- data.frame(
    Module_X = character(),
    Module_Y = character(),
    Overlap_Rate = numeric()
  )
  
  # Compute overlap rate for every pair of modules
  for (i in 1:(length(current_module) - 1)) {
    for (j in (i + 1):length(current_module)) {
      module_x <- current_module[[i]]  # 基因模块 X
      module_y <- current_module[[j]]  # 基因模块 Y
      
      # Count overlapping genes
      overlap <- length(intersect(module_x, module_y))
      
      # Compute overlap rate
      overlap_rate <- overlap / length(module_x)
      
      module_overlap <- rbind(
        module_overlap,
        data.frame(
          Module_X = names(current_module)[i],  
          Module_Y = names(current_module)[j],  
          Overlap_Rate = overlap_rate           
        )
      )
    }
  }
  
  # For each module, keep only the highest overlap with any other module
  max_overlap <- module_overlap %>%
    group_by(Module_X) %>%
    slice_max(Overlap_Rate, with_ties = FALSE) %>%
    ungroup()
  
  # Identify highly recurrent modules at different thresholds
  recurrence70modules <- max_overlap[max_overlap$Overlap_Rate >= 0.7, ]
  recurrence75modules <- max_overlap[max_overlap$Overlap_Rate >= 0.75, ]
  recurrence80modules <- max_overlap[max_overlap$Overlap_Rate >= 0.80, ]
  recurrence85modules <- max_overlap[max_overlap$Overlap_Rate >= 0.85, ]
  
  # Store results for this cell type
  result_list[[module_var]] <- list(
    overlap = max_overlap,
    recurrence70modules = recurrence70modules,
    recurrence75modules = recurrence75modules,
    recurrence80modules = recurrence80modules,
    recurrence85modules = recurrence85modules
  )
  
  cat("Processing completed for:", module_var, "\n")
}
