#============================================Dataset Splitting============================================

# Split the expression data for each cell type into training and test sets with a 7:3 ratio, repeated 10 times
# This results in 10 pairs of training and test datasets

combat_L23 <- readRDS("combat_L23.Rdata")
n_splits <- 10

train_list <- list()
test_list <- list()
sample_list <- list()
test_sample_list <- list()

for (i in 1:n_splits) {
  set.seed(i)
  train_id <- sample(seq_len(ncol(combat_L23 )),size=floor(0.7*ncol(combat_L23)))
  L23_train <- combat_L23[,train_id]
  L23_train_sample <- sample_combat[match(colnames(L23_train),sample_combat$cell),]
  L23_test <- combat_L23[,-train_id]
  L23_test_sample <- sample_combat[match(colnames(L23_test),sample_combat$cell),]
  
  train_list[[i]] <- L23_train
  names(train_list)[i] <- paste0("IN_PV_train_", i)
  
  test_list[[i]] <- L23_test
  names(test_list)[i] <- paste0("IN_PV_test_", i)
  
  sample_list[[i]] <- L23_train_sample
  names(sample_list)[i] <- paste0("IN_PV_train_sample_", i)
  
  test_sample_list[[i]] <- L23_test_sample
  names(test_sample_list)[i] <- paste0("IN_PV_test_sample_", i)
  
  write.csv(train_list[[i]], file=paste0(names(train_list)[i], ".csv"))
  write.csv(test_list[[i]], file=paste0(names(test_list)[i], ".csv"))
  write.csv(sample_list[[i]], file=paste0(names(sample_list)[i], ".csv"))
  write.csv(test_sample_list[[i]], file=paste0(names(test_sample_list)[i], ".csv"))
}