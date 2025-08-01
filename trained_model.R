library(xgboost)

# Define all cell types
cell_types <- c("AST_FB", "AST_PP", "Endothelial", "IN_PV", "IN_SST", "IN_SV2C", "IN_VIP",
                "L23", "L4", "L56", "L56_CC", "Microglia", "Neu_mat", "Neu_NRGN_I", "Neu_NRGN_II",
                "Oligodendrocytes", "OPC")

# Initialize results list
results <- list()

# Loop through each module in xgboost_list
for (module in names(xgboost_list)) {
  
  # Match cell type by prefix
  matches <- sapply(cell_types, function(ct) substr(module, 1, nchar(ct)) == ct)
  cell_type <- names(matches)[matches]
  if (length(cell_type) == 0) {
    warning(paste("Unknown cell type for module:", module))
    next
  }
  cell_type <- cell_type[1]
  
  # Extract partition number (assuming module name ends with digits)
  partition_match <- regmatches(module, regexec("\\d+$", module))
  if (is.na(partition_match) || length(partition_match) == 0) {
    warning(paste("Cannot extract partition from module:", module))
    next
  }
  partition <- partition_match[1]
  
  # Construct file paths
  train_path <- paste0("F:/apa_rud/", cell_type, "/", cell_type, "_train_", partition, ".csv")
  test_path <- paste0("F:/apa_rud/", cell_type, "/", cell_type, "_test_", partition, ".csv")
  train_sample_path <- paste0("F:/apa_rud/", cell_type, "/", cell_type, "_train_sample_", partition, ".csv")
  test_sample_path <- paste0("F:/apa_rud/", cell_type, "/", cell_type, "_test_sample_", partition, ".csv")
  
  # Check if all required files exist
  if (!all(file.exists(c(train_path, test_path, train_sample_path, test_sample_path)))) {
    warning(paste("Missing files for module:", module))
    next
  }
  
  # Read training data
  train <- read.csv(train_path)
  rownames(train) <- make.unique(as.character(train$X))
  train <- train[, -1]  # Remove first column (assumed to be row index)
  
  # Read test data
  test <- read.csv(test_path)
  rownames(test) <- make.unique(as.character(test$X))
  test <- test[, -1]  # Remove first column
  
  # Read training sample labels
  train_sample <- read.csv(train_sample_path, row.names = 1)
  
  # Read test sample labels
  test_sample <- read.csv(test_sample_path, row.names = 1)
  
  # Extract features common to both the module and the dataset
  features <- intersect(xgboost_list[[module]], rownames(train))
  if (length(features) == 0) {
    warning(paste("No overlapping features for module:", module))
    next
  }
  
  # Subset data to selected features
  train <- train[features, ]
  test <- test[features, ]
  
  # Transpose data to have samples as rows and genes as columns
  train <- t(train)
  test <- t(test)
  
  # Encode diagnosis labels numerically: ASD -> 1, Control -> 0
  train_sample$diagnosis <- ifelse(train_sample$diagnosis == "ASD", 1, 0)
  test_sample$diagnosis <- ifelse(test_sample$diagnosis == "ASD", 1, 0)
  
  # Train XGBoost classifier
  xgb_model <- xgboost(data = as.matrix(train), 
                       label = train_sample$diagnosis,
                       max_depth = 6, 
                       eta = 0.01, 
                       subsample = 0.5, 
                       nrounds = 15,
                       verbose = 0)
  
  # Predict on test set
  pred_prob <- predict(xgb_model, as.matrix(test))
  pred <- ifelse(pred_prob > 0.5, 1, 0)
  
  # Calculate evaluation metrics
  tp <- sum(pred == 1 & test_sample$diagnosis == 1)  # True positives
  tn <- sum(pred == 0 & test_sample$diagnosis == 0)  # True negatives
  fp <- sum(pred == 1 & test_sample$diagnosis == 0)  # False positives
  fn <- sum(pred == 0 & test_sample$diagnosis == 1)  # False negatives
  
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))  # Same as sensitivity
  specificity <- ifelse(tn + fp == 0, 0, tn / (tn + fp))
  f1_score <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  
  # Store results for this module
  results[[module]] <- list(
    accuracy = accuracy,
    precision = precision,
    sensitivity = recall,
    specificity = specificity,
    f1_score = f1_score,
    n_features = length(features)  # Number of genes used
  )
}