library(Seurat)
library(tidyverse)

out.path <- "./"

# Coarse cell types tier
# Data not provided in repository
message("Coarse label transfer")
message("\tReading data...")

# For splitting data during development
# This is saved from the SVM split
val_idx <- py_load_object('val_idx.pkl')
val_idx <- val_idx

data <- readRDS("data/final_cluster_assignment.RDS")
reference <- data[, -val_idx]
message("\tReference data:")
print(reference)
query <- data[, val_idx]
message("\tQuery data:")
print(query)

Idents(reference) <- "coarse_clusters"

message("\tLearning anchors...")
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "LogNormalize",
  reference.assay = "integrated",
  query.assay = "integrated",
  reduction = "pcaproject",
  features = VariableFeatures(object = reference),
  npcs = NULL,
  dims = 1:28
)
message("\tTransferring labels...")
predictions <- TransferData(
  anchorset = anchors,
  refdata = reference$coarse_clusters,
  weight.reduction = "pcaproject",
  dims = 1:28
)
query <- AddMetaData(query, metadata = predictions)

# Save queried objects
message("\tSaving coarse results...")
saveRDS(query, paste0(out.path, 'predicted_coarse_types.rds'))

# For passing on to neural network in 5e, subset out predicted neurons and doublets
query.neural <- subset(query, subset = predicted.id %in% c('Neuron', 'Doublets', 'Motoneuron'))
saveRDS(query.neural, paste0(out.path, 'predicted_neurons_doublets.rds'))
