library(Seurat)
library(tidyverse)

out.path <- "./"

# First neural tier
# Data not provided in repository
# The reference is formed from all neurons that were in the coarse cell types reference
# The query is all remaining neurons
message("Tier 1")
message("\tReading data...")
reference <- readRDS('clean_neuro_train.rds')
message("\tReference data:")
print(reference)
query <- readRDS("clean_neuro_test.rds")
message("\tQuery data:")
print(query)

Idents(reference) <- "coarse_neuron_clusters"

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
  dims = 1:100
)
message("\tTransferring labels...")
predictions <- TransferData(
  anchorset = anchors,
  refdata = reference$coarse_neuron_clusters,
  weight.reduction = "pcaproject",
  dims = 1:100
)
query <- AddMetaData(query, metadata = predictions)

# Save queried objects
message("\tSaving coarse neuron results")
saveRDS(query, paste0(out.path, 'coarse_neuron_predictions.rds'))

# Second neural tier
# Data not provided in repository
# The reference is all mid/ventral neurons that were in the coarse cell types reference
# The query is all cells predicted to be mid/ventral in the first tier
reference <- readRDS('clean_midventral_train.rds')
query.ventral <- subset(query, subset = predicted.id =='Mid/Ventral')

Idents(reference) <- "final_cluster_assignment"

message("\tLearning anchors...")
anchors <- FindTransferAnchors(
  reference = reference,
  query = query.ventral,
  normalization.method = "LogNormalize",
  reference.assay = "integrated",
  query.assay = "integrated",
  reduction = "pcaproject",
  features = VariableFeatures(object = reference),
  npcs = NULL,
  dims = 1:30
)
message("\tTransferring labels...")
predictions <- TransferData(
  anchorset = anchors,
  refdata = reference$final_cluster_assignment,
  weight.reduction = "pcaproject",
  dims = 1:30
)
query.ventral <- AddMetaData(query.ventral, metadata = predictions)

# Save queried objects
message("\tSaving Mid/Ventral results...")
saveRDS(query.neural, paste0(out.path, "mid_ventral_predictions.rds"))
