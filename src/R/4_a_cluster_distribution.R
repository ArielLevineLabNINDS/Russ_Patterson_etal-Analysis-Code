library(tidyverse)
library(Seurat)

# Data not provided in repository
obj <- readRDS('all_cells.rds')

sathyamurthy.data <- subset(obj, subset = dataset == 'Sathyamurthy')
haring.data <- subset(obj, subset = datset == 'Haring')

# Sort to desired order
Idents(sathyamurthy.data) <- 'orig_sathyamurthy_clusters'
clusters <- sort(levels(Idents(sathyamurthy.data)))
clusters <- clusters[c(37, 36, 29, 1, 27, 34, 35, 2, 9:16, 3:4, 17, 5:8, 18:26, 28, 30:33, 40:48, 38, 39, 49, 50)]
Idents(sathyamurthy.data) <- ordered(Idents(sathyamurthy.data), levels = clusters)

Idents(haring.data) <- 'orig_haring_clusters'

# Correct typo in original labels
clusters[19] <- 'DE-12'
levels(Idents(sathyamurthy.data)) <- clusters

# Store changes
sathyamurthy.data[['orig_sathyamurthy_clusters']] <- Idents(sathyamurthy.data)

# Count co-labels
sathyamurthy.data <- as_tibble(
  sathyamurthy.data[[
    c('final_cluster_assignment', 'orig_sathyamurthy_clusters')
  ]]
  ) %>%
  group_by(final_cluster_assignment, orig_sathyamurthy_clusters) %>%
  summarise('co_count' = n()) %>%
  ungroup() %>%
  group_by(orig_sathyamurthy_clusters) %>%
  mutate('orig_count' = sum(co_count)) %>%
  ungroup() %>%
  mutate('co_frequency' = co_count / orig_count)

haring.data <- as_tibble(
  haring.data[[
    c('final_cluster_assignment', 'orig_haring_clusters')
  ]]
  ) %>%
  group_by(final_cluster_assignment, orig_haring_clusters) %>%
  summarise('co_count' = n()) %>%
  ungroup() %>%
  group_by(orig_haring_clusters) %>%
  mutate('orig_count' = sum(co_count)) %>%
  ungroup() %>%
  mutate('co_frequency' = co_count / orig_count)
