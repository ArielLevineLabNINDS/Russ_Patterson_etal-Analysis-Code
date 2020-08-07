library(Seurat)
library(tidyverse)

# Data not provided in the repository.
obj <- readRDS('clean_neurons.rds')

centroids <- as_tibble(Embeddings(obj, reduction = 'integrated_PCA')) %>%
  add_column('cluster' = as.characer(
    obj[['final_cluster_assignment']][[1]]
  )) %>%
  pivot_longer(-cluster, names_to = 'PC', values_to = 'coord') %>%
  mutate_at('PC', str_remove_all, pattern = '\\D') %>%
  mutate_at('PC', as.numeric) %>%
  group_by(cluster, PC) %>%
  summarise('avg.pc' = mean(coord)) %>%
  ungroup() %>%
  pivot_wider(names_from = cluster, values_from = avg.pc) %>%
  select(-PC) %>%
  slice(1:50)

distance <- dist(t(centroids))
hc <- hclust(distance, method = 'complete')
plot(hc, hang = -1)
