library(tidyverse)
library(Seurat)
library(lineup)

nth.euclidean.distance <- function(a, b, n = 50) {
  sq.dif <- function(x, y) {
    z <- (x - y)**2
    return(z)
  }

  a <- a[1:n]
  b <- b[1:n]

  c <- unlist(map2(a, b, sq.dif))

  dist <- sqrt(sum(c))
  return(dist)
}

# Data not provided in repository
obj <- readRDS('all_cells.rds')
sathyamurthy.data <- subset(obj, subset = dataset == 'Sathyamurthy')
haring.data <- subset(obj, subset = dataset == 'Haring')

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

# Reduce to neural clusters only
sathyamurthy.data <- subset(sathyamurthy.data, subset = orig_reduced %in% clusters[8:37])

# Get embeddings
sathyamurthy.pca <- as_tibble(
    Embeddings(sathyamurthy.data, reduc = "integrated_PCA")
  ) %>%
  add_column(cluster = as.character(
    sathyamurthy.data[["orig_sathyamurthy_clusters"]][[1]]
  )) %>%
  pivot_longer(-cluster, names_to = "PC", values_to = "coord") %>%
  mutate_at("PC", str_remove_all, pattern = "\\D") %>%
  mutate_at("PC", as.numeric) %>%
  group_by(cluster, PC) %>%
  summarise(avg.pc = mean(coord)) %>%
  ungroup() %>%
  pivot_wider(names_from = cluster, values_from = avg.pc) %>%
  select(-PC) %>%
  select(`DE-1`, `DE-2`:`DE-9`, `DE-10`:`DE-16`, `DI-1`:`MI-4`) %>%
  slice(1:50)

haring.pca <- as_tibble(Embeddings(haring.data, reduc = "pca")) %>%
  add_column(cluster = haring.data[["orig_haring_clusters"]][[1]]) %>%
  pivot_longer(-cluster, names_to = "PC", values_to = "coord") %>%
  mutate_at("PC", str_remove_all, pattern = "\\D") %>%
  mutate_at("PC", as.numeric) %>%
  group_by(cluster, PC) %>%
  summarise(avg.pc = mean(coord)) %>%
  ungroup() %>%
  pivot_wider(names_from = cluster, values_from = avg.pc) %>%
  select(-PC) %>%
  select(Glut1, Glut2:Glut9, Glut10:Glut15, Gaba1,Gaba2:Gaba9, Gaba10:Gaba15) %>%
  slice(1:50)

# Calculate distance
# Custom function used to allow calculations between matrices
distance <- as_tibble(
  outer(anu.pca, haring.pca, 
    Vectorize(nth.euclidean.distance, vectorize.args = c("a", "b"))
  ), 
  rownames = "cluster")

# Get HVGs
hvg.500 <- VariableFeatures(obj)[1:500]

# Get Expression
sathyamurthy.exp <- AverageExpression(
  sathyamurthy.data, assays = 'RNA', features = hvg.500
)[[1]]
haring.exp <- AverageExpression(
  haring.data, assays = 'RNA', features = hvg.500
)[[1]]

# Get expression correlation
correlation <- as_tibble(
  corbetw2mat(sathyamurthy.exp, haring.expression, what = "all"),
  rownames = 'clusters'
)
