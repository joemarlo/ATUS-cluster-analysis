library(tidyverse)
library(mclust)
library(NbClust)
set.seed(44)
theme_set(theme_minimal())

atus_long <- read_tsv('Data/atus.tsv')
demographics <- read_tsv('Data/demographic.tsv')

# analysis ----------------------------------------------------------------

# check for NAs
dim(na.omit(atus_long)) == dim(atus_long)

# look at the data
atus_long %>% 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~description, scales = "free") + 
  labs(title = 'Feature densities (un-transformed)', 
       x = NULL,
       y = NULL)

# log transform the data
cats_to_cube <- unique(atus_long$description)[
  !(unique(atus_long$description) %in% c('Socializing, Relaxing, and Leisure', 'Sleep'))]
atus_cube_root <- atus_long %>%
  filter(description %in% cats_to_cube) %>% 
  group_by(description) %>% 
  mutate(value = value^(1/3)) %>%
  ungroup() %>% 
  bind_rows(
    atus_long %>%
      filter(!(description %in% cats_to_cube))
  )

atus_cube_root %>%
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap( ~ description, scales = "free") +
  labs(title = 'Feature densities (cube-root-transformed)',
       x = NULL,
       y = NULL)

# check spread of the data
atus_cube_root %>% 
  group_by(description) %>% 
  summarize(var = var(value))

# scale the data ??
atus_scaled <- atus_cube_root %>% 
  group_by(description) %>% 
  mutate(value = scale(value)) %>% 
  ungroup()

# check spread of the data
atus_scaled %>% 
  group_by(description) %>% 
  summarize(var = var(value))

# look at the data
atus_scaled %>% 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~description, scales = "free") + 
  labs(title = 'Feature densities (cube and scaled transformed)', 
       x = NULL)

# pivot wider and turn into matrix
atus_wide <- atus_scaled %>% 
  pivot_wider(values_from = value, names_from = description) %>% 
  select(-ID) %>% 
  as.matrix()


# PCA ---------------------------------------------------------------------

# run PCA
atus_pca <- prcomp(atus_wide)
summary(atus_pca)

pca_plot <- atus_pca$x[ , 1:3]
rgl::plot3d(pca_plot)


# resample data using survey weights --------------------------------------

# simple sample
# atus_only <- sample_n(atus_wide, 100000)

# function to scale [0.1]
scale_01 <- function(x) (x - min(x)) / (max(x) - min(x))

# sample using survey weights
total_rows <- nrow(atus_wide)
sample_size <- 10000
rows_to_keep <- sample(1:total_rows, size = sample_size, prob = scale_01(demographics$survey_weight), replace = TRUE)
IDs_kept <- atus_scaled %>% 
  pivot_wider(values_from = value, names_from = description) %>% 
  select(ID) %>% 
  .[rows_to_keep,]
atus_resampled <- atus_wide[rows_to_keep,]

# pairs plot of resampled data
as_tibble(atus_resampled) %>% 
  GGally::ggpairs(mapping = aes(alpha = 0.2))

# correlations
corr_matrix <- cor(atus_resampled)
corr_matrix[lower.tri(corr_matrix)] %>% density() %>% plot(main = 'Density of correlations between activities')

# run PCA
atus_pca <- prcomp(atus_resampled)
summary(atus_pca)


# memory management -------------------------------------------------------

rm(atus_scaled, atus_wide, atus_cube_root, corr_matrix, pca_plot, cats_to_cube, total_rows)
gc()


# hierarchical -------------------------------------------------------------

# distance matrix for features
dist_sc <- dist(atus_resampled, method = 'euclidean')
# try single, centroid, and ward (D2) linkage hier clustering
# hcl_single <- hclust(d = dist_sc, method = 'single')
# hcl_centroid <- hclust(d = dist_sc, method = 'centroid')
hcl_ward <- hclust(d = dist_sc, method = 'ward.D2')

library(dendextend)

dev.off()
# par(mfrow = c(3, 1))
# # nearest neighbors method
# plot(hcl_single, hang = -1, main = 'Single Linkage', 
#      labels = FALSE, xlab = '', sub = '')
# # groups centroid
# plot(hcl_centroid, hang = -1, main = 'Centroid Linkage', 
#      labels = FALSE, xlab = '',  sub = '')
# Ward’s minimum variance method, 
# with dissimilarities are squared before clustering
dend <- as.dendrogram(hcl_ward)
hcl_k <- 3
dend_col <- color_branches(dend, k = hcl_k)
plot(dend_col, main = paste0('Ward (D2) Linkage: K = ', hcl_k))

# examine the demographics through the clusters
groupings <- cutree(hcl_ward, 3)
demographics %>% 
  right_join(bind_cols(IDs_kept, group = groupings)) %>%
  select(-c('ID', 'survey_weight')) %>%
  mutate_all(as.character) %>% 
  pivot_longer(cols = -group) %>% 
  ggplot(aes(x = value)) +
  geom_histogram(stat = 'count') +
  facet_grid(group ~ name, scales = 'free') +
  labs(title = 'Demographics split by cluster', 
       x = NULL,
       y = NULL)
  

# memory management -------------------------------------------------------

rm(dist_sc, dend, dend_col, hcl_k, dend_col, groupings, sample_size)
gc()


# optimizing hierarchical cluster sizes ------------------------------------------------

# get optimal cluster sizes 
hcl_ch <- NbClust(
  data = atus_resampled,
  max.nc = 14,
  method = 'ward.D2',
  index = 'ch'
)
hcl_ch$All.index %>% plot(type = 'l', x = 2:14, main = 'Hierarchical: C(g) index', xlab = 'n clusters')

# get optimal cluster sizes 
hcl_silhouette <- NbClust(
  data = atus_resampled,
  max.nc = 14,
  method = 'ward.D2',
  index = 'silhouette'
)
hcl_silhouette$All.index %>% plot(type = 'l', x = 2:14, main = 'Hierarchical: Silhouette', xlab = 'n clusters')



# kmeans ------------------------------------------------------------------

# get optimal cluster sizes 
km_ch <- NbClust(
  data = atus_resampled,
  max.nc = 20,
  method = 'kmeans',
  index = 'ch'
)
km_ch$All.index %>% plot(type = 'l', x = 2:20, main = 'Kmeans: C(g) index', xlab = 'n clusters')

# get optimal cluster sizes 
km_silhouette <- NbClust(
  data = atus_resampled,
  max.nc = 20,
  method = 'kmeans',
  index = 'silhouette'
)
km_silhouette$All.index %>% plot(type = 'l', x = 2:20, main = 'Kmeans: Silhouette', xlab = 'n clusters')

# run the final kmeans algo with optimal number of clusters
km_eleven <- kmeans(x = atus_resampled,
                  centers = km_ch$Best.nc[[1]],
                  nstart = 100,
                  iter.max = 30,
                  algorithm = 'Hartigan-Wong')

# plot the clusters in PC space
atus_pca$x[, 1:2] %>% 
  as_tibble() %>% 
  mutate(Cluster = as.factor(km_eleven$cluster)) %>% 
  ggplot(aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.3) +
  labs(main = 'K-means cluster solution in PC space')

# examine the demographics through the clusters
demographics %>% 
  right_join(bind_cols(IDs_kept, cluster = as.factor(km_eleven$cluster))) %>%
  select(-c('ID', 'survey_weight')) %>%
  mutate_all(as.character) %>% 
  pivot_longer(cols = -cluster) %>% 
  ggplot(aes(x = value)) +
  geom_histogram(stat = 'count') +
  facet_grid(cluster ~ name, scales = 'free') +
  labs(title = 'Demographics split by cluster', 
       subtitle = 'K-means cluster solution',
       x = NULL,
       y = NULL)


# model based  clustering -------------------------------------------------

# run model
mcl <- Mclust(atus_resampled)
summary(mcl)

# plot 
plot(mcl, what = "BIC")
factoextra::fviz_mclust(mcl, "classification", geom = "point", alpha = 0.3)

# examine the demographics through the clusters
demographics %>% 
  right_join(bind_cols(IDs_kept, cluster = as.factor(mcl$classification))) %>%
  select(-c('ID', 'survey_weight')) %>%
  mutate_all(as.character) %>% 
  pivot_longer(cols = -cluster) %>% 
  ggplot(aes(x = value)) +
  geom_histogram(stat = 'count') +
  facet_grid(cluster ~ name, scales = 'free') +
  labs(title = 'Demographics split by cluster', 
       subtitle = 'Model-based cluster solution',
       x = NULL,
       y = NULL)

