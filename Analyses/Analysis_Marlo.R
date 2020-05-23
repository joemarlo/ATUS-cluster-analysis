library(tidyverse)
library(mclust)
library(NbClust)
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
cats_to_log <- unique(atus_scaled$description)[
  !(unique(atus_scaled$description) %in% c('Socializing, Relaxing, and Leisure', 'Sleep'))]
atus_log <- atus_long %>%
  filter(description %in% cats_to_log) %>% 
  group_by(description) %>% 
  mutate(value = log(jitter(value + 0.5, factor = 1))) %>% # added noise to account for log 0
  ungroup() %>% 
  bind_rows(
    atus_long %>%
      filter(!(description %in% cats_to_log))
  )

atus_log %>% 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~description, scales = "free") + 
  labs(title = 'Feature densities (un-transformed)', 
       x = NULL)

atus_cube_root <- atus_long %>%
  filter(description %in% cats_to_log) %>% 
  group_by(description) %>% 
  mutate(value = value^(1/3)) %>%
  ungroup() %>% 
  bind_rows(
    atus_long %>%
      filter(!(description %in% cats_to_log))
  )

atus_cube_root %>%
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap( ~ description, scales = "free") +
  labs(title = 'Feature densities (cube-transformed)',
       x = NULL)

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
  select(-TUCASEID) %>% 
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
rows_to_keep <- sample(1:total_rows, size = sample_size, prob = scale_01(weights$TUFNWGTP), replace = TRUE)

resampled_atus <- atus_wide[rows_to_keep,]

# pairs plot of resampled data
as_tibble(resampled_atus) %>% 
  GGally::ggpairs(mapping = aes(alpha = 0.2))

# run PCA
atus_pca <- prcomp(resampled_atus)
summary(atus_pca)


# memory management -------------------------------------------------------

rm(atus_long, atus_scaled, weights)
gc()


# hierarchical -------------------------------------------------------------

# distance matrix for features
dist_sc <- dist(resampled_atus, method = 'euclidean')
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
hcl_k <- 4
dend_col <- color_branches(dend, k = hcl_k)
plot(dend_col, main = paste0('Ward (D2) Linkage: K = ', hcl_k))


# memory management -------------------------------------------------------

rm(dend, hcl_k, dend_col)
gc()




# Naive Bayes and SVM -----------------------------------------------------




# clustering --------------------------------------------------------------

# get optimal cluster sizes 
cluster_sizes_hcl_ch <- NbClust(data = nba_feat_sc,
                                # it will likely be harder to interpret clusters
                                # past this amount
                                max.nc = 6,
                                method = 'ward.D2',
                                index = 'ch')

# get optimal cluster sizes 
cluster_sizes_hcl_s <- NbClust(data = nba_feat_sc,
                               # it will likely be harder to interpret clusters
                               # past this amount
                               max.nc = 6,
                               method = 'ward.D2',
                               index = 'silhouette')

