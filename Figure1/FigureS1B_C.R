## ===========================================
## K-means & GMM clustering for Visium and Visium HD data
## (micro vs macro metastasis identification)
## ===========================================

set.seed(1234)

## Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(matrixStats)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(stringr)
  library(plyr)
  library(factoextra)
  library(mclust)
})

## -------------------------------
## 0) Load data
## -------------------------------
load("liver_foci_mtx.RData")  # contains liver_foci_agg$Spatial
load("foci_dat.RData")        # contains annotation info

## Expression matrix: genes x samples
dats <- as.matrix(liver_foci_agg$Spatial)
message(dim(dats)[1], " genes × ", dim(dats)[2], " samples loaded")

## Create unique cell/spot ID
foci_dat$bar2 <- paste0(foci_dat$sample, "~", foci_dat$Barcode)

## -------------------------------
## 1) Select highly variable genes
## -------------------------------
vars <- rowVars(dats)
names(vars) <- rownames(dats)
top_n <- min(5000, nrow(dats))
hvg <- names(sort(vars, decreasing = TRUE))[seq_len(top_n)]
mat_hvg <- dats[hvg, , drop = FALSE]

## Log transform & transpose (samples x genes)
expr_log <- t(log2(mat_hvg + 1))

## PCA (top 20 PCs for clustering input)
pca <- prcomp(expr_log, center = TRUE, scale. = TRUE)
pcs <- as.data.frame(pca$x[, 1:20, drop = FALSE])
pcs$sample <- rownames(pcs)

## True labels (based on name pattern)
pcs$true_label <- ifelse(grepl("micro", pcs$sample, ignore.case = TRUE), "micro",
                         ifelse(grepl("macro", pcs$sample, ignore.case = TRUE), "macro", "NA"))

## -------------------------------
## 2) Add group size as covariate
## -------------------------------
ref_size <- as.data.frame(table(foci_dat$group))
colnames(ref_size) <- c("group", "size")
ref_size$group <- gsub("_", "-", ref_size$group)

pcs$group <- str_split(pcs$sample, "~", simplify = TRUE)[, 1]
pcs <- pcs %>%
  left_join(ref_size, by = "group") %>%
  mutate(size = ifelse(is.na(size), 0, size),
         size_log = log2(pmax(size, 1)))

## ===========================================
## 3) K-means clustering
## ===========================================
km_input <- pcs[, 1:20, drop = FALSE]
km <- kmeans(km_input, centers = 2, nstart = 30)
pcs$km_cluster <- factor(km$cluster)

## K-means visualization
p1 <- ggplot(pcs, aes(PC1, PC2, color = km_cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "K-means (k = 2)", color = "Cluster") +
  theme_classic()
print(p1)

## Visualize true labels
p2 <- ggplot(pcs, aes(PC1, PC2, color = true_label)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c(micro = "#6A3D9A", macro = "#E31A1C", `NA` = "grey70")) +
  labs(title = "PC1–PC2 colored by true labels", color = "") +
  theme_classic()
print(p2)

## ===========================================
## 4) GMM clustering (mclust)
## ===========================================
gmm_input <- cbind(as.matrix(pcs[, 1:20, drop = FALSE]),
                   size = pcs$size_log)

gmm <- Mclust(gmm_input, G = 2)
print(summary(gmm))

pcs$gmm_class <- factor(gmm$classification)
pcs$gmm_post1 <- gmm$z[, 1]

## GMM visualization
p3 <- ggplot(pcs, aes(size_log, gmm_post1, color = true_label)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c(micro = "#6A3D9A", macro = "#E31A1C", `NA` = "grey70")) +
  labs(title = "GMM: Posterior probability vs log2(size)",
       x = "log2(group size + 1)",
       y = "Posterior(Class 1)") +
  theme_classic()
print(p3)

## Optional export
# write.csv(pcs, "visium_micro_macro_clustering_results.csv", row.names = FALSE)
