# 1. Define Jaccard similarity function
jaccard <- function(x, y) {
  inter <- length(intersect(x, y))
  union_len <- length(union(x, y))
  if (union_len == 0) return(0)
  inter / union_len
}

# 2. Convert cluster -> gene list into a list
cluster_list <- split(fig3_visium_hd$gene, fig3_visium_hd$cluster)

# 3. Initialize Jaccard matrix (rows = fig3_visium, columns = clusters)
jac_mat <- matrix(
  NA,
  nrow = length(fig3_visium),
  ncol = length(cluster_list),
  dimnames = list(
    names(fig3_visium),
    names(cluster_list)
  )
)

# 4. Fill Jaccard matrix
for (i in seq_along(fig3_visium)) {
  for (j in seq_along(cluster_list)) {
    jac_mat[i, j] <- jaccard(fig3_visium[[i]], cluster_list[[j]])
  }
}

# 5. Color palette
cols <- colorRampPalette(brewer.pal(12, "PuOr"))(50) %>% rev()

# 6. Heatmap visualization
library(pheatmap)

pheatmap(
  jac_mat[
    c(
      "liver_iCAF", "liver_myCAF", "liver_SmoothMuscle", "liver_Endo",
      "liver_Hepatocyte", "liver_Mac_CLEC4G", "liver_Mac_CD68",
      "liver_Neutrophil", "liver_Inflammatory", "liver_T_B",
      "liver_EP_KRT7", "liver_Tum_CEACAM5", "liver_Tum_SLC2A1", "liver_Cycling"
    ),
    c(
      "iCAF", "myCAF", "SMC", "Endo",
      "Hepatocyte1", "Hepatocyte2",
      "Mac_Kupffer", "Mac_CD68", "Neutrophil",
      "Tcell", "Bcell", "Plasma",
      "EP_KRT7", "TU_CEACAM5", "TU_SLC2A1", "TU_Cycling"
    )
  ],
  cluster_row   = FALSE,
  cluster_col   = FALSE,
  cellwidth     = 12,
  cellheight    = 12,
  scale         = "row",
  color         = cols,
  border        = FALSE,
  display_numbers = FALSE,
  number_format = "%.2f"
)
