
##tree

library(ape)
library(phangorn)
library(ggtree)
library(ggsci)
library(RColorBrewer)
library(scales)

## avg.cnv.df: samples x features matrix/data.frame
## qq: distance method (e.g. "euclidean", "manhattan", etc.)
## i:  index or sample group id used in path/output

# 1. Distance and NJ tree
data_dist <- dist(avg.cnv.df, method = qq)
tree_nj   <- NJ(data_dist)

# 2. Root the tree at "Normal"
rooted_tree <- root(tree_nj, outgroup = "Normal", resolve.root = TRUE)

# 3. Define group information (here: each tip as its own group)
tip_labels <- rooted_tree$tip.label
cluster    <- tip_labels
names(cluster) <- tip_labels

groupInfo <- lapply(
  levels(as.factor(cluster)),
  function(cc, cluster_vec) {
    names(cluster_vec)[cluster_vec == cc]
  },
  cluster_vec = cluster
)
names(groupInfo) <- names(table(cluster))

# 4. Colors
cols <- c(
  pal_npg()(10),
  pal_igv()(9),
  pal_uchicago("light")(9),
  pal_futurama()(12),
  pal_aaas()(10)
)[-8]

names(cols) <- paste0("Clone_", seq_len(20))
my_cols_used_v1 <- c(cols, "Normal" = "lightgrey")

# 5. Group tree and plot
tree_grouped <- groupOTU(rooted_tree, groupInfo)

gp <- ggtree(tree_grouped, aes(color = group), size = 2) +
  theme_tree() +
  scale_color_manual(values = my_cols_used_v1) +
  geom_tippoint(size = 15) +
  geom_tiplab(hjust = 0.55, vjust = 0.45, color = "black", size = 3) +
  theme(legend.position = "none")

###heatmap

library(ComplexHeatmap)
library(ggplot2)
library(grid)  # for gpar, unit

## ---------------- Chromosome factor and colors ----------------

# Set chromosome order based on human genome convention
chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY")

# Keep only chromosomes present in anno.col$V2
n_chr <- length(unique(anno.col$V2))
chr_levels <- chr_levels[seq_len(n_chr)]

# Factor chromosome annotation with defined levels
anno.col$chr <- factor(anno.col$V2, levels = chr_levels)

# Alternating colors for chromosomes
chr_colors <- rep(c("lightgrey", "darkgrey"), length.out = n_chr)

mat_colors <- list(chr = chr_colors)
names(mat_colors$chr) <- levels(anno.col$chr)

## ---------------- CNV color scales ----------------

cols_gain <- c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "white")
colfunc_gain <- colorRampPalette(cols_gain)

cols_loss <- c("white", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061")
colfunc_loss <- colorRampPalette(cols_loss)

# Build CNV color vector based on avg.cnv.df range
cnv.cols <- c(
  colfunc_gain(20),
  "white",
  colfunc_loss(
    round(
      20 * (1 - min(avg.cnv.df)) / (max(avg.cnv.df) - 1),
      0
    )
  )
)

## ---------------- Column annotation (chromosome) ----------------

dat <- data.frame(chr = anno.col[,"chr"])

top_ann <- HeatmapAnnotation(
  df = dat,
  which = "column",
  annotation_name_side = "left",
  col = mat_colors,
  show_legend = FALSE
)

## ---------------- Prepare matrix and row order ----------------

consensus_df <- avg.cnv.df
clusters     <- meta.data
ploidy_VAL   <- 1


## ---------------- Heatmap ----------------

ht_list <- Heatmap(
  avg.cnv.df,
  name = "inferCNV",
  col = rev(cnv.cols),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  top_annotation = top_ann,
  row_split = factor(rownames(avg.cnv.df), levels = rownames(avg.cnv.df)),
  row_title = NULL,
  column_split = anno.col$chr,
  column_title = NULL,
  row_gap = unit(0, "mm"),
  column_gap = unit(0, "mm"),
  border = TRUE,
  heatmap_legend_param = list(
    title_gp   = gpar(fontsize = 12, fontface = "bold"),
    labels_gp  = gpar(fontsize = 12),
    legend_height = unit(5, "cm"),
    title_position = "topleft"
  )
)


###spatial plot


library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(stringr)
library(ggpubr)

## ---------------------------- configuration ----------------------------

sample_groups <- list(
  c("M-ST-13", "M-ST-05", "M-ST-14"),
  c("M-ST-07", "M-ST-08"),
  c("M-ST-06"),
  c("M-ST-09", "M-ST-10"),
  c("IX_56688", "IX_56689", "IX_56690",
    "M-ST-17", "M-ST-19", "M-ST-16", "M-ST-20"),
  c("M-ST-15", "M-ST-11", "M-ST-12", "M-ST-18"),
  c("M-ST-23", "M-ST-24"),
  c("M-ST-27", "M-ST-38", "M-ST-28"),
  c("M-ST-21", "M-ST-22", "M-ST-31", "M-ST-32", "M-ST-37"),
  c("M-ST-35", "M-ST-33", "M-ST-34", "M-ST-36"),
  c("M-ST-25", "M-ST-26", "M-ST-29", "M-ST-30")
)

patients <- c(
  "PAT1", "PAT2", "PAT3", "PAT4", "PAT5",
  "PAT6", "PAT7", "PAT8", "PAT9", "PAT10", "PAT11"
)


setwd(base_dir)

## sce_list, crc_merge_filt_pathway, dat_cyto_merge
load("crc_merge_filt_pathway.RData")

## ---------------------------- helper functions ----------------------------

get_range <- function(x) {
  c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
}

make_feature_plot <- function(
    sce,
    values,
    low,
    mid = NULL,
    high,
    limits,
    use_mid = TRUE
) {
  names(values) <- colnames(sce)
  
  p <- featurePlot(sce, values, color = NA)
  
  if (use_mid) {
    midpoint <- mean(limits, na.rm = TRUE)
    p <- p +
      scale_fill_gradient2(
        low = low,
        mid = mid,
        high = high,
        na.value = "lightgrey",
        limits = limits,
        midpoint = midpoint
      )
  } else {
    p <- p +
      scale_fill_gradient(
        low = low,
        high = high,
        na.value = "lightgrey",
        limits = limits
      )
  }
  
  p + theme(legend.position = "none")
}

plot_quiescence <- function(sce, sample_id, limits, feature_name = "Quiescence") {
  vals <- as.numeric(sce_list[[sample_id]][[feature_name]])
  make_feature_plot(
    sce = sce,
    values = vals,
    low  = "#80CDC1",
    mid  = "#F5F5F5",
    high = "#BF812D",
    limits = limits,
    use_mid = TRUE
  )
}

plot_differentiation <- function(sce, sample_id, limits, feature_name = "Differentiation") {
  vals <- as.numeric(sce_list[[sample_id]][[feature_name]])
  make_feature_plot(
    sce = sce,
    values = vals,
    low  = "#7FBC41",
    mid  = "#F7F7F7",
    high = "#DE77AE",
    limits = limits,
    use_mid = TRUE
  )
}

plot_cellcycle <- function(sce, sample_id, limits, feature_name = "CellCycle") {
  vals <- as.numeric(sce_list[[sample_id]][[feature_name]])
  make_feature_plot(
    sce = sce,
    values = vals,
    low  = "#8073AC",
    mid  = "#F7F7F7",
    high = "#E08214",
    limits = limits,
    use_mid = TRUE
  )
}

plot_cytotrace <- function(sce, sample_id, limits, feature_name = "cytotrace") {
  vals <- as.numeric(sce_list[[sample_id]][[feature_name]])
  make_feature_plot(
    sce = sce,
    values = vals,
    low  = "#2166AC",
    mid  = "#F7F7F7",
    high = "#B2182B",
    limits = limits,
    use_mid = TRUE
  )
}

plot_mki67 <- function(sce, sample_id, limits, feature_name = "MKI67") {
  vals <- as.numeric(sce_list[[sample_id]][[feature_name]])
  make_feature_plot(
    sce = sce,
    values = vals,
    low  = "#FFF5F0",
    mid  = NULL,
    high = "red",
    limits = limits,
    use_mid = FALSE
  )
}

## ---------------------------- CytoTRACE sample column ----------------------------

dat_cyto_merge$sample <- str_split(dat_cyto_merge$bar2, pattern = "~", simplify = TRUE)[, 1]

## ---------------------------- main loop: per patient ----------------------------

for (i in seq_along(sample_groups)) {
  
  group_samples <- sample_groups[[i]]
  cat("Processing", patients[i], "with samples:", paste(group_samples, collapse = ", "), "\n")
  
  # ranges across all samples of this patient group (for consistent color scale)
  range_quiescence <- crc_merge_filt_pathway$Quiescence1[
    crc_merge_filt_pathway$sample %in% group_samples
  ] |> get_range()
  
  range_diff <- crc_merge_filt_pathway$Differentiation1[
    crc_merge_filt_pathway$sample %in% group_samples
  ] |> get_range()
  
  range_cc <- crc_merge_filt_pathway$CellCycle1[
    crc_merge_filt_pathway$sample %in% group_samples
  ] |> get_range()
  
  range_cyto <- dat_cyto_merge$val[
    dat_cyto_merge$sample %in% group_samples
  ] |> get_range()
  
  # MKI67 range: computed from sce_list for all samples in this group
  mki_vals_all <- c()
  for (s_id in group_samples) {
    mki_vals_all <- c(mki_vals_all, as.numeric(sce_list[[s_id]][["MKI67"]]))
  }
  range_mki <- get_range(mki_vals_all)
  
  # collect plots for this patient
  plots <- list()
  
  for (sample_id in group_samples) {
    obj <- sce_list[[sample_id]]
    
    p_cyto <- plot_cytotrace(
      sce = obj,
      sample_id = sample_id,
      limits = range_cyto
    )
    plots <- c(plots, list(p_cyto))
    
    p_quies <- plot_quiescence(
      sce = obj,
      sample_id = sample_id,
      limits = range_quiescence
    )
    plots <- c(plots, list(p_quies))
    
    p_diff <- plot_differentiation(
      sce = obj,
      sample_id = sample_id,
      limits = range_diff
    )
    plots <- c(plots, list(p_diff))
    
    p_cc <- plot_cellcycle(
      sce = obj,
      sample_id = sample_id,
      limits = range_cc
    )
    plots <- c(plots, list(p_cc))
    
    p_mki <- plot_mki67(
      sce = obj,
      sample_id = sample_id,
      limits = range_mki
    )
    plots <- c(plots, list(p_mki))
  }
  
  n_samples <- length(group_samples)
  combined.gp <- do.call(
    ggarrange,
    c(plots, ncol = 5, nrow = n_samples)
  )
  
  pdf_file <- paste0("Spatial_", patients[i], ".pdf")
  pdf(pdf_file, height = 3 * n_samples, width = 15)
  print(combined.gp)
  dev.off()
}

