## ===================== Libraries =====================

library(plyr)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)

## ===================== Helpers =====================

paired_boxplot <- function(df, value_col, title, color_hex) {
  ggplot(df, aes(x = anno, y = .data[[value_col]], fill = anno, color = anno)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_point(aes(group = sample), width = 0.2, size = 2, alpha = 0.7) +
    geom_line(aes(group = sample), color = "gray", linetype = "dashed", alpha = 0.5) +
    stat_compare_means(
      method = "t.test",
      paired = TRUE,
      label = "p.format"
    ) +
    labs(title = title, x = "", y = "") +
    theme_classic2() +
    scale_color_manual(values = c("lightgrey", color_hex)) +
    scale_fill_manual(values = alpha(c("lightgrey", color_hex), 0.5)) +
    theme(legend.position = "none")
}

## ===================== 1. Pathway scores (Quiescence, Differentiation, CellCycle) =====================



# Map annotation and patient to pathway object
crc_merge_filt_pathway$anno <- mapvalues(
  x    = rownames(crc_merge_filt_pathway),
  from = rownames(crc_merge_filt_meta4),
  to   = crc_merge_filt_meta4$new2
)
crc_merge_filt_pathway$patient <- mapvalues(
  x    = rownames(crc_merge_filt_pathway),
  from = rownames(crc_merge_filt_meta4),
  to   = crc_merge_filt_meta4$patient
)

# Keep macro/micro metastasis only
crc_merge_filt_pathway_filt <- crc_merge_filt_pathway[
  grep("macrometastasis|micrometastasis", crc_merge_filt_pathway$anno),
]

# Keep selected columns: (1, 8, 15:30) as in original code
crc_merge_filt_pathway_filt <- crc_merge_filt_pathway_filt[, c(1, 8, 15:30)]

# Group label: patient~anno
crc_merge_filt_pathway_filt$group <- paste0(
  crc_merge_filt_pathway_filt$patient, "~", crc_merge_filt_pathway_filt$anno
)

# Summarize mean per group
df_summary <- crc_merge_filt_pathway_filt[, c(3:17, 19)] %>%
  group_by(group) %>%
  summarize(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

# Split group into sample & anno
df_summary$sample <- str_split(df_summary$group, pattern = "~", simplify = TRUE)[, 1]
df_summary$anno   <- str_split(df_summary$group, pattern = "~", simplify = TRUE)[, 2]

# Melt to long format
df_summary_melt <- melt(df_summary)

# Keep liver tumor only
df_liver <- df_summary_melt %>%
  filter(grepl("Liver", anno), grepl("tumor", anno))

df_liver$anno <- factor(
  df_liver$anno,
  levels = c("Liver macrometastasis tumor", "Liver micrometastasis tumor")
)

# Split for each pathway
df_quies <- df_liver[df_liver$variable == "Quiescence1", ]
df_diff  <- df_liver[df_liver$variable == "Differentiation1", ]
df_cc    <- df_liver[df_liver$variable == "CellCycle1", ]

# Plots
p1 <- paired_boxplot(df_quies, "value", "Quiescence",      "#2171B5")
p2 <- paired_boxplot(df_diff,  "value", "Differentiation", "#238B45")
p3 <- paired_boxplot(df_cc,    "value", "CellCycle",       "#D94801")

## ===================== 2. CytoTRACE score =====================

dat_cyto_merge$anno <- mapvalues(
  x    = dat_cyto_merge$bar,
  from = rownames(crc_merge_filt_meta4),
  to   = crc_merge_filt_meta4$new2
)
dat_cyto_merge$sample <- mapvalues(
  x    = dat_cyto_merge$bar,
  from = rownames(crc_merge_filt_meta4),
  to   = crc_merge_filt_meta4$patient
)
dat_cyto_merge$group <- paste0(dat_cyto_merge$sample, "~", dat_cyto_merge$anno)

df_cyto_summary <- dat_cyto_merge[, c("group", "val")] %>%
  group_by(group) %>%
  summarize(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

df_cyto_summary$sample <- str_split(df_cyto_summary$group, pattern = "~", simplify = TRUE)[, 1]
df_cyto_summary$anno   <- str_split(df_cyto_summary$group, pattern = "~", simplify = TRUE)[, 2]

df_cyto_liver <- df_cyto_summary[grep("Liver", df_cyto_summary$anno), ]
df_cyto_liver$anno <- factor(
  df_cyto_liver$anno,
  levels = c("Liver macrometastasis tumor", "Liver micrometastasis tumor")
)

p_cyto <- paired_boxplot(df_cyto_liver, "val", "Cytotrace Score", "#AE017E")

## ===================== 3. CNV score =====================


score1_dat_filt$group <- paste0(score1_dat_filt$patient, "~", score1_dat_filt$anno)

df_cnv_summary <- score1_dat_filt[, c("group", "cnv.score")] %>%
  group_by(group) %>%
  summarize(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

df_cnv_summary$sample <- str_split(df_cnv_summary$group, pattern = "~", simplify = TRUE)[, 1]
df_cnv_summary$anno   <- str_split(df_cnv_summary$group, pattern = "~", simplify = TRUE)[, 2]

df_cnv_liver <- df_cnv_summary[grep("Liver", df_cnv_summary$anno), ]
df_cnv_liver$anno <- factor(
  df_cnv_liver$anno,
  levels = c("Liver macrometastasis tumor", "Liver micrometastasis tumor")
)

p_cnv <- paired_boxplot(df_cnv_liver, "cnv.score", "CNV Score", "#E74C3C")

## ===================== 4. Combine plots =====================

plot_list <- list(p_cnv, p_cyto, p1, p2, p3)

combined.gp <- do.call(ggarrange, c(plot_list, nrow = 1))
print(combined.gp)
