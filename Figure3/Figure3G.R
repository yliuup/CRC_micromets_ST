
library(pheatmap)

# Marker genes (ordered)
genes <- c(
  "TRBC2", "CD3E", "TRAC",                     # T cells
  "CD79A", "MS4A1", "IGHM", "IGHA1", "IGHG1",  # B cells
  "MZB1",                                      # Plasma
  "CXCL14", "C3",                              # iCAF
  "CXCL1",                                     # 
  "THBS2", "POSTN", "COL6A3",                  # myCAF
  "TAGLN", "MYH11", "MYL9",                    # SMC
  "PECAM1", "PLVAP", "RGCC",                   # Endothelial
  "CSF3R", "S100A8", "S100A9",                 # Neutrophil
  "SPP1", "CD68", "TREM2", "CD5L", "C1QA", "MARCO",  # Macrophages
  "ALB", "APOE", "APOB",                       # Hepatocytes
  "EPCAM", "CEACAM5", "CD24",                  # TU-CEACAM5
  "SLC2A1", "ENO2", "KRT19",                   # TU-SLC2A1
  "CDK1", "MKI67", "STMN1",                    # Cycling
  "KRT7", "MUC6", "AQP1"                       # EP-KRT7
)

# Cell clusters (in desired order)
cell_clusters <- c(
  "Tcell", "Bcell", "Plasma",
  "iCAF", "myCAF", "SMC", "Endo",
  "Neutrophil", "Mac-CD68", "Mac-Kupffer",
  "Hepatocyte1", "Hepatocyte2",
  "TU-CEACAM5", "TU-SLC2A1", "TU-Cycling",
  "EP-KRT7"
)

# Keep only clusters present in dat
cell_clusters <- cell_clusters[cell_clusters %in% colnames(dat)]

# Heatmap
pheatmap(
  dat[genes, cell_clusters] %>% t(),
  scale         = "column",
  cluster_row   = FALSE,
  cluster_col   = FALSE,
  color         = cols,
  border        = FALSE,
  cellwidth     = 12,
  cellheight    = 10,
  treeheight_row = 0,
  border_color  = "white",
  angle_col     = 45
)
