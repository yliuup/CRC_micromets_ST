############################################################
## CellChat spatial pipeline (Visium HD)
## - English comments only, de-duplicated imports
## - Anonymized paths with placeholders
## - Distance-constrained + contact-dependent inference
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(future)
})

## -------------------------------
## 0) User config (EDIT)
## -------------------------------
data_dir   <- "<DATA_DIR>"        # e.g., ".../3.Seurat/"
rdata_file <- file.path(data_dir, "dat_split.RData")
out_dir    <- "<OUTPUT_DIR>"      # e.g., ".../cellchat_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Parallel settings
plan("multicore", workers = 5)
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB; adjust as needed

## Spatial parameters (tune for your platform)
## - distance.use=TRUE enables distance constraint
## - interaction.range: max diffusion distance for diffusible ligands (in same unit as coordinates)
## - contact.dependent + contact.range: threshold for contact-based signaling
interaction_range <- 200
contact_range     <- 100
k_min_groups      <- 3   # min groups for truncatedMean

## -------------------------------
## 1) Load Seurat object
## -------------------------------
stopifnot(file.exists(rdata_file))
load(rdata_file)  # expects an object 'dat_split'
stopifnot(exists("dat_split"))
dat <- dat_split[[2]]              # choose your subset; adjust if needed
stopifnot(inherits(dat, "Seurat"))

## Normalized expression (RNA assay)
data.input <- Seurat::GetAssayData(dat, slot = "data", assay = "RNA")

## -------------------------------
## 2) Build meta & spatial coordinates
## -------------------------------
## Meta must provide: 'anno' (group labels), 'sample' (replicate labels)
## Here we reuse the active identity as anno and ensure a 'sample' column exists.
if (!"sample" %in% colnames(dat@meta.data)) {
  stop("meta.data must contain a 'sample' column for spatial aggregation.")
}
dat$anno <- as.character(Idents(dat))  # or your curated cluster labels
meta <- dat@meta.data[, c("anno", "sample")]
meta$sample <- factor(meta$sample)
meta$anno   <- paste0("C_", meta$anno) # optional prefix to avoid DB name collision

## Spatial coordinates must be named 'imagerow' and 'imagecol'
## If you already have imagerow/imagecol, just select them; otherwise map row/col.
coord_cols <- c("row", "col")
if (all(coord_cols %in% colnames(dat@meta.data))) {
  spatial.locs <- dat@meta.data[, coord_cols, drop = FALSE]
  colnames(spatial.locs) <- c("imagerow", "imagecol")
} else if (all(c("imagerow", "imagecol") %in% colnames(dat@meta.data))) {
  spatial.locs <- dat@meta.data[, c("imagerow", "imagecol"), drop = FALSE]
} else {
  stop("No spatial coordinates found. Expect 'row/col' or 'imagerow/imagecol' in meta.data.")
}

## Spatial factors:
## - ratio: unit conversion factor from coordinate unit to microns (set to 1 if already ~um)
## - tol  : tolerance for nearest-neighbor graph construction
spatial.factors <- data.frame(ratio = 1, tol = 5)  # DO NOT use the 'factor' type here

## -------------------------------
## 3) Create CellChat object
## -------------------------------
cellchat <- createCellChat(
  object           = data.input,
  meta             = meta,
  group.by         = "anno",
  datatype         = "spatial",
  coordinates      = spatial.locs,
  spatial.factors  = spatial.factors
)

## Database
CellChatDB      <- CellChatDB.human   # or CellChatDB.mouse
cellchat@DB     <- CellChatDB         # use default v2 DB (includes non-protein signaling)

## -------------------------------
## 4) Subset signaling genes & prefilter interactions
## -------------------------------
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = FALSE)

## -------------------------------
## 5) Compute communication probability (spatial)
## -------------------------------
ptm <- Sys.time()
cellchat <- computeCommunProb(
  cellchat,
  type              = "truncatedMean",
  trim              = 0.01,
  k.min             = k_min_groups,
  distance.use      = TRUE,
  interaction.range = interaction_range,
  scale.distance    = 0.5,
  contact.dependent = TRUE,
  contact.range     = contact_range
)
elapsed <- difftime(Sys.time(), ptm, units = "secs")
message(sprintf("[OK] computeCommunProb finished in %.2f sec", as.numeric(elapsed)))

## -------------------------------
## 6) Filter, pathway-level aggregation & centrality
## -------------------------------
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)                       # aggregate net (cell-cell, pathway)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

## -------------------------------
## 7) Save outputs
## -------------------------------
saveRDS(cellchat, file = file.path(out_dir, "cellchat_spatial.rds"))

## Optionally export key tables
df.net      <- subsetCommunication(cellchat)             # LR pairs (all/pathway/aggregate)
df.net.path <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.net,      file.path(out_dir, "communication_pairs.csv"), row.names = FALSE)
write.csv(df.net.path, file.path(out_dir, "communication_pairs_pathway.csv"), row.names = FALSE)

## -------------------------------
## 8) Quick QA plots (optional)
## -------------------------------
# pdf(file.path(out_dir, "net_overview.pdf"), width = 9, height = 7)
# netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)),
#                  weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)),
#                  weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength")
# dev.off()

message("[DONE] CellChat spatial pipeline completed.")
