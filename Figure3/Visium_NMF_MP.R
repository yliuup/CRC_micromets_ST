############################################################
## Spatial NMF → Metaprograms (MPs)
############################################################

set.seed(1234)

############################################################
## 0) Configuration (EDIT THESE)
############################################################
data_dir   <- "<DATA_DIRECTORY_PATH>"         # directory containing *_nmf.RData and /processed matrices
script_dir <- "<CUSTOM_FUNCTIONS_DIRECTORY>"  # directory containing helper scripts (if available)
out_dir    <- "<OUTPUT_DIRECTORY_PATH>"       # e.g., "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Thresholds / parameters
TOP_N_HVG              <- 50   # top genes per NMF program (e.g., 50)
MIN_INTERSECT_INITIAL  <- 15   # min overlap to seed a cluster
MIN_INTERSECT_CLUSTER  <- 15   # min overlap to add to cluster
MIN_GROUP_SIZE         <- 4    # minimal group size for a founder cluster
EXPR_CENTER            <- TRUE # scalop::sigScores parameter
EXPR_NBIN              <- 15   # scalop::sigScores parameter

############################################################
## 1) Libraries
############################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(reshape2)
  library(ggplot2)
  library(scales)
  library(data.table)
})

## Optional: scoring (requires scalop)
have_scalop <- requireNamespace("scalop", quietly = TRUE)

############################################################
## 2) Optional helper scripts 
############################################################
robust_filter_available <- FALSE
try({
  src1 <- file.path(script_dir, "custom_magma.R")
  src2 <- file.path(script_dir, "dbscan_get_sig.R")
  src3 <- file.path(script_dir, "robust_nmf_programs.R")
  src4 <- file.path(script_dir, "nmf_cell_class.R")
  if (file.exists(src3)) {
    source(src3)
    robust_filter_available <- TRUE
  }
  if (file.exists(src1)) source(src1)
  if (file.exists(src2)) source(src2)
  if (file.exists(src4)) source(src4)
}, silent = TRUE)

## Fallback robust filter: keep all columns
if (!robust_filter_available) {
  message("[INFO] 'robust_nmf_programs.R' not found. Using identity filter (keep all NMF programs).")
  robust_nmf_programs <- function(x, intra_min=35, intra_max=10, inter_filter=TRUE, inter_min=10) {
    unique(unlist(lapply(x, colnames)))
  }
}

############################################################
## 3) Collect NMF runs (expects *_nmf.RData with objects w,h)
############################################################
nmf_files <- list.files(path = data_dir, pattern = "_nmf\\.RData$", full.names = TRUE)
stopifnot(length(nmf_files) > 0)

nmf_genes <- list()   # list of W (genes x programs)
nmf_cells <- list()   # list of H (programs x cells/spots)

for (f in nmf_files) {
  env <- new.env()
  load(f, envir = env)  # expects env$w and env$h
  if (!all(c("w","h") %in% ls(env))) {
    stop("File ", f, " does not contain objects 'w' and 'h'.")
  }
  key <- gsub("_nmf\\.RData$", "", basename(f))
  nmf_genes[[key]] <- env$w
  nmf_cells[[key]] <- env$h
}

save(nmf_genes, nmf_cells, file = file.path(out_dir, "nmf_collected.RData"))
message("[OK] Collected NMF runs: ", length(nmf_genes))

############################################################
## 4) Extract top-N (e.g., 50) genes for each NMF program
############################################################
nmf_topN <- lapply(nmf_genes, function(mat) {
  apply(mat, 2, function(v) names(sort(v, decreasing = TRUE))[seq_len(min(TOP_N_HVG, length(v)))])
})
## Filter robust metaprograms (if helper is available)
nmf_filter <- robust_nmf_programs(
  nmf_topN, intra_min = 35, intra_max = 10, inter_filter = TRUE, inter_min = 10
)
nmf_topN <- lapply(nmf_topN, function(x) x[, colnames(x) %in% nmf_filter, drop = FALSE])

## Merge all programs to one matrix (genes per column/program)
nmf_topN_merged <- do.call(cbind, nmf_topN)
stopifnot(ncol(nmf_topN_merged) > 0)
message("[OK] NMF programs after robust filter: ", ncol(nmf_topN_merged))

save(nmf_topN_merged, file = file.path(out_dir, "nmf_topN_merged.RData"))

############################################################
## 5) Build intersection matrix (overlap counts among programs)
############################################################
message("[..] Computing intersection matrix ...")
intersect_mat <- apply(nmf_topN_merged, 2, function(g1)
  apply(nmf_topN_merged, 2, function(g2) length(intersect(g1, g2)))
)
save(intersect_mat, file = file.path(out_dir, "nmf_intersect_mat.RData"))
message("[OK] Intersection matrix dimension: ", paste(dim(intersect_mat), collapse = " x "))

############################################################
## 6) Iterative clustering to identify MP clusters
############################################################
Cluster_list <- list()
MP_list <- list()
k <- 1

M <- intersect_mat
NMF <- nmf_topN_merged

while (TRUE) {
  if (ncol(M) == 0 || nrow(M) == 0) break
  ## number of neighbors above initial cutoff
  neighbor_counts <- sort(apply(M, 2, function(x) sum(x >= MIN_INTERSECT_INITIAL) - 1), decreasing = TRUE)
  if (length(neighbor_counts) == 0 || neighbor_counts[1] <= MIN_GROUP_SIZE) break
  
  seed <- names(neighbor_counts)[1]
  curr <- seed
  MP_genes <- NMF[, seed, drop = TRUE]
  
  ## remove seed from pool for candidate extension
  NMF <- NMF[, colnames(NMF) != seed, drop = FALSE]
  
  ## compute intersection with current MP genes
  if (ncol(NMF) > 0) {
    add_scores <- sort(apply(NMF, 2, function(x) length(intersect(MP_genes, x))), decreasing = TRUE)
  } else {
    add_scores <- c()
  }
  
  while (length(add_scores) > 0 && add_scores[1] >= MIN_INTERSECT_CLUSTER) {
    nxt <- names(add_scores)[1]
    curr <- c(curr, nxt)
    ## update pool
    NMF <- NMF[, colnames(NMF) != nxt, drop = FALSE]
    if (ncol(NMF) > 0) {
      add_scores <- sort(apply(NMF, 2, function(x) length(intersect(MP_genes, x))), decreasing = TRUE)
    } else {
      add_scores <- c()
    }
  }
  
  Cluster_list[[paste0("Cluster_", k)]] <- curr
  MP_list[[paste0("MP_", k)]] <- MP_genes
  ## remove cluster programs from M
  keep_rows <- !(rownames(M) %in% curr)
  keep_cols <- !(colnames(M) %in% curr)
  M <- M[keep_rows, keep_cols, drop = FALSE]
  k <- k + 1
}

save(Cluster_list, MP_list, file = file.path(out_dir, "mp_clusters.RData"))
message("[OK] Identified MP clusters: ", length(MP_list))

############################################################
## 7) (Optional) Heatmap of program similarity
############################################################
try({
  m_melt <- reshape2::melt(intersect_mat)
  p <- ggplot(m_melt, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black") +
    coord_equal() +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = "NMF program overlap (gene intersections)")
  ggsave(file.path(out_dir, "nmf_similarity_heatmap.pdf"), p, width = 6, height = 6)
  message("[OK] Saved heatmap: nmf_similarity_heatmap.pdf")
}, silent = TRUE)

############################################################
## 8) Score MPs on spatial samples (if matrices available)
##    Expect per-sample matrices under:
##    file.path(data_dir, "processed/") as RDS files
##    rows = spots/bins; cols = genes; normalized/logCPM recommended
############################################################
processed_dir <- file.path(data_dir, "processed")
if (dir.exists(processed_dir) && have_scalop) {
  files <- list.files(processed_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) > 0) {
    message("[..] Scoring MPs on ", length(files), " processed matrices ...")
    ## Use MP_list (character vectors of genes) as signatures
    meta_program <- MP_list
    
    score_list <- lapply(files, function(f) {
      m <- readRDS(f)  # rows = spots/bins, cols = genes
      ## Ensure matrix orientation: scalop::sigScores expects rows=samples/spots, cols=genes
      if (!is.matrix(m)) m <- as.matrix(m)
      ## intersect by gene names
      ## sigScores will internally handle missing genes; we pass the list of vectors
      sigs <- scalop::sigScores(m, meta_program, expr.center = EXPR_CENTER, expr.nbin = EXPR_NBIN)
      df <- data.frame(spot = rownames(sigs), sigs, check.names = FALSE)
      df$sample <- basename(f)
      df
    })
    
    scores <- do.call(rbind, score_list)
    saveRDS(scores, file.path(out_dir, "mp_scores.rds"))
    message("[OK] Saved MP scores: mp_scores.rds")
    
    ## Assign max-scoring MP per spot
    score_max <- lapply(score_list, function(df){
      rn <- df$spot
      mat <- as.matrix(df[, setdiff(colnames(df), c("spot","sample")), drop = FALSE])
      max_idx <- max.col(mat, ties.method = "first")
      max_lab <- colnames(mat)[max_idx]
      data.frame(spot = rn, cluster = max_lab, sample = unique(df$sample), check.names = FALSE)
    })
    nmf_clu_dat <- do.call(rbind, score_max)
    save(nmf_clu_dat, file = file.path(out_dir, "nmf_clu_dat.rdata"))
    message("[OK] Saved MP max-cluster assignments: nmf_clu_dat.rdata")
  } else {
    message("[INFO] No processed matrices found under: ", processed_dir)
  }
} else {
  if (!have_scalop) {
    message("[INFO] Package 'scalop' not available; skipping MP scoring.")
  } else {
    message("[INFO] Directory not found for processed matrices: ", processed_dir)
  }
}

############################################################
## 9) Session info (for reproducibility)
############################################################
if (requireNamespace("sessioninfo", quietly = TRUE)) {
  writeLines(capture.output(sessioninfo::session_info()),
             con = file.path(out_dir, "session_info.txt"))
} else {
  writeLines(capture.output(sessionInfo()),
             con = file.path(out_dir, "session_info.txt"))
}
message("[DONE] Pipeline completed.")
