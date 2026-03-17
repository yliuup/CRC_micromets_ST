############################################################
## Hex-like neighborhood counting (6 neighbors) on Visium grid
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(plyr)     # for mapvalues
  library(reshape2)
})

## -------------------------------
## 0) Parameters / paths (EDIT)
## -------------------------------
data_dir   <- "<DATA_DIR>"   # e.g., ".../9.nb/data"
result_dir <- "<RESULT_DIR>" # e.g., ".../9.nb/res"
meta_file  <- file.path(data_dir, "nmf_clu_dat_merge.RData")  # contains `nmf_clu_dat_merge`

## Spot location file pattern: "<sample_id>_loc_dat.csv"
## The sample id will be passed via commandArgs(trailingOnly=TRUE)[1]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
sample_id <- args[1]

loc_file  <- file.path(data_dir, paste0(sample_id, "_loc_dat.csv"))
out_file  <- file.path(result_dir, paste0(sample_id, "nb_raw.RData"))

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

## -------------------------------
## 1) Helpers
## -------------------------------

# Return an integer vector of neighbor indices for a single (row, col)
# using the 6-neighbor offset pattern observed in the provided script.
get_neighbor_indices <- function(position_df, r, c) {
  idx <- c(
    which(position_df$row == r     & position_df$col == (c - 2)),
    which(position_df$row == r     & position_df$col == (c + 2)),
    which(position_df$row == r - 1 & position_df$col == (c - 1)),
    which(position_df$row == r - 1 & position_df$col == (c + 1)),
    which(position_df$row == r + 1 & position_df$col == (c - 1)),
    which(position_df$row == r + 1 & position_df$col == (c + 1))
  )
  unique(idx)
}

# Build a neighbors table (rows = focal spots by bar; columns = neighbor slots)
# Each cell holds the neighbor's cluster label (or NA if missing)
build_neighbors_table <- function(position_df) {
  pos_list <- vector("list", nrow(position_df))
  for (j in seq_len(nrow(position_df))) {
    rj <- position_df$row[j]
    cj <- position_df$col[j]
    n_idx <- get_neighbor_indices(position_df, rj, cj)
    pos_list[[j]] <- position_df$anno[n_idx]
  }
  max_len <- max(lengths(pos_list))
  pos_padded <- lapply(pos_list, function(x) { length(x) <- max_len; x })
  neigh_tab <- do.call(rbind, pos_padded)
  rownames(neigh_tab) <- position_df$bar
  as.data.frame(neigh_tab, stringsAsFactors = FALSE)
}

# Count neighbor cluster frequencies per row (spot), then normalize
# by the number of non-NA neighbors (i.e., up to 6).
count_and_normalize <- function(neighbors_table, valid_levels) {
  # count table per row with fixed levels
  count_mat <- t(apply(neighbors_table, 1, function(x) {
    as.integer(table(factor(x, levels = valid_levels)))
  }))
  colnames(count_mat) <- valid_levels
  na_per_row <- apply(neighbors_table, 1, function(x) sum(is.na(x)))
  denom <- (ncol(neighbors_table) - na_per_row)
  denom[denom == 0] <- NA_real_
  norm_mat <- sweep(count_mat, 1, denom, "/")
  norm_mat[is.na(norm_mat)] <- 0
  as.data.frame(norm_mat, stringsAsFactors = FALSE)
}

## -------------------------------
## 2) Load metadata and locations
## -------------------------------
if (!file.exists(meta_file)) stop("Meta file not found: ", meta_file)
if (!file.exists(loc_file))  stop("Location file not found: ", loc_file)

load(meta_file)  # expects `nmf_clu_dat_merge`
if (!exists("nmf_clu_dat_merge")) stop("`nmf_clu_dat_merge` not in meta file")

anno_df <- nmf_clu_dat_merge[nmf_clu_dat_merge$sample == sample_id, ]
if (nrow(anno_df) == 0) stop("No annotations for sample: ", sample_id)

position <- read.csv(loc_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
req_cols <- c("row", "col", "bar")
if (!all(req_cols %in% colnames(position))) {
  stop("Position file must contain columns: ", paste(req_cols, collapse = ", "))
}

## Map annotation from anno_df to position by bar/spot_names
position$anno <- plyr::mapvalues(
  x    = position$bar,
  from = anno_df$spot_names,
  to   = as.character(anno_df$cluster),
  warn_missing = FALSE
)

## Optional: drop spots without annotation (safer than comparing to bar names)
position <- position[!is.na(position$anno), , drop = FALSE]

## -------------------------------
## 3) Neighborhood table
## -------------------------------
neighbors_table <- build_neighbors_table(position)

## -------------------------------
## 4) Count + normalize
## -------------------------------
valid_levels <- unique(position$anno)
cluster_counts_df_nor <- count_and_normalize(neighbors_table, valid_levels)

## Attach identifiers
cluster_counts_df_nor$sample <- sample_id
cluster_counts_df_nor$bar    <- rownames(cluster_counts_df_nor)

## -------------------------------
## 5) Save results
## -------------------------------
save(cluster_counts_df_nor, file = out_file)
message("[OK] Neighbor-normalized cluster frequencies saved: ", out_file)


######

############################################################
## Minimal distance from each spot to each cluster (per sample)
############################################################

suppressPackageStartupMessages({
  library(stringr)
  library(plyr)     # for mapvalues
  library(dplyr)
})

## -------------------------------
## 0) Parameters / paths (EDIT)
## -------------------------------
data_dir   <- "<DATA_DIR>"   # e.g., ".../10.dis/data"
result_dir <- "<RESULT_DIR>" # e.g., ".../10.dis.res"
meta_file  <- file.path(data_dir, "nmf_clu_dat_merge.RData")  # provides `nmf_clu_dat_merge`

## The sample id is passed via command line: Rscript script.R <sample_id>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
sample_id <- args[1]

loc_file  <- file.path(data_dir, paste0(sample_id, "_loc_dat.csv"))
out_file  <- file.path(result_dir, paste0(sample_id, "_min_distances_df.RData"))

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

## -------------------------------
## 1) Load metadata and locations
## -------------------------------
if (!file.exists(meta_file)) stop("Meta file not found: ", meta_file)
if (!file.exists(loc_file))  stop("Location file not found: ", loc_file)

load(meta_file)  # expects `nmf_clu_dat_merge`
if (!exists("nmf_clu_dat_merge")) stop("`nmf_clu_dat_merge` not found in meta file")

anno_df <- nmf_clu_dat_merge[nmf_clu_dat_merge$sample == sample_id, ]
if (nrow(anno_df) == 0) stop("No annotations for sample: ", sample_id)

position <- read.csv(loc_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
req_cols <- c("row", "col", "bar")
if (!all(req_cols %in% colnames(position))) {
  stop("Position file must contain columns: ", paste(req_cols, collapse = ", "))
}

## Map cluster labels from annotation table to position table by bar/spot_names
position$anno <- plyr::mapvalues(
  x    = position$bar,
  from = anno_df$spot_names,
  to   = as.character(anno_df$cluster),
  warn_missing = FALSE
)

## Keep only spots with a valid annotation
position <- position[!is.na(position$anno), , drop = FALSE]
if (nrow(position) == 0) stop("All positions lack annotations after mapping.")

## Ensure numeric coordinates
position$row <- as.numeric(position$row)
position$col <- as.numeric(position$col)

## -------------------------------
## 2) Distance helpers
## -------------------------------
euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

## Compute the minimal distance from each focal spot to each cluster
## Returns a data.frame: rows = spots (bar), cols = clusters (levels of `position$anno`)
calculate_min_distances <- function(position_df) {
  clusters <- unique(position_df$anno)
  n_spots  <- nrow(position_df)
  # Pre-allocate matrix with Inf
  out_mat  <- matrix(Inf, nrow = n_spots, ncol = length(clusters))
  colnames(out_mat) <- clusters
  rownames(out_mat) <- position_df$bar
  
  ## Pre-split by cluster to avoid repeated filtering
  cluster_index <- split(seq_len(n_spots), position_df$anno)
  
  for (i in seq_len(n_spots)) {
    r0 <- position_df$row[i]
    c0 <- position_df$col[i]
    ## For each cluster, compute minimal distance to any spot in that cluster
    for (cl in clusters) {
      idx <- cluster_index[[cl]]
      rr  <- position_df$row[idx]
      cc  <- position_df$col[idx]
      # Vectorized distance to all spots in this cluster
      dists <- sqrt((rr - r0)^2 + (cc - c0)^2)
      # Minimal distance (including zero if same cluster and same coordinates)
      out_mat[i, cl] <- min(dists)
    }
  }
  as.data.frame(out_mat, stringsAsFactors = FALSE)
}

## -------------------------------
## 3) Run computation
## -------------------------------
min_distances_df <- calculate_min_distances(position)

## Attach identifiers for convenience
min_distances_df$bar    <- rownames(min_distances_df)
min_distances_df$sample <- sample_id

## -------------------------------
## 4) Save
## -------------------------------
save(min_distances_df, file = out_file)
message("[OK] Saved: ", out_file)

## -------------------------------
## 5) (Optional) sanity prints
## -------------------------------
message("[INFO] Spots: ", nrow(position),
        " | Clusters: ", length(unique(position$anno)))
