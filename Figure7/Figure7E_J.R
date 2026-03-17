## -------------------------------
## Six-gene signature survival pipeline (clean & robust)
## -------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
  library(survival)
  library(survminer)
  library(ggplot2)
})

# --- helpers ---
ZScore <- function(x){
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x,  na.rm = TRUE)
  if (is.na(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}
as_binary01 <- function(x){
  if (is.logical(x)) return(as.integer(x))
  x <- as.numeric(x); x[is.na(x)] <- 0; x[x != 0] <- 1L; as.integer(x)
}

# --- inputs expected ---
# dat:    gene x sample normalized matrix (rownames=genes, colnames=samples)
# clinical: data.frame with columns: sample, RFS_months, Relapse_event, covariates...
stopifnot(exists("dat"), is.matrix(dat) || is.data.frame(dat))
dat <- as.matrix(dat)
stopifnot(!is.null(rownames(dat)), !is.null(colnames(dat)))

stopifnot(exists("clinical"), is.data.frame(clinical))
stopifnot(all(c("sample","RFS_months","Relapse_event") %in% colnames(clinical)))

# --- genes & reference gene ---
gene_used <- c("RNF40","AEN","WEE1","BCL7B","YME1L1","COX17")
normalize.by.gene <- "EPCAM"   # set to NULL to skip normalization

# --- subset genes present ---
gene_used <- intersect(gene_used, rownames(dat))
if (length(gene_used) == 0) stop("None of the six signature genes are present in 'dat'.")
if (length(gene_used) < 6) warning("Missing signature genes: ",
                                   paste(setdiff(c("RNF40","AEN","WEE1","BCL7B","YME1L1","COX17"), gene_used), collapse=", "))

# --- sample alignment (expression ↔ clinical) ---
common_samples <- intersect(colnames(dat), clinical$sample)
if (length(common_samples) < 20) warning("Only ", length(common_samples), " overlapping samples.")
dat <- dat[, common_samples, drop = FALSE]
clinical_filt <- clinical %>% filter(sample %in% common_samples) %>% distinct(sample, .keep_all = TRUE)
# reorder clinical to match expression columns
clinical_filt <- clinical_filt[match(colnames(dat), clinical_filt$sample), ]

# --- z-score per gene, then mean per sample ---
dat_nor <- dat[gene_used, , drop = FALSE]
dat_nor <- dat_nor[rowSums(is.na(dat_nor)) < ncol(dat_nor), , drop = FALSE]  # drop all-NA genes if any
expression_dataset_zscore <- t(apply(dat_nor, 1, ZScore))   # gene x sample
expression_dataset_zscore_mean <- colMeans(expression_dataset_zscore, na.rm = TRUE)

# --- optional normalization by a reference gene (e.g., EPCAM) ---
if (!is.null(normalize.by.gene)) {
  if (!(normalize.by.gene %in% rownames(dat))) {
    warning("normalize.by.gene '", normalize.by.gene, "' not in 'dat'; skipping normalization.")
  } else {
    ref_vec_z <- ZScore(dat[normalize.by.gene, colnames(dat), drop = TRUE])
    expression_dataset_zscore_mean <- expression_dataset_zscore_mean - as.vector(ref_vec_z)
  }
}

# --- ROC to get best threshold (use clinical_filt) ---
clinical_filt$Relapse_event <- as_binary01(clinical_filt$Relapse_event)

roc.module <- pROC::roc(response = clinical_filt$Relapse_event,
                        predictor = as.numeric(expression_dataset_zscore_mean),
                        quiet = TRUE)
auc_val <- as.numeric(pROC::auc(roc.module))

coords.modult <- pROC::coords(roc.module, "best",
                              input = "threshold",
                              ret   = c("threshold","specificity","sensitivity"),
                              as.list = FALSE, drop = TRUE,
                              best.method = "closest.topleft")

thr  <- as.numeric(coords.modult["threshold"])
sens <- as.numeric(coords.modult["sensitivity"])
spec <- as.numeric(coords.modult["specificity"])

message(sprintf("[ROC] AUC = %.3f | threshold = %.4f | sens = %.3f | spec = %.3f",
                auc_val, thr, sens, spec))

# --- dichotomize score by threshold ---
binary_expression <- as.integer(expression_dataset_zscore_mean > thr)

# attach to clinical_filt for KM/COX
clinical_filt$genes_exp <- binary_expression
clinical_filt$exp       <- as.numeric(expression_dataset_zscore_mean)

# --- KM (RFS) ---
fit_km <- with(clinical_filt,
               survfit(Surv(RFS_months, Relapse_event) ~ genes_exp,
                       conf.type = "log-log"))
plt_km <- ggsurvplot(fit_km,
                     data = clinical_filt,
                     conf.int = FALSE, pval = TRUE, risk.table = TRUE,
                     legend.labs = c("low","high"),
                     legend.title = "signature",
                     palette = c("blue","red"),
                     title = "Kaplan-Meier: RFS by six-gene signature",
                     risk.table.height = .20)
print(plt_km)

# --- Cox (multivariable) ---
# NOTE: make sure these covariates exist in clinical_filt; rename if needed.
cox_formula <- as.formula(
  paste(
    "Surv(RFS_months, Relapse_event) ~ exp +",
    "Age + sex +",
    "Primary..1.Right..2.Left. +",
    "Primary.Node.Positive..Yes.1..No.0. +",
    "Synchronous.1..Metachronous.2 +",
    "Bilateral.Liver.Mets..Yes.1..No.0. +",
    "Exrahepatic.disease..Yes.1..No.0. +",
    "Hepatectomy..Major.1..minor.oor.other.2. +",
    "R1.resection..Yes.1..No..0. +",
    "Pathology.Response..Complete.0..Major.1..Minor.2. +",
    "RAS...braf"
  )
)

vars_needed <- all.vars(cox_formula)
cox_df <- clinical_filt[, vars_needed, drop = TRUE] %>% na.omit()
message("[Cox] N used (complete cases) = ", nrow(cox_df))

a <- coxph(cox_formula, data = cox_df)
plt_forest <- ggforest(a, data = cox_df, fontsize = 1.3)
print(plt_forest)

# --- Optional: save objects ---
# saveRDS(list(
#   roc = roc.module, auc = auc_val, thr = thr, sens = sens, spec = spec,
#   km_fit = fit_km, km_plot = plt_km,
#   cox_fit = a, cox_plot = plt_forest,
#   clinical_used = clinical_filt,
#   signature_score = expression_dataset_zscore_mean
# ), file = "six_gene_signature_survival_results.rds")
