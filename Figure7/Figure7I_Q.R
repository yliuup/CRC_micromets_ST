# ============================================================
# Survival analysis for in-house RFS cohort
# Signature: RNF40, AEN, WEE1, BCL7B, YME1L1, COX17
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

# -----------------------------
# 1. Set paths and load data
# -----------------------------
work_dir <- "/Users/yliu"
setwd(work_dir)

load("survival_data.RData")

output_dir <- file.path(work_dir, "survival_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Define helper functions
# -----------------------------
z_score <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# -----------------------------
# 3. Define gene signature
# -----------------------------
signature_genes <- c(
  "RNF40", "AEN", "WEE1",
  "BCL7B", "YME1L1", "COX17"
)

normalize_by_gene <- "EPCAM"

# Check gene availability
missing_genes <- setdiff(signature_genes, rownames(ff2))
if (length(missing_genes) > 0) {
  warning("Missing signature genes: ", paste(missing_genes, collapse = ", "))
}

signature_genes <- intersect(signature_genes, rownames(ff2))

if (!normalize_by_gene %in% rownames(ff2)) {
  stop("Normalization gene not found in expression matrix: ", normalize_by_gene)
}

# -----------------------------
# 4. Calculate signature score
# -----------------------------
expr_signature <- ff2[signature_genes, , drop = FALSE] %>%
  na.omit()

expr_signature_z <- t(apply(expr_signature, 1, z_score))

signature_score <- colMeans(expr_signature_z, na.rm = TRUE)

# Normalize signature score by EPCAM expression
epcam_score <- as.numeric(z_score(ff2[normalize_by_gene, ]))

signature_score_norm <- signature_score - epcam_score

clinical_filt2$signature_score <- signature_score_norm

# -----------------------------
# 5. Ternary grouping
# 0 = low, 1 = intermediate, 2 = high
# -----------------------------
q40 <- quantile(signature_score_norm, 0.40, na.rm = TRUE)
q60 <- quantile(signature_score_norm, 0.60, na.rm = TRUE)

clinical_filt2$signature_group <- cut(
  signature_score_norm,
  breaks = c(-Inf, q40, q60, Inf),
  labels = c("Low", "Intermediate", "High"),
  include.lowest = TRUE
)

clinical_filt2$signature_group <- factor(
  clinical_filt2$signature_group,
  levels = c("Low", "Intermediate", "High")
)

# -----------------------------
# 6. Kaplan-Meier survival analysis
# -----------------------------
km_fit <- survfit(
  Surv(RFS_months, Relapse_event) ~ signature_group,
  data = clinical_filt2,
  conf.type = "log-log"
)

km_plot <- ggsurvplot(
  km_fit,
  data = clinical_filt2,
  conf.int = FALSE,
  pval = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.20,
  legend.title = "Signature",
  legend.labs = c("Low", "Intermediate", "High"),
  palette = c("blue", "green", "red"),
  title = "Kaplan-Meier Curve for Survival",
  xlab = "Time (months)",
  ylab = "Survival probability"
)

pdf(
  file.path(output_dir, "survival_signature_KM_ternary.pdf"),
  height = 4.9,
  width = 3.1
)
print(km_plot)
dev.off()

# -----------------------------
# 7. Cox proportional hazards model
# -----------------------------
cox_fit <- coxph(
  Surv(RFS_months, Relapse_event) ~
    signature_score +
    Age +
    Pathology.Response..Complete.0..Major.1..Minor.2. +
    RAS...braf,
  data = clinical_filt2
)

cox_summary <- summary(cox_fit)

# Save Cox summary
capture.output(
  cox_summary,
  file = file.path(output_dir, "survival_signature_Cox_summary.txt")
)

# -----------------------------
# 8. Forest plot
# -----------------------------
pdf(
  file.path(output_dir, "survival_signature_Cox_forest.pdf"),
  height = 4.8,
  width = 6
)
ggforest(cox_fit, data = clinical_filt2, fontsize = 1.1)
dev.off()

# -----------------------------
# 9. Optional: proportional hazards test
# -----------------------------
ph_test <- cox.zph(cox_fit)

capture.output(
  ph_test,
  file = file.path(output_dir, "survival_signature_Cox_PH_test.txt")
)