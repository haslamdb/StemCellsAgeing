#-------------------------------------------------------------------------------
#
# Script: Batch-Corrected ALDEx2 and DESeq2 Pathway Analysis
#
# Author:  David Haslam
# Date:    February 2026
#
# Description:
#   Re-run ALDEx2 and DESeq2 pathway analyses for DY vs DO (and Y vs O)
#   with Experiment (batch) as a covariate to match the batch-corrected
#   MaAsLin2 analysis.
#
#   ALDEx2: Uses aldex.clr() + aldex.glm() with model ~ Experiment + Groups
#   DESeq2: Uses design ~ Experiment + Groups
#
#-------------------------------------------------------------------------------

#===============================================================================
# SECTION 0: SETUP
#===============================================================================

library(ALDEx2)
library(DESeq2)
library(ggplot2)
library(dplyr)

# Load workspace (has Metadata with Experiment column)
if (file.exists("GeigerData_ReviewerRevisions_20251225")) {
  load("GeigerData_ReviewerRevisions_20251225")
} else {
  load("GeigerData20250512")
}

cat("\n========================================\n")
cat("Batch-Corrected ALDEx2 + DESeq2 Analysis\n")
cat("========================================\n\n")

#===============================================================================
# SECTION 1: PREPARE DATA
#===============================================================================

cat("=== Section 1: Preparing Data ===\n\n")

# Load pathway data
pathway_file <- "GeigerFiles_pathabundance20220910-cpm_unstratified.tsv"
PathwayTable <- read.csv(pathway_file, stringsAsFactors = FALSE,
                          header = TRUE, sep = "\t")
names(PathwayTable)[1] <- "Pathway"
colnames(PathwayTable) <- gsub(".paired_Abundance.CPM", "", colnames(PathwayTable))

PathwayTable <- PathwayTable[!grepl("^#|UNMAPPED|UNINTEGRATED", PathwayTable$Pathway), ]
row.names(PathwayTable) <- PathwayTable$Pathway
PathwayTable$Pathway <- NULL
colnames(PathwayTable) <- gsub("^X", "", colnames(PathwayTable))

# Prepare metadata
sample_meta <- Metadata[Metadata$SampleID %in% colnames(PathwayTable), ]
sample_meta$Groups <- as.character(sample_meta$Groups)
sample_meta$Groups <- gsub("Young_C57BL-Untransplanted", "Y", sample_meta$Groups)
sample_meta$Groups <- gsub("Old_C57BL-Untransplanted", "O", sample_meta$Groups)
sample_meta$Groups <- gsub("Young_RAG1-Untransplanted", "RAG1", sample_meta$Groups)
sample_meta$Groups <- gsub("RAG1-Young_HSCs", "DY", sample_meta$Groups)
sample_meta$Groups <- gsub("RAG1-Old_HSCs", "DO", sample_meta$Groups)
sample_meta$Groups <- gsub("RAG1-Rejuvenated_HSCs", "Rejuv", sample_meta$Groups)

sample_meta <- sample_meta[sample_meta$Groups %in% c("Y", "O", "RAG1", "DY", "DO"), ]
rownames(sample_meta) <- sample_meta$SampleID
sample_meta$Experiment <- factor(sample_meta$Experiment)

# Match samples
common_samples <- intersect(colnames(PathwayTable), rownames(sample_meta))
PathwayTable <- PathwayTable[, common_samples]
sample_meta <- sample_meta[common_samples, ]

cat("Batch structure:\n")
print(table(sample_meta$Groups, sample_meta$Experiment))

#===============================================================================
# SECTION 2: ALDEx2 BATCH-CORRECTED - DY vs DO
#===============================================================================

cat("\n=== Section 2: ALDEx2 Batch-Corrected - DY vs DO ===\n\n")

# Subset to DY and DO
dy_do_samples <- sample_meta$SampleID[sample_meta$Groups %in% c("DY", "DO")]
dy_do_meta <- sample_meta[dy_do_samples, ]
dy_do_meta$Groups <- factor(dy_do_meta$Groups, levels = c("DO", "DY"))

# Prepare counts for ALDEx2 (features as rows, samples as columns)
scale_factor <- 100
pathway_counts <- as.data.frame(lapply(PathwayTable[, dy_do_samples], function(x) {
  as.integer(round(as.numeric(x) * scale_factor))
}))
rownames(pathway_counts) <- rownames(PathwayTable)

# Filter low-count pathways
min_total <- 100
pathway_counts <- pathway_counts[rowSums(pathway_counts, na.rm = TRUE) >= min_total, ]
pathway_counts[is.na(pathway_counts)] <- 0

cat("ALDEx2 input:", nrow(pathway_counts), "pathways x", ncol(pathway_counts), "samples\n")
cat("DY:", sum(dy_do_meta$Groups == "DY"), "  DO:", sum(dy_do_meta$Groups == "DO"), "\n\n")

# --- 2a: ALDEx2 without batch (original, for comparison) ---
cat("Running ALDEx2 DY vs DO (NO batch correction)...\n")
conds_DY_DO <- as.character(dy_do_meta$Groups)
aldex_DY_DO_noBatch <- aldex(pathway_counts, conds_DY_DO,
                              mc.samples = 128, test = "t",
                              effect = TRUE, verbose = TRUE, denom = "all")
aldex_DY_DO_noBatch$Pathway <- rownames(aldex_DY_DO_noBatch)

# --- 2b: ALDEx2 with batch correction using glm ---
cat("\nRunning ALDEx2 DY vs DO (WITH batch correction via glm)...\n")

# Build model matrix: ~ Experiment + Groups
mm_DY_DO <- model.matrix(~ Experiment + Groups, data = dy_do_meta)
cat("Model matrix columns:", colnames(mm_DY_DO), "\n")

# Generate CLR instances - must pass model matrix as conditions for aldex.glm
clr_DY_DO <- aldex.clr(pathway_counts, mm_DY_DO,
                         mc.samples = 128, denom = "all", verbose = TRUE)

# Run GLM on CLR instances
glm_DY_DO <- aldex.glm(clr_DY_DO, mm_DY_DO)
glm_DY_DO$Pathway <- rownames(glm_DY_DO)

# The Groups effect column name depends on model matrix
# It should be something like "GroupsDY.pval" or "model.GroupsDY Pr(>|t|).BH"
cat("\nGLM result columns:\n")
print(colnames(glm_DY_DO))

# Find the Groups DY coefficient columns
groups_cols <- grep("GroupsDY", colnames(glm_DY_DO), value = TRUE)
cat("\nColumns for Groups effect:\n")
print(groups_cols)

# Extract p-values for the Groups effect
# ALDEx2 GLM column names: "GroupsDY:Est", "GroupsDY:pval", "GroupsDY:pval.padj"
est_col <- "GroupsDY:Est"
pval_col <- "GroupsDY:pval"
bh_col <- "GroupsDY:pval.padj"

cat("\nEstimate column:", est_col, "\n")
cat("P-value column:", pval_col, "\n")
cat("BH column:", bh_col, "\n")

# Create clean results dataframe
aldex_batch_DY_DO <- data.frame(
  Pathway = rownames(glm_DY_DO),
  Estimate = as.numeric(glm_DY_DO[, est_col]),
  pval = as.numeric(glm_DY_DO[, pval_col]),
  pval_BH = as.numeric(glm_DY_DO[, bh_col]),
  stringsAsFactors = FALSE
)

aldex_batch_DY_DO <- aldex_batch_DY_DO[order(aldex_batch_DY_DO$pval), ]

cat("\nTop 20 pathways (ALDEx2 batch-corrected, DY vs DO):\n")
print(head(aldex_batch_DY_DO[, c("Pathway", "Estimate", "pval", "pval_BH")], 20))

# B6 pathways
b6_patterns <- c("PYRIDOX", "pyridox", "pyridoxal", "PLP", "PWY0.845", "PWY0-845")
b6_idx <- grepl(paste(b6_patterns, collapse = "|"), aldex_batch_DY_DO$Pathway, ignore.case = TRUE)

cat("\n--- B6 Pathways: ALDEx2 DY vs DO (batch-corrected) ---\n")
if (any(b6_idx)) {
  print(aldex_batch_DY_DO[b6_idx, ])
  write.csv(aldex_batch_DY_DO[b6_idx, ],
            "ALDEx2_B6_Pathways_DY_vs_DO_BatchCorrected.csv", row.names = FALSE)
}

# Also show the original (no batch) B6 results for comparison
b6_orig <- aldex_DY_DO_noBatch[grepl(paste(b6_patterns, collapse = "|"),
                                       aldex_DY_DO_noBatch$Pathway, ignore.case = TRUE), ]
cat("\n--- B6 Pathways: ALDEx2 DY vs DO (NO batch correction, for comparison) ---\n")
if (nrow(b6_orig) > 0) {
  print(b6_orig[, c("Pathway", "effect", "we.ep", "we.eBH", "wi.ep", "wi.eBH")])
}

# Save full results
write.csv(aldex_batch_DY_DO, "ALDEx2_Pathways_DY_vs_DO_BatchCorrected_Full.csv", row.names = FALSE)
write.csv(glm_DY_DO, "ALDEx2_GLM_DY_vs_DO_BatchCorrected_AllColumns.csv")


#===============================================================================
# SECTION 3: DESeq2 BATCH-CORRECTED - DY vs DO
#===============================================================================

cat("\n=== Section 3: DESeq2 Batch-Corrected - DY vs DO ===\n\n")

# Prepare count matrix for DESeq2
count_matrix <- round(as.matrix(PathwayTable[, dy_do_samples]) * 100)
count_matrix[is.na(count_matrix)] <- 0
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

# Filter low-count pathways
keep <- rowSums(count_matrix >= 10) >= 3
count_matrix <- count_matrix[keep, ]

cat("DESeq2 input:", nrow(count_matrix), "pathways x", ncol(count_matrix), "samples\n\n")

# --- 3a: Without batch (original) ---
cat("Running DESeq2 DY vs DO (NO batch correction)...\n")
dds_noBatch <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = dy_do_meta,
  design = ~ Groups
)
dds_noBatch <- DESeq(dds_noBatch)
res_noBatch <- results(dds_noBatch, contrast = c("Groups", "DY", "DO"))
res_noBatch <- as.data.frame(res_noBatch[order(res_noBatch$pvalue), ])
res_noBatch$Pathway <- rownames(res_noBatch)

# --- 3b: With batch correction ---
cat("Running DESeq2 DY vs DO (WITH batch correction: ~ Experiment + Groups)...\n")
dds_batch <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = dy_do_meta,
  design = ~ Experiment + Groups
)
dds_batch <- DESeq(dds_batch)
res_batch <- results(dds_batch, contrast = c("Groups", "DY", "DO"))
res_batch <- as.data.frame(res_batch[order(res_batch$pvalue), ])
res_batch$Pathway <- rownames(res_batch)

cat("\nDESeq2 DY vs DO (batch-corrected) summary:\n")
cat("  Significant (padj < 0.1):", sum(res_batch$padj < 0.1, na.rm = TRUE), "\n")
cat("  Significant (padj < 0.05):", sum(res_batch$padj < 0.05, na.rm = TRUE), "\n")
cat("  Nominally significant (p < 0.05):", sum(res_batch$pvalue < 0.05, na.rm = TRUE), "\n")

cat("\nTop 20 pathways (DESeq2 batch-corrected, DY vs DO):\n")
print(head(res_batch[, c("Pathway", "log2FoldChange", "pvalue", "padj")], 20))

# B6 pathways
b6_idx_deseq <- grepl(paste(b6_patterns, collapse = "|"), res_batch$Pathway, ignore.case = TRUE)

cat("\n--- B6 Pathways: DESeq2 DY vs DO (batch-corrected) ---\n")
if (any(b6_idx_deseq)) {
  b6_deseq_batch <- res_batch[b6_idx_deseq, ]
  print(b6_deseq_batch[, c("Pathway", "log2FoldChange", "lfcSE", "pvalue", "padj")])
  write.csv(b6_deseq_batch, "DESeq2_B6_Pathways_DY_vs_DO_BatchCorrected.csv", row.names = FALSE)
}

# Original (no batch) B6 results for comparison
b6_deseq_orig <- res_noBatch[grepl(paste(b6_patterns, collapse = "|"),
                                    res_noBatch$Pathway, ignore.case = TRUE), ]
cat("\n--- B6 Pathways: DESeq2 DY vs DO (NO batch correction, for comparison) ---\n")
if (nrow(b6_deseq_orig) > 0) {
  print(b6_deseq_orig[, c("Pathway", "log2FoldChange", "lfcSE", "pvalue", "padj")])
}

# Save full results
write.csv(res_batch, "DESeq2_Pathways_DY_vs_DO_BatchCorrected_Full.csv", row.names = FALSE)
write.csv(res_noBatch, "DESeq2_Pathways_DY_vs_DO_NoBatch_Full.csv", row.names = FALSE)


#===============================================================================
# SECTION 4: UPDATED CROSS-METHOD COMPARISON
#===============================================================================

cat("\n=== Section 4: Updated Cross-Method Comparison for B6 Pathways (DY vs DO) ===\n\n")

# Helper to extract p-values
get_b6_p <- function(df, pathway_col, pval_col, pattern) {
  idx <- grep(pattern, df[[pathway_col]], ignore.case = TRUE)
  if (length(idx) > 0) return(df[[pval_col]][idx[1]])
  return(NA)
}

# Build comparison table
cat("Cross-method comparison for B6 pathways (DY vs DO):\n\n")
cat(sprintf("%-30s | %-18s | %-18s\n", "Method", "PYRIDOXSYN-PWY p", "PWY0-845 p"))
cat(paste(rep("-", 70), collapse = ""), "\n")

# Wilcoxon (from original analysis)
cat(sprintf("%-30s | %-18s | %-18s\n", "Wilcoxon (original)", "0.024*", "0.027*"))

# ALDEx2 no batch
p1 <- get_b6_p(aldex_DY_DO_noBatch, "Pathway", "we.ep", "PYRIDOXSYN")
p2 <- get_b6_p(aldex_DY_DO_noBatch, "Pathway", "we.ep", "PWY0.845")
cat(sprintf("%-30s | %-18s | %-18s\n", "ALDEx2 Welch (no batch)",
            paste0(signif(p1, 3), ifelse(p1 < 0.05, "*", "")),
            paste0(signif(p2, 3), ifelse(p2 < 0.05, "*", ""))))

# ALDEx2 with batch
p1b <- get_b6_p(aldex_batch_DY_DO, "Pathway", "pval", "PYRIDOXSYN")
p2b <- get_b6_p(aldex_batch_DY_DO, "Pathway", "pval", "PWY0.845")
cat(sprintf("%-30s | %-18s | %-18s\n", "ALDEx2 GLM (batch-corrected)",
            paste0(signif(p1b, 3), ifelse(!is.na(p1b) && p1b < 0.05, "*", "")),
            paste0(signif(p2b, 3), ifelse(!is.na(p2b) && p2b < 0.05, "*", ""))))

# DESeq2 no batch
p1d <- get_b6_p(res_noBatch, "Pathway", "pvalue", "PYRIDOXSYN")
p2d <- get_b6_p(res_noBatch, "Pathway", "pvalue", "PWY0.845")
cat(sprintf("%-30s | %-18s | %-18s\n", "DESeq2 (no batch)",
            paste0(signif(p1d, 3), ifelse(!is.na(p1d) && p1d < 0.05, "*", "")),
            paste0(signif(p2d, 3), ifelse(!is.na(p2d) && p2d < 0.05, "*", ""))))

# DESeq2 with batch
p1db <- get_b6_p(res_batch, "Pathway", "pvalue", "PYRIDOXSYN")
p2db <- get_b6_p(res_batch, "Pathway", "pvalue", "PWY0.845")
cat(sprintf("%-30s | %-18s | %-18s\n", "DESeq2 (batch-corrected)",
            paste0(signif(p1db, 3), ifelse(!is.na(p1db) && p1db < 0.05, "*", "")),
            paste0(signif(p2db, 3), ifelse(!is.na(p2db) && p2db < 0.05, "*", ""))))

# MaAsLin2 (from previous run)
cat(sprintf("%-30s | %-18s | %-18s\n", "MaAsLin2 (no batch)", "0.0988", "0.131"))
cat(sprintf("%-30s | %-18s | %-18s\n", "MaAsLin2 (batch-corrected)", "0.00703*", "0.0114*"))

# Save comparison table
comparison_df <- data.frame(
  Method = c("Wilcoxon (original)",
             "ALDEx2 Welch (no batch)", "ALDEx2 GLM (batch-corrected)",
             "DESeq2 (no batch)", "DESeq2 (batch-corrected)",
             "MaAsLin2 (no batch)", "MaAsLin2 (batch-corrected)"),
  PYRIDOXSYN_PWY_p = c(0.024, p1, p1b, p1d, p1db, 0.0988, 0.00703),
  PWY0_845_p = c(0.027, p2, p2b, p2d, p2db, 0.131, 0.0114),
  Batch_Corrected = c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)
comparison_df$PYRIDOXSYN_sig <- comparison_df$PYRIDOXSYN_PWY_p < 0.05
comparison_df$PWY0_845_sig <- comparison_df$PWY0_845_p < 0.05

write.csv(comparison_df, "B6_CrossMethod_Comparison_BatchEffect.csv", row.names = FALSE)

cat("\n\nMethods with nominally significant B6 pathways (p < 0.05):\n")
cat("  Without batch correction: ")
no_batch <- comparison_df[!comparison_df$Batch_Corrected, ]
cat(sum(no_batch$PYRIDOXSYN_sig), "/ ", nrow(no_batch), " methods\n")
cat("  With batch correction:    ")
with_batch <- comparison_df[comparison_df$Batch_Corrected, ]
cat(sum(with_batch$PYRIDOXSYN_sig, na.rm = TRUE), "/", nrow(with_batch), "methods\n")


#===============================================================================
# SECTION 5: SUMMARY
#===============================================================================

cat("\n\n========================================\n")
cat("Batch-Corrected Analysis Complete!\n")
cat("========================================\n\n")

cat("Output files:\n")
cat("  ALDEx2_B6_Pathways_DY_vs_DO_BatchCorrected.csv\n")
cat("  ALDEx2_Pathways_DY_vs_DO_BatchCorrected_Full.csv\n")
cat("  ALDEx2_GLM_DY_vs_DO_BatchCorrected_AllColumns.csv\n")
cat("  DESeq2_B6_Pathways_DY_vs_DO_BatchCorrected.csv\n")
cat("  DESeq2_Pathways_DY_vs_DO_BatchCorrected_Full.csv\n")
cat("  DESeq2_Pathways_DY_vs_DO_NoBatch_Full.csv\n")
cat("  B6_CrossMethod_Comparison_BatchEffect.csv\n")

sessionInfo()
