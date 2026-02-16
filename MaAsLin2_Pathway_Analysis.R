#-------------------------------------------------------------------------------
#
# Script: MaAsLin2 Pathway Differential Abundance Analysis
#
# Author:  David Haslam
# Date:    February 2026
#
# Description:
#   MaAsLin2 analysis for pathway differential abundance (DY vs DO and Y vs O)
#   with proper FDR correction. MaAsLin2 is designed for microbiome data and
#   works directly with HUMAnN3 CPM-normalized pathway abundance.
#
#   Advantages over DESeq2 for this data:
#   - Designed for microbiome abundance data (not RNA-seq)
#   - Works directly with CPM values (no pseudo-count conversion needed)
#   - Can include covariates (e.g., batch/experiment)
#   - Multiple normalization options (TSS, CLR, none for pre-normalized data)
#   - Linear model framework well-suited to continuous abundance data
#
#-------------------------------------------------------------------------------

#===============================================================================
# SECTION 0: SETUP
#===============================================================================

# Install Maaslin2 if needed
if (!require("Maaslin2", quietly = TRUE)) {
  BiocManager::install("Maaslin2")
}

library(Maaslin2)
library(dplyr)
library(ggplot2)

# Load workspace
if (file.exists("GeigerData_ReviewerRevisions_20251225")) {
  load("GeigerData_ReviewerRevisions_20251225")
} else {
  load("GeigerData20250512")
}

cat("\n========================================\n")
cat("MaAsLin2 Pathway Analysis\n")
cat("========================================\n\n")

#===============================================================================
# SECTION 1: PREPARE DATA
#===============================================================================

cat("=== Section 1: Preparing Data ===\n\n")

# Load HUMAnN3 pathway abundance (CPM-normalized, unstratified)
pathway_file <- "GeigerFiles_pathabundance20220910-cpm_unstratified.tsv"
PathwayTable <- read.csv(pathway_file, stringsAsFactors = FALSE,
                          header = TRUE, sep = "\t")
names(PathwayTable)[1] <- "Pathway"
colnames(PathwayTable) <- gsub(".paired_Abundance.CPM", "", colnames(PathwayTable))

# Clean: remove UNMAPPED, UNINTEGRATED, and comment lines
PathwayTable <- PathwayTable[!grepl("^#|UNMAPPED|UNINTEGRATED", PathwayTable$Pathway), ]
rownames(PathwayTable) <- PathwayTable$Pathway
PathwayTable$Pathway <- NULL
colnames(PathwayTable) <- gsub("^X", "", colnames(PathwayTable))

# Prepare metadata
sample_meta <- Metadata[Metadata$SampleID %in% colnames(PathwayTable), ]

# Standardize group names (match what's used in the analysis)
sample_meta$Groups <- as.character(sample_meta$Groups)
sample_meta$Groups <- gsub("Young_C57BL-Untransplanted", "Y", sample_meta$Groups)
sample_meta$Groups <- gsub("Old_C57BL-Untransplanted", "O", sample_meta$Groups)
sample_meta$Groups <- gsub("Young_RAG1-Untransplanted", "RAG1", sample_meta$Groups)
sample_meta$Groups <- gsub("RAG1-Young_HSCs", "DY", sample_meta$Groups)
sample_meta$Groups <- gsub("RAG1-Old_HSCs", "DO", sample_meta$Groups)
sample_meta$Groups <- gsub("RAG1-Rejuvenated_HSCs", "Rejuv", sample_meta$Groups)

# Filter to main groups
sample_meta <- sample_meta[sample_meta$Groups %in% c("Y", "O", "RAG1", "DY", "DO"), ]
rownames(sample_meta) <- sample_meta$SampleID

# Match columns
common_samples <- intersect(colnames(PathwayTable), rownames(sample_meta))
PathwayTable <- PathwayTable[, common_samples]
sample_meta <- sample_meta[common_samples, ]

cat("Data dimensions:", nrow(PathwayTable), "pathways x", ncol(PathwayTable), "samples\n")
cat("Sample sizes:\n")
print(table(sample_meta$Groups))

# MaAsLin2 needs: samples as ROWS, features as COLUMNS
pathway_for_maaslin <- as.data.frame(t(PathwayTable))


#===============================================================================
# SECTION 2: MaAsLin2 ANALYSIS - DY vs DO
#===============================================================================

cat("\n=== Section 2: MaAsLin2 - DY vs DO ===\n\n")

# Subset to DY and DO
dy_do_samples <- sample_meta$SampleID[sample_meta$Groups %in% c("DY", "DO")]
dy_do_pathways <- pathway_for_maaslin[dy_do_samples, ]
dy_do_meta <- sample_meta[dy_do_samples, ]
dy_do_meta$Groups <- factor(dy_do_meta$Groups, levels = c("DO", "DY"))  # DO = reference

# Run MaAsLin2 - basic model (no covariates)
cat("Running MaAsLin2 DY vs DO (basic model)...\n")
maaslin_DY_DO <- Maaslin2(
  input_data = dy_do_pathways,
  input_metadata = dy_do_meta,
  output = "MaAsLin2_DY_vs_DO",
  fixed_effects = c("Groups"),
  reference = c("Groups,DO"),
  normalization = "NONE",         # Data is already CPM-normalized
  transform = "LOG",              # Log-transform CPM values
  analysis_method = "LM",         # Linear model
  min_abundance = 0.0,            # Include all pathways
  min_prevalence = 0.1,           # Present in at least 10% of samples
  max_significance = 0.25,        # Report results up to q=0.25
  correction = "BH",              # Benjamini-Hochberg FDR
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  cores = 1
)

# Extract and examine B6 results
maaslin_results_DY_DO <- maaslin_DY_DO$results
cat("\nDY vs DO: Total features tested:", nrow(maaslin_results_DY_DO), "\n")
cat("Significant (q < 0.25):", sum(maaslin_results_DY_DO$qval < 0.25, na.rm = TRUE), "\n")
cat("Significant (q < 0.1):", sum(maaslin_results_DY_DO$qval < 0.1, na.rm = TRUE), "\n")
cat("Nominally significant (p < 0.05):", sum(maaslin_results_DY_DO$pval < 0.05, na.rm = TRUE), "\n")

# B6-related pathways
b6_patterns <- c("PYRIDOX", "pyridox", "pyridoxal", "PLP")
b6_idx <- grepl(paste(b6_patterns, collapse = "|"), maaslin_results_DY_DO$feature,
                ignore.case = TRUE)

cat("\n--- B6 Pathways in DY vs DO ---\n")
if (any(b6_idx)) {
  b6_results_DY_DO <- maaslin_results_DY_DO[b6_idx, ]
  print(b6_results_DY_DO[, c("feature", "coef", "stderr", "pval", "qval")])
  write.csv(b6_results_DY_DO, "MaAsLin2_B6_Pathways_DY_vs_DO.csv", row.names = FALSE)
} else {
  cat("No B6 pathways found in results.\n")
}

# Top 20 results by p-value
cat("\n--- Top 20 Pathways (DY vs DO) ---\n")
top_DY_DO <- head(maaslin_results_DY_DO[order(maaslin_results_DY_DO$pval), ], 20)
print(top_DY_DO[, c("feature", "coef", "pval", "qval")])

# Save full results
write.csv(maaslin_results_DY_DO, "MaAsLin2_Pathways_DY_vs_DO_Full.csv", row.names = FALSE)

# Significant results
sig_DY_DO <- maaslin_results_DY_DO[maaslin_results_DY_DO$qval < 0.25, ]
write.csv(sig_DY_DO, "MaAsLin2_Pathways_DY_vs_DO_Significant.csv", row.names = FALSE)


#===============================================================================
# SECTION 2B: MaAsLin2 - DY vs DO WITH BATCH COVARIATE
#===============================================================================

cat("\n=== Section 2B: MaAsLin2 - DY vs DO (with batch covariate) ===\n\n")

# Check if Experiment column exists for batch correction
if ("Experiment" %in% colnames(dy_do_meta)) {
  dy_do_meta$Experiment <- factor(dy_do_meta$Experiment)

  cat("Running MaAsLin2 DY vs DO (with batch covariate)...\n")
  maaslin_DY_DO_batch <- Maaslin2(
    input_data = dy_do_pathways,
    input_metadata = dy_do_meta,
    output = "MaAsLin2_DY_vs_DO_BatchCorrected",
    fixed_effects = c("Groups", "Experiment"),
    reference = c("Groups,DO"),
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    min_abundance = 0.0,
    min_prevalence = 0.1,
    max_significance = 0.25,
    correction = "BH",
    plot_heatmap = TRUE,
    plot_scatter = TRUE,
    cores = 1
  )

  maaslin_results_DY_DO_batch <- maaslin_DY_DO_batch$results
  # Filter to Groups effect only (not Experiment)
  group_results <- maaslin_results_DY_DO_batch[maaslin_results_DY_DO_batch$metadata == "Groups", ]

  cat("\nDY vs DO (batch-corrected): Total features:", nrow(group_results), "\n")
  cat("Significant (q < 0.25):", sum(group_results$qval < 0.25, na.rm = TRUE), "\n")
  cat("Nominally significant (p < 0.05):", sum(group_results$pval < 0.05, na.rm = TRUE), "\n")

  b6_idx_batch <- grepl(paste(b6_patterns, collapse = "|"), group_results$feature,
                         ignore.case = TRUE)
  if (any(b6_idx_batch)) {
    cat("\n--- B6 Pathways (batch-corrected) ---\n")
    print(group_results[b6_idx_batch, c("feature", "coef", "stderr", "pval", "qval")])
    write.csv(group_results[b6_idx_batch, ],
              "MaAsLin2_B6_Pathways_DY_vs_DO_BatchCorrected.csv", row.names = FALSE)
  }

  write.csv(group_results, "MaAsLin2_Pathways_DY_vs_DO_BatchCorrected_Full.csv", row.names = FALSE)
} else {
  cat("Experiment column not found. Skipping batch-corrected analysis.\n")
  cat("Check column names in metadata:\n")
  print(colnames(dy_do_meta))
}


#===============================================================================
# SECTION 3: MaAsLin2 ANALYSIS - Y vs O
#===============================================================================

cat("\n=== Section 3: MaAsLin2 - Y vs O ===\n\n")

y_o_samples <- sample_meta$SampleID[sample_meta$Groups %in% c("Y", "O")]
y_o_pathways <- pathway_for_maaslin[y_o_samples, ]
y_o_meta <- sample_meta[y_o_samples, ]
y_o_meta$Groups <- factor(y_o_meta$Groups, levels = c("O", "Y"))

cat("Running MaAsLin2 Y vs O...\n")
maaslin_Y_O <- Maaslin2(
  input_data = y_o_pathways,
  input_metadata = y_o_meta,
  output = "MaAsLin2_Y_vs_O",
  fixed_effects = c("Groups"),
  reference = c("Groups,O"),
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM",
  min_abundance = 0.0,
  min_prevalence = 0.1,
  max_significance = 0.25,
  correction = "BH",
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  cores = 1
)

maaslin_results_Y_O <- maaslin_Y_O$results
cat("\nY vs O: Total features tested:", nrow(maaslin_results_Y_O), "\n")
cat("Significant (q < 0.25):", sum(maaslin_results_Y_O$qval < 0.25, na.rm = TRUE), "\n")
cat("Significant (q < 0.1):", sum(maaslin_results_Y_O$qval < 0.1, na.rm = TRUE), "\n")

b6_idx_yo <- grepl(paste(b6_patterns, collapse = "|"), maaslin_results_Y_O$feature,
                    ignore.case = TRUE)
if (any(b6_idx_yo)) {
  cat("\n--- B6 Pathways in Y vs O ---\n")
  print(maaslin_results_Y_O[b6_idx_yo, c("feature", "coef", "stderr", "pval", "qval")])
  write.csv(maaslin_results_Y_O[b6_idx_yo, ],
            "MaAsLin2_B6_Pathways_Y_vs_O.csv", row.names = FALSE)
}

write.csv(maaslin_results_Y_O, "MaAsLin2_Pathways_Y_vs_O_Full.csv", row.names = FALSE)


#===============================================================================
# SECTION 4: CROSS-METHOD COMPARISON TABLE FOR B6 PATHWAYS
#===============================================================================

cat("\n=== Section 4: Cross-Method Comparison for B6 Pathways ===\n\n")

# Build a comparison table across all methods for B6 pathways (DY vs DO)
b6_comparison <- data.frame(
  Pathway = c("PYRIDOXSYN-PWY: PLP biosynthesis I",
              "PWY0-845: Superpathway PLP biosynthesis"),
  stringsAsFactors = FALSE
)

# Add MaAsLin2 results
if (any(b6_idx)) {
  b6_maas <- maaslin_results_DY_DO[b6_idx, ]
  for (i in 1:nrow(b6_comparison)) {
    pattern <- ifelse(i == 1, "PYRIDOXSYN", "PWY0.845|PWY0-845")
    match_idx <- grep(pattern, b6_maas$feature, ignore.case = TRUE)
    if (length(match_idx) > 0) {
      b6_comparison$MaAsLin2_coef[i] <- round(b6_maas$coef[match_idx[1]], 3)
      b6_comparison$MaAsLin2_p[i] <- signif(b6_maas$pval[match_idx[1]], 3)
      b6_comparison$MaAsLin2_q[i] <- signif(b6_maas$qval[match_idx[1]], 3)
    }
  }
}

# Print summary
cat("Cross-method comparison for B6 pathways (DY vs DO):\n\n")
cat("Method         | PYRIDOXSYN-PWY p | PWY0-845 p | Notes\n")
cat("Wilcoxon       | 0.024*           | 0.027*     | Original analysis, nominally sig\n")
cat("ALDEx2 Welch   | 0.011*           | 0.012*     | CLR-based, effect 0.72-0.76\n")
cat("ALDEx2 Wilcox  | 0.020*           | 0.024*     | CLR-based\n")
cat("DESeq2         | 0.367            | 0.512      | NOT significant (pseudo-count issue?)\n")

if (exists("b6_comparison") && "MaAsLin2_p" %in% colnames(b6_comparison)) {
  cat(sprintf("MaAsLin2       | %-16s | %-10s | Linear model on log(CPM)\n",
              paste0(b6_comparison$MaAsLin2_p[1],
                     ifelse(b6_comparison$MaAsLin2_p[1] < 0.05, "*", "")),
              paste0(b6_comparison$MaAsLin2_p[2],
                     ifelse(b6_comparison$MaAsLin2_p[2] < 0.05, "*", ""))))
}

write.csv(b6_comparison, "B6_CrossMethod_Comparison_DY_vs_DO.csv", row.names = FALSE)


#===============================================================================
# SECTION 5: SUMMARY
#===============================================================================

cat("\n\n========================================\n")
cat("MaAsLin2 Analysis Complete!\n")
cat("========================================\n\n")

cat("Output files:\n")
cat("  MaAsLin2_Pathways_DY_vs_DO_Full.csv\n")
cat("  MaAsLin2_Pathways_DY_vs_DO_Significant.csv\n")
cat("  MaAsLin2_Pathways_Y_vs_O_Full.csv\n")
cat("  MaAsLin2_B6_Pathways_DY_vs_DO.csv\n")
cat("  MaAsLin2_B6_Pathways_Y_vs_O.csv\n")
cat("  B6_CrossMethod_Comparison_DY_vs_DO.csv\n")
if ("Experiment" %in% colnames(dy_do_meta)) {
  cat("  MaAsLin2_Pathways_DY_vs_DO_BatchCorrected_Full.csv\n")
  cat("  MaAsLin2_B6_Pathways_DY_vs_DO_BatchCorrected.csv\n")
}
cat("\n  MaAsLin2_DY_vs_DO/ (directory with plots)\n")
cat("  MaAsLin2_Y_vs_O/ (directory with plots)\n")
if ("Experiment" %in% colnames(dy_do_meta)) {
  cat("  MaAsLin2_DY_vs_DO_BatchCorrected/ (directory with plots)\n")
}

cat("\n--- Suggested M&M Text ---\n")
cat("Differential pathway abundance was additionally assessed using MaAsLin2\n")
cat("(Mallick et al., 2021), which employs generalized linear models on\n")
cat("log-transformed CPM-normalized pathway abundances with Benjamini-Hochberg\n")
cat("FDR correction. Experimental batch was included as a covariate where\n")
cat("applicable to account for potential confounding across sequencing runs.\n")

sessionInfo()
