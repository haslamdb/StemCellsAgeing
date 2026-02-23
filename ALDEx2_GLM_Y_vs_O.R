#-------------------------------------------------------------------------------
#
# Script: ALDEx2 GLM for Y vs O (to get proper SEs for log2FC plots)
#
# Author:  David Haslam
# Date:    February 2026
#
# Description:
#   Re-run the Y vs O ALDEx2 analysis using aldex.glm() instead of the
#   classic aldex() t-test approach. The GLM provides proper standard errors
#   for the group effect coefficient, which are needed for bar plots with
#   error bars. If Y and O span multiple batches, includes Experiment as
#   a covariate.
#
#-------------------------------------------------------------------------------

library(ALDEx2)
library(dplyr)

# Load workspace
if (file.exists("GeigerData_ReviewerRevisions_20251225")) {
  load("GeigerData_ReviewerRevisions_20251225")
} else {
  load("GeigerData20250512")
}

cat("\n========================================\n")
cat("ALDEx2 GLM: Y vs O\n")
cat("========================================\n\n")

#===============================================================================
# SECTION 1: PREPARE DATA
#===============================================================================

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

# Subset to Y and O only
yo_samples <- sample_meta$SampleID[sample_meta$Groups %in% c("Y", "O")]
yo_meta <- sample_meta[sample_meta$SampleID %in% yo_samples, ]
rownames(yo_meta) <- yo_meta$SampleID
yo_meta$Groups <- factor(yo_meta$Groups, levels = c("O", "Y"))  # O as reference
yo_meta$Experiment <- factor(yo_meta$Experiment)

cat("Y vs O sample sizes:\n")
print(table(yo_meta$Groups))

cat("\nBatch structure (Y vs O):\n")
batch_table <- table(yo_meta$Groups, yo_meta$Experiment)
print(batch_table)

# Check if batch correction is needed
n_batches <- sum(colSums(batch_table > 0) > 0)
both_groups_in_batch <- sum(apply(batch_table > 0, 2, all))
cat("\nNumber of batches with Y or O samples:", n_batches, "\n")
cat("Batches containing BOTH Y and O:", both_groups_in_batch, "\n")

#===============================================================================
# SECTION 2: PREPARE COUNTS
#===============================================================================

scale_factor <- 100
pathway_counts <- as.data.frame(lapply(PathwayTable[, yo_samples], function(x) {
  as.integer(round(as.numeric(x) * scale_factor))
}))
rownames(pathway_counts) <- rownames(PathwayTable)

# Filter low-count pathways
min_total <- 100
pathway_counts <- pathway_counts[rowSums(pathway_counts, na.rm = TRUE) >= min_total, ]
pathway_counts[is.na(pathway_counts)] <- 0

cat("\nALDEx2 input:", nrow(pathway_counts), "pathways x", ncol(pathway_counts), "samples\n")

#===============================================================================
# SECTION 3: ALDEx2 GLM
#===============================================================================

# Decide whether to include batch
use_batch <- both_groups_in_batch >= 2  # need at least 2 batches with both groups

if (use_batch) {
  cat("\nUsing model: ~ Experiment + Groups (with batch correction)\n")
  mm_YO <- model.matrix(~ Experiment + Groups, data = yo_meta)
} else {
  cat("\nUsing model: ~ Groups (no batch correction — groups not overlapping across batches)\n")
  mm_YO <- model.matrix(~ Groups, data = yo_meta)
}

cat("Model matrix columns:", paste(colnames(mm_YO), collapse = ", "), "\n")
cat("Model matrix dimensions:", nrow(mm_YO), "x", ncol(mm_YO), "\n\n")

# Generate CLR instances
cat("Running ALDEx2 CLR transformation (128 Monte Carlo instances)...\n")
clr_YO <- aldex.clr(pathway_counts, mm_YO, mc.samples = 128, denom = "all", verbose = TRUE)

# Run GLM
cat("Running ALDEx2 GLM...\n")
glm_YO <- aldex.glm(clr_YO, mm_YO)
glm_YO$Pathway <- rownames(glm_YO)

cat("\nGLM result columns:\n")
print(colnames(glm_YO))

# Find the Groups Y coefficient columns
groups_cols <- grep("GroupsY", colnames(glm_YO), value = TRUE)
cat("\nColumns for Groups Y effect:\n")
print(groups_cols)

# Extract relevant columns
est_col <- "GroupsY:Est"
se_col <- "GroupsY:SE"
tval_col <- "GroupsY:t.val"
pval_col <- "GroupsY:pval"
bh_col <- "GroupsY:pval.padj"

#===============================================================================
# SECTION 4: RESULTS
#===============================================================================

# Create clean results dataframe
aldex_glm_YO <- data.frame(
  Pathway = rownames(glm_YO),
  Estimate_CLR = as.numeric(glm_YO[, est_col]),
  SE = as.numeric(glm_YO[, se_col]),
  t.val = as.numeric(glm_YO[, tval_col]),
  pval = as.numeric(glm_YO[, pval_col]),
  pval_BH = as.numeric(glm_YO[, bh_col]),
  stringsAsFactors = FALSE
)

aldex_glm_YO <- aldex_glm_YO[order(aldex_glm_YO$pval), ]

cat("\nTop 20 pathways (ALDEx2 GLM, Y vs O):\n")
print(head(aldex_glm_YO[, c("Pathway", "Estimate_CLR", "SE", "pval", "pval_BH")], 20))

# B6 pathways
b6_patterns <- c("PYRIDOX", "pyridox", "pyridoxal", "PLP", "PWY0.845", "PWY0-845")
b6_idx <- grepl(paste(b6_patterns, collapse = "|"), aldex_glm_YO$Pathway, ignore.case = TRUE)

cat("\n--- B6 Pathways: ALDEx2 GLM Y vs O ---\n")
if (any(b6_idx)) {
  b6_results <- aldex_glm_YO[b6_idx, ]
  print(b6_results)

  # Save B6-only results
  write.csv(b6_results,
            "response_to_reviewers/ALDEx2_B6_Pathways_Y_vs_O_GLM.csv",
            row.names = FALSE)
  cat("\nSaved: response_to_reviewers/ALDEx2_B6_Pathways_Y_vs_O_GLM.csv\n")
}

# Save full results
write.csv(aldex_glm_YO,
          "ALDEx2_Pathways_Y_vs_O_GLM_Full.csv",
          row.names = FALSE)
cat("Saved: ALDEx2_Pathways_Y_vs_O_GLM_Full.csv\n")

# Save all GLM columns for reference
write.csv(glm_YO,
          "ALDEx2_GLM_Y_vs_O_AllColumns.csv")
cat("Saved: ALDEx2_GLM_Y_vs_O_AllColumns.csv\n")

#===============================================================================
# SECTION 5: COMPARISON WITH CLASSIC ALDEx2
#===============================================================================

cat("\n--- Comparison: Classic ALDEx2 vs GLM for B6 Pathways (Y vs O) ---\n\n")

# Load previously saved classic results if available
classic_file <- "response_to_reviewers_archive/ALDEx2_Pathways_Y_vs_O_Full.csv"
if (file.exists(classic_file)) {
  classic <- read.csv(classic_file)
  b6_classic <- classic[grepl(paste(b6_patterns, collapse = "|"),
                               classic$Pathway, ignore.case = TRUE), ]

  cat("Classic ALDEx2 (Welch t-test):\n")
  print(b6_classic[, c("Pathway", "diff.btw", "effect", "we.ep", "we.eBH")])

  cat("\nALDEx2 GLM:\n")
  if (any(b6_idx)) {
    print(b6_results[, c("Pathway", "Estimate_CLR", "SE", "pval", "pval_BH")])
  }
}

cat("\n\nDone!\n")
sessionInfo()
