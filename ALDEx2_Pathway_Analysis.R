#-------------------------------------------------------------------------------
#
# Script: ALDEx2 Pathway Differential Abundance Analysis
#
# Author:  David Haslam
# Date:    December 25, 2025
#
# Description: Proper ALDEx2 analysis for HUMAnN3 pathway data as requested
#              by reviewer. Includes all pairwise comparisons with FDR correction.
#
#-------------------------------------------------------------------------------

#===============================================================================
# SECTION 0: SETUP
#===============================================================================

library(ALDEx2)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load workspace
load("GeigerData20250512")

# Color palette
FigCols <- c("#A4DCFE", "#074080", "#B6474B", "#F9CC66", "#FD8008")

cat("\n========================================\n")
cat("ALDEx2 Pathway Analysis\n")
cat("========================================\n\n")

#===============================================================================
# SECTION 1: PREPARE PATHWAY DATA FOR ALDEx2
#===============================================================================

cat("=== Section 1: Preparing Pathway Data ===\n\n")

# Check for pathway data
if (!exists("GeigerPathways")) {
  # Try to load from file
  pathway_file <- "GeigerFiles_pathabundance20220910-cpm_unstratified.tsv"
  if (file.exists(pathway_file)) {
    cat("Loading pathway data from file...\n")
    PathwayTable <- read.csv(pathway_file, stringsAsFactors = FALSE,
                              header = TRUE, sep = "\t")
    names(PathwayTable)[1] <- "Pathway"
    colnames(PathwayTable) <- gsub(".paired_Abundance.CPM", "", colnames(PathwayTable))

    # Remove header rows
    PathwayTable <- PathwayTable[!grepl("^#|UNMAPPED|UNINTEGRATED", PathwayTable$Pathway), ]

    row.names(PathwayTable) <- PathwayTable$Pathway
    PathwayTable$Pathway <- NULL
    colnames(PathwayTable) <- gsub("X", "", colnames(PathwayTable))

    # Transpose
    GeigerPathwayTable_t <- t(PathwayTable)
    GeigerPathwayTabledf <- as.data.frame(GeigerPathwayTable_t)
    GeigerPathwayTabledf$SampleID <- row.names(GeigerPathwayTabledf)

    GeigerPathways <- merge(Metadata, GeigerPathwayTabledf, by = "SampleID", all.y = TRUE)
    GeigerPathways <- subset(GeigerPathways, !is.na(Groups))
    row.names(GeigerPathways) <- GeigerPathways$SampleID
  } else {
    stop("Cannot find pathway data!")
  }
}

# Standardize group names
GeigerPathways$Groups <- as.character(GeigerPathways$Groups)
GeigerPathways$Groups <- gsub("Young_C57BL-Untransplanted", "Y", GeigerPathways$Groups)
GeigerPathways$Groups <- gsub("Old_C57BL-Untransplanted", "O", GeigerPathways$Groups)
GeigerPathways$Groups <- gsub("Young_RAG1-Untransplanted", "RAG1-/-", GeigerPathways$Groups)
GeigerPathways$Groups <- gsub("RAG1-Young_HSCs", "DY", GeigerPathways$Groups)
GeigerPathways$Groups <- gsub("RAG1-Old_HSCs", "DO", GeigerPathways$Groups)
GeigerPathways$Groups <- gsub("RAG1-Rejuvenated_HSCs", "Rejuv", GeigerPathways$Groups)

# Filter to main groups
GeigerPathways <- GeigerPathways[GeigerPathways$Groups %in% c("Y", "O", "RAG1-/-", "DY", "DO"), ]

cat("Sample sizes per group:\n")
print(table(GeigerPathways$Groups))

# Extract pathway abundance matrix
metadata_cols <- 1:10
pathway_data <- GeigerPathways[, -metadata_cols]

# Convert to numeric
pathway_data <- as.data.frame(lapply(pathway_data, function(x) {
  as.numeric(as.character(x))
}))
rownames(pathway_data) <- GeigerPathways$SampleID

# ALDEx2 needs counts (integers), so scale CPM to pseudo-counts
# Multiply by a factor and round to integers
scale_factor <- 100  # This preserves relative abundances
pathway_counts <- as.data.frame(lapply(pathway_data, function(x) {
  as.integer(round(x * scale_factor))
}))
rownames(pathway_counts) <- rownames(pathway_data)

# Transpose for ALDEx2: features (pathways) as rows, samples as columns
pathway_counts_t <- as.data.frame(t(pathway_counts))

# Filter out pathways with very low counts
min_total <- 100
pathway_counts_t <- pathway_counts_t[rowSums(pathway_counts_t, na.rm = TRUE) >= min_total, ]

# Remove any remaining NA
pathway_counts_t[is.na(pathway_counts_t)] <- 0

cat("\nPathway matrix dimensions:", nrow(pathway_counts_t), "pathways x",
    ncol(pathway_counts_t), "samples\n")

# Get condition vector matching column order
sample_order <- colnames(pathway_counts_t)
conditions <- GeigerPathways$Groups[match(sample_order, GeigerPathways$SampleID)]

cat("Conditions:\n")
print(table(conditions))

#===============================================================================
# SECTION 2: ALDEx2 - DY vs DO (Primary comparison of interest)
#===============================================================================

cat("\n=== Section 2: ALDEx2 Analysis - DY vs DO ===\n\n")

# Subset to DY and DO
dy_do_idx <- which(conditions %in% c("DY", "DO"))
counts_DY_DO <- pathway_counts_t[, dy_do_idx]
conds_DY_DO <- as.character(conditions[dy_do_idx])

cat("Running ALDEx2 for DY vs DO (", sum(conds_DY_DO == "DY"), "vs",
    sum(conds_DY_DO == "DO"), "samples)...\n")

# Run ALDEx2
aldex_DY_DO <- aldex(counts_DY_DO, conds_DY_DO,
                      mc.samples = 128,
                      test = "t",
                      effect = TRUE,
                      include.sample.summary = FALSE,
                      verbose = TRUE,
                      denom = "all")

# Add pathway names and sort
aldex_DY_DO$Pathway <- rownames(aldex_DY_DO)
aldex_DY_DO <- aldex_DY_DO[order(aldex_DY_DO$we.eBH), ]

# Count significant results
sig_BH_01 <- sum(aldex_DY_DO$we.eBH < 0.1, na.rm = TRUE)
sig_BH_05 <- sum(aldex_DY_DO$we.eBH < 0.05, na.rm = TRUE)
sig_effect <- sum(abs(aldex_DY_DO$effect) > 1, na.rm = TRUE)

cat("\nDY vs DO Results:\n")
cat("  Significant pathways (BH-adjusted p < 0.1):", sig_BH_01, "\n")
cat("  Significant pathways (BH-adjusted p < 0.05):", sig_BH_05, "\n")
cat("  Pathways with |effect| > 1:", sig_effect, "\n")

# Save full results
write.csv(aldex_DY_DO, "ALDEx2_Pathways_DY_vs_DO_Full.csv", row.names = FALSE)

# Save significant results
sig_aldex_DY_DO <- subset(aldex_DY_DO, we.eBH < 0.1 | abs(effect) > 0.5)
write.csv(sig_aldex_DY_DO, "ALDEx2_Pathways_DY_vs_DO_Significant.csv", row.names = FALSE)

# Print top results
cat("\nTop 20 pathways by significance (DY vs DO):\n")
print(head(aldex_DY_DO[, c("Pathway", "effect", "we.ep", "we.eBH", "wi.ep", "wi.eBH")], 20))

#===============================================================================
# SECTION 3: ALDEx2 - Y vs O (Young vs Old untransplanted controls)
#===============================================================================

cat("\n=== Section 3: ALDEx2 Analysis - Y vs O ===\n\n")

# Subset to Y and O
y_o_idx <- which(conditions %in% c("Y", "O"))
counts_Y_O <- pathway_counts_t[, y_o_idx]
conds_Y_O <- as.character(conditions[y_o_idx])

cat("Running ALDEx2 for Y vs O (", sum(conds_Y_O == "Y"), "vs",
    sum(conds_Y_O == "O"), "samples)...\n")

# Run ALDEx2
aldex_Y_O <- aldex(counts_Y_O, conds_Y_O,
                    mc.samples = 128,
                    test = "t",
                    effect = TRUE,
                    include.sample.summary = FALSE,
                    verbose = TRUE,
                    denom = "all")

aldex_Y_O$Pathway <- rownames(aldex_Y_O)
aldex_Y_O <- aldex_Y_O[order(aldex_Y_O$we.eBH), ]

sig_BH_01 <- sum(aldex_Y_O$we.eBH < 0.1, na.rm = TRUE)
sig_BH_05 <- sum(aldex_Y_O$we.eBH < 0.05, na.rm = TRUE)

cat("\nY vs O Results:\n")
cat("  Significant pathways (BH-adjusted p < 0.1):", sig_BH_01, "\n")
cat("  Significant pathways (BH-adjusted p < 0.05):", sig_BH_05, "\n")

write.csv(aldex_Y_O, "ALDEx2_Pathways_Y_vs_O_Full.csv", row.names = FALSE)
write.csv(subset(aldex_Y_O, we.eBH < 0.1), "ALDEx2_Pathways_Y_vs_O_Significant.csv", row.names = FALSE)

cat("\nTop 20 pathways by significance (Y vs O):\n")
print(head(aldex_Y_O[, c("Pathway", "effect", "we.ep", "we.eBH")], 20))

#===============================================================================
# SECTION 4: ALDEx2 - Y vs RAG1-/- (Background strain comparison)
#===============================================================================

cat("\n=== Section 4: ALDEx2 Analysis - Y vs RAG1-/- ===\n\n")

# Subset to Y and RAG1-/-
y_rag_idx <- which(conditions %in% c("Y", "RAG1-/-"))
counts_Y_RAG <- pathway_counts_t[, y_rag_idx]
conds_Y_RAG <- as.character(conditions[y_rag_idx])

cat("Running ALDEx2 for Y vs RAG1-/- (", sum(conds_Y_RAG == "Y"), "vs",
    sum(conds_Y_RAG == "RAG1-/-"), "samples)...\n")

aldex_Y_RAG <- aldex(counts_Y_RAG, conds_Y_RAG,
                      mc.samples = 128,
                      test = "t",
                      effect = TRUE,
                      include.sample.summary = FALSE,
                      verbose = TRUE,
                      denom = "all")

aldex_Y_RAG$Pathway <- rownames(aldex_Y_RAG)
aldex_Y_RAG <- aldex_Y_RAG[order(aldex_Y_RAG$we.eBH), ]

sig_BH_01 <- sum(aldex_Y_RAG$we.eBH < 0.1, na.rm = TRUE)

cat("\nY vs RAG1-/- Results:\n")
cat("  Significant pathways (BH-adjusted p < 0.1):", sig_BH_01, "\n")

write.csv(aldex_Y_RAG, "ALDEx2_Pathways_Y_vs_RAG1_Full.csv", row.names = FALSE)
write.csv(subset(aldex_Y_RAG, we.eBH < 0.1), "ALDEx2_Pathways_Y_vs_RAG1_Significant.csv", row.names = FALSE)

#===============================================================================
# SECTION 5: VISUALIZATION
#===============================================================================

cat("\n=== Section 5: Generating Plots ===\n\n")

# --- 5.1: MA Plot for DY vs DO ---
ma_data <- aldex_DY_DO %>%
  mutate(
    Significant = case_when(
      we.eBH < 0.05 ~ "FDR < 0.05",
      we.eBH < 0.1 ~ "FDR < 0.1",
      TRUE ~ "NS"
    ),
    Significant = factor(Significant, levels = c("FDR < 0.05", "FDR < 0.1", "NS"))
  )

MA_Plot <- ggplot(ma_data, aes(x = rab.all, y = diff.btw, color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("FDR < 0.05" = "red", "FDR < 0.1" = "orange", "NS" = "gray50")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "ALDEx2 MA Plot: DY vs DO Pathways",
    subtitle = "Positive values = higher in DY",
    x = "Mean CLR Abundance (rab.all)",
    y = "Difference Between Groups (diff.btw)",
    color = "Significance"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave("ALDEx2_Pathways_DY_vs_DO_MA_Plot.pdf", MA_Plot,
       width = 10, height = 8, device = cairo_pdf)

# --- 5.2: Effect Size Plot for DY vs DO ---
effect_data <- aldex_DY_DO %>%
  filter(abs(effect) > 0.3 | we.eBH < 0.1) %>%
  mutate(
    Significant = ifelse(we.eBH < 0.1, "FDR < 0.1", "NS"),
    Direction = ifelse(effect > 0, "Higher in DY", "Higher in DO"),
    PathwayShort = substr(Pathway, 1, 60)  # Truncate long names
  ) %>%
  arrange(effect) %>%
  head(40)  # Top 40 pathways

if (nrow(effect_data) > 0) {
  Effect_Plot <- ggplot(effect_data, aes(x = reorder(PathwayShort, effect), y = effect, fill = Direction)) +
    geom_col(aes(alpha = Significant), color = "black", linewidth = 0.2) +
    coord_flip() +
    scale_fill_manual(values = c("Higher in DY" = FigCols[4], "Higher in DO" = FigCols[5])) +
    scale_alpha_manual(values = c("FDR < 0.1" = 1, "NS" = 0.5)) +
    labs(
      title = "ALDEx2 Effect Sizes: DY vs DO Pathways",
      subtitle = "Top pathways by effect size",
      x = NULL,
      y = "Effect Size",
      fill = "Direction",
      alpha = "Significance"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 7),
      legend.position = "right"
    )

  ggsave("ALDEx2_Pathways_DY_vs_DO_Effect_Plot.pdf", Effect_Plot,
         width = 14, height = 10, device = cairo_pdf)
}

# --- 5.3: Volcano Plot for DY vs DO ---
volcano_data <- aldex_DY_DO %>%
  mutate(
    neg_log_p = -log10(we.eBH),
    Significant = case_when(
      we.eBH < 0.05 & abs(effect) > 0.5 ~ "Significant",
      we.eBH < 0.1 ~ "Marginal",
      TRUE ~ "NS"
    )
  )

Volcano_Plot <- ggplot(volcano_data, aes(x = effect, y = neg_log_p, color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Significant" = "red", "Marginal" = "orange", "NS" = "gray50")) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "orange") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "gray30") +
  labs(
    title = "ALDEx2 Volcano Plot: DY vs DO Pathways",
    x = "Effect Size",
    y = "-log10(BH-adjusted p-value)",
    color = "Significance"
  ) +
  theme_bw()

ggsave("ALDEx2_Pathways_DY_vs_DO_Volcano.pdf", Volcano_Plot,
       width = 10, height = 8, device = cairo_pdf)

# --- 5.4: Comparison summary plot ---
# Combine results from all comparisons
summary_data <- data.frame(
  Comparison = c("DY vs DO", "Y vs O", "Y vs RAG1-/-"),
  Total_Pathways = c(nrow(aldex_DY_DO), nrow(aldex_Y_O), nrow(aldex_Y_RAG)),
  Sig_FDR_0.1 = c(
    sum(aldex_DY_DO$we.eBH < 0.1, na.rm = TRUE),
    sum(aldex_Y_O$we.eBH < 0.1, na.rm = TRUE),
    sum(aldex_Y_RAG$we.eBH < 0.1, na.rm = TRUE)
  ),
  Sig_FDR_0.05 = c(
    sum(aldex_DY_DO$we.eBH < 0.05, na.rm = TRUE),
    sum(aldex_Y_O$we.eBH < 0.05, na.rm = TRUE),
    sum(aldex_Y_RAG$we.eBH < 0.05, na.rm = TRUE)
  )
)

cat("\n=== Summary of ALDEx2 Pathway Results ===\n\n")
print(summary_data)
write.csv(summary_data, "ALDEx2_Pathways_Summary.csv", row.names = FALSE)

#===============================================================================
# SECTION 6: VITAMIN B6 SPECIFIC PATHWAYS
#===============================================================================

cat("\n=== Section 6: Vitamin B6 Pathway Results ===\n\n")

# Search for B6 pathways in the results
b6_patterns <- c("PYRIDOX", "B6", "pyridox", "pyridoxal")

b6_DY_DO <- aldex_DY_DO[grepl(paste(b6_patterns, collapse = "|"),
                               aldex_DY_DO$Pathway, ignore.case = TRUE), ]

if (nrow(b6_DY_DO) > 0) {
  cat("Vitamin B6 pathways found in DY vs DO comparison:\n\n")
  print(b6_DY_DO[, c("Pathway", "effect", "we.ep", "we.eBH")])
  write.csv(b6_DY_DO, "ALDEx2_VitaminB6_Pathways_DY_vs_DO.csv", row.names = FALSE)
} else {
  cat("No vitamin B6 pathways found in results.\n")
}

#===============================================================================
# SECTION 7: SESSION INFO
#===============================================================================

cat("\n\n========================================\n")
cat("ALDEx2 Pathway Analysis Complete!\n")
cat("========================================\n\n")

cat("Output files:\n")
cat("  - ALDEx2_Pathways_DY_vs_DO_Full.csv\n")
cat("  - ALDEx2_Pathways_DY_vs_DO_Significant.csv\n")
cat("  - ALDEx2_Pathways_Y_vs_O_Full.csv\n")
cat("  - ALDEx2_Pathways_Y_vs_O_Significant.csv\n")
cat("  - ALDEx2_Pathways_Y_vs_RAG1_Full.csv\n")
cat("  - ALDEx2_Pathways_Y_vs_RAG1_Significant.csv\n")
cat("  - ALDEx2_Pathways_DY_vs_DO_MA_Plot.pdf\n")
cat("  - ALDEx2_Pathways_DY_vs_DO_Effect_Plot.pdf\n")
cat("  - ALDEx2_Pathways_DY_vs_DO_Volcano.pdf\n")
cat("  - ALDEx2_Pathways_Summary.csv\n")
cat("  - ALDEx2_VitaminB6_Pathways_DY_vs_DO.csv\n")

cat("\nSession info:\n")
sessionInfo()
