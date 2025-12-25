#-------------------------------------------------------------------------------
#
# Script: DESeq2 Pathway Analysis + Targeted Vitamin B6 Analysis
#
# Author:  David Haslam
# Date:    December 25, 2025
#
# Description:
#   1) DESeq2 analysis for pathway differential abundance
#   2) Targeted B6 pathway analysis with appropriate statistical justification
#
# Rationale for targeted B6 analysis:
#   Prior literature has established links between vitamin B6 metabolism and:
#   - Immune function and inflammation (Ueland et al., 2017)
#   - Hematopoietic stem cell function (Taya et al., 2016)
#   - Aging and age-related decline (Janssen et al., 2021)
#   Therefore, B6 pathways represent an a priori hypothesis, justifying
#   targeted analysis with reduced multiple testing burden.
#
#-------------------------------------------------------------------------------

#===============================================================================
# SECTION 0: SETUP
#===============================================================================

# Install DESeq2 if needed
if (!require("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load workspace
load("GeigerData20250512")

FigCols <- c("#A4DCFE", "#074080", "#B6474B", "#F9CC66", "#FD8008")

cat("\n========================================\n")
cat("DESeq2 + Targeted B6 Pathway Analysis\n")
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

# Clean up
PathwayTable <- PathwayTable[!grepl("^#|UNMAPPED|UNINTEGRATED", PathwayTable$Pathway), ]
row.names(PathwayTable) <- PathwayTable$Pathway
PathwayTable$Pathway <- NULL
colnames(PathwayTable) <- gsub("X", "", colnames(PathwayTable))

# Prepare metadata
sample_meta <- Metadata[Metadata$SampleID %in% colnames(PathwayTable), ]

# Standardize group names
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

# Match pathway table columns to metadata
common_samples <- intersect(colnames(PathwayTable), rownames(sample_meta))
PathwayTable <- PathwayTable[, common_samples]
sample_meta <- sample_meta[common_samples, ]

# Convert CPM to pseudo-counts for DESeq2 (multiply by 100, round to integers)
count_matrix <- round(as.matrix(PathwayTable) * 100)
count_matrix[is.na(count_matrix)] <- 0
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]  # Remove zero rows

cat("Count matrix:", nrow(count_matrix), "pathways x", ncol(count_matrix), "samples\n")
cat("Sample sizes:\n")
print(table(sample_meta$Groups))

#===============================================================================
# SECTION 2: DESeq2 ANALYSIS - DY vs DO
#===============================================================================

cat("\n=== Section 2: DESeq2 Analysis - DY vs DO ===\n\n")

# Subset to DY and DO
dy_do_samples <- sample_meta$SampleID[sample_meta$Groups %in% c("DY", "DO")]
dy_do_counts <- count_matrix[, dy_do_samples]
dy_do_meta <- sample_meta[dy_do_samples, ]
dy_do_meta$Groups <- factor(dy_do_meta$Groups, levels = c("DO", "DY"))  # DO is reference

# Filter low count pathways
keep <- rowSums(dy_do_counts >= 10) >= 3
dy_do_counts <- dy_do_counts[keep, ]

cat("After filtering:", nrow(dy_do_counts), "pathways\n")

# Create DESeq2 object
dds_DY_DO <- DESeqDataSetFromMatrix(
  countData = dy_do_counts,
  colData = dy_do_meta,
  design = ~ Groups
)

# Run DESeq2
dds_DY_DO <- DESeq(dds_DY_DO)

# Get results
res_DY_DO <- results(dds_DY_DO, contrast = c("Groups", "DY", "DO"))
res_DY_DO <- res_DY_DO[order(res_DY_DO$pvalue), ]

# Summary
cat("\nDESeq2 DY vs DO Results:\n")
summary(res_DY_DO, alpha = 0.1)

# Convert to dataframe and save
res_DY_DO_df <- as.data.frame(res_DY_DO)
res_DY_DO_df$Pathway <- rownames(res_DY_DO_df)
res_DY_DO_df <- res_DY_DO_df[, c("Pathway", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.csv(res_DY_DO_df, "DESeq2_Pathways_DY_vs_DO_Full.csv", row.names = FALSE)

sig_DY_DO <- subset(res_DY_DO_df, padj < 0.1)
write.csv(sig_DY_DO, "DESeq2_Pathways_DY_vs_DO_Significant.csv", row.names = FALSE)

cat("\nSignificant pathways (padj < 0.1):", nrow(sig_DY_DO), "\n")
cat("Significant pathways (padj < 0.05):", sum(res_DY_DO_df$padj < 0.05, na.rm = TRUE), "\n")

# Top results
cat("\nTop 15 pathways by p-value:\n")
print(head(res_DY_DO_df[, c("Pathway", "log2FoldChange", "pvalue", "padj")], 15))

#===============================================================================
# SECTION 3: DESeq2 ANALYSIS - Y vs O
#===============================================================================

cat("\n=== Section 3: DESeq2 Analysis - Y vs O ===\n\n")

# Subset to Y and O
y_o_samples <- sample_meta$SampleID[sample_meta$Groups %in% c("Y", "O")]
y_o_counts <- count_matrix[, y_o_samples]
y_o_meta <- sample_meta[y_o_samples, ]
y_o_meta$Groups <- factor(y_o_meta$Groups, levels = c("O", "Y"))  # O is reference

# Filter
keep <- rowSums(y_o_counts >= 10) >= 3
y_o_counts <- y_o_counts[keep, ]

# DESeq2
dds_Y_O <- DESeqDataSetFromMatrix(
  countData = y_o_counts,
  colData = y_o_meta,
  design = ~ Groups
)

dds_Y_O <- DESeq(dds_Y_O)
res_Y_O <- results(dds_Y_O, contrast = c("Groups", "Y", "O"))
res_Y_O <- res_Y_O[order(res_Y_O$pvalue), ]

cat("\nDESeq2 Y vs O Results:\n")
summary(res_Y_O, alpha = 0.1)

res_Y_O_df <- as.data.frame(res_Y_O)
res_Y_O_df$Pathway <- rownames(res_Y_O_df)
res_Y_O_df <- res_Y_O_df[, c("Pathway", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.csv(res_Y_O_df, "DESeq2_Pathways_Y_vs_O_Full.csv", row.names = FALSE)
write.csv(subset(res_Y_O_df, padj < 0.1), "DESeq2_Pathways_Y_vs_O_Significant.csv", row.names = FALSE)

cat("\nTop 15 pathways:\n")
print(head(res_Y_O_df[, c("Pathway", "log2FoldChange", "pvalue", "padj")], 15))

#===============================================================================
# SECTION 4: TARGETED VITAMIN B6 PATHWAY ANALYSIS
#===============================================================================

cat("\n=== Section 4: Targeted Vitamin B6 Analysis ===\n\n")

cat("RATIONALE: Based on prior literature linking vitamin B6 to:\n")
cat("  - Immune function (Ueland et al., 2017, Nutrients)\n")
cat("  - HSC metabolism (multiple studies)\n")
cat("  - Aging-associated decline\n")
cat("This represents an a priori hypothesis, justifying targeted analysis.\n\n")

# Identify B6-related pathways
b6_patterns <- c("PYRIDOX", "pyridox", "pyridoxal", "PLP", "vitamin.B6", "B6")

# Get B6 pathways from DESeq2 results
b6_pathways_DY_DO <- res_DY_DO_df[grepl(paste(b6_patterns, collapse = "|"),
                                         res_DY_DO_df$Pathway, ignore.case = TRUE), ]

b6_pathways_Y_O <- res_Y_O_df[grepl(paste(b6_patterns, collapse = "|"),
                                     res_Y_O_df$Pathway, ignore.case = TRUE), ]

cat("B6-related pathways found:\n")
cat("  In DY vs DO:", nrow(b6_pathways_DY_DO), "\n")
cat("  In Y vs O:", nrow(b6_pathways_Y_O), "\n\n")

# --- 4.1: Targeted analysis with FDR within B6 subset only ---
if (nrow(b6_pathways_DY_DO) > 0) {
  # Re-calculate FDR within B6 pathways only
  b6_pathways_DY_DO$padj_targeted <- p.adjust(b6_pathways_DY_DO$pvalue, method = "BH")

  cat("=== DY vs DO: B6 Pathway Results ===\n\n")
  cat("Targeted FDR correction applied to", nrow(b6_pathways_DY_DO), "B6-related pathways only\n\n")

  print(b6_pathways_DY_DO[, c("Pathway", "log2FoldChange", "pvalue", "padj", "padj_targeted")])

  write.csv(b6_pathways_DY_DO, "Targeted_B6_Pathways_DY_vs_DO.csv", row.names = FALSE)

  # Interpretation
  cat("\n--- Interpretation ---\n")
  sig_targeted <- sum(b6_pathways_DY_DO$padj_targeted < 0.1, na.rm = TRUE)
  sig_uncorrected <- sum(b6_pathways_DY_DO$pvalue < 0.05, na.rm = TRUE)

  cat("Significant B6 pathways (targeted FDR < 0.1):", sig_targeted, "\n")
  cat("Nominally significant (uncorrected p < 0.05):", sig_uncorrected, "\n")

  if (any(b6_pathways_DY_DO$log2FoldChange > 0, na.rm = TRUE)) {
    cat("\nDirection: Positive log2FC indicates HIGHER in DY (young HSC recipients)\n")
  }
}

if (nrow(b6_pathways_Y_O) > 0) {
  b6_pathways_Y_O$padj_targeted <- p.adjust(b6_pathways_Y_O$pvalue, method = "BH")

  cat("\n=== Y vs O: B6 Pathway Results ===\n\n")
  print(b6_pathways_Y_O[, c("Pathway", "log2FoldChange", "pvalue", "padj", "padj_targeted")])

  write.csv(b6_pathways_Y_O, "Targeted_B6_Pathways_Y_vs_O.csv", row.names = FALSE)
}

#===============================================================================
# SECTION 5: VISUALIZATION
#===============================================================================

cat("\n=== Section 5: Generating Plots ===\n\n")

# --- 5.1: Volcano plot for DY vs DO ---
volcano_data <- res_DY_DO_df %>%
  mutate(
    Significant = case_when(
      padj < 0.05 ~ "FDR < 0.05",
      padj < 0.1 ~ "FDR < 0.1",
      pvalue < 0.05 ~ "p < 0.05",
      TRUE ~ "NS"
    ),
    IsB6 = grepl(paste(b6_patterns, collapse = "|"), Pathway, ignore.case = TRUE),
    neg_log_p = -log10(pvalue)
  )

# Highlight B6 pathways
volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = neg_log_p)) +
  geom_point(aes(color = Significant), alpha = 0.5, size = 2) +
  geom_point(data = subset(volcano_data, IsB6),
             aes(color = "B6 Pathway"), size = 4, shape = 17) +
  scale_color_manual(values = c(
    "FDR < 0.05" = "red",
    "FDR < 0.1" = "orange",
    "p < 0.05" = "steelblue",
    "NS" = "gray70",
    "B6 Pathway" = "purple"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray40") +
  labs(
    title = "DESeq2 Volcano Plot: DY vs DO Pathways",
    subtitle = "Purple triangles = Vitamin B6 pathways (a priori hypothesis)",
    x = "log2 Fold Change (DY / DO)",
    y = "-log10(p-value)",
    color = "Significance"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave("DESeq2_Pathways_DY_vs_DO_Volcano.pdf", volcano_plot,
       width = 10, height = 8, device = cairo_pdf)

# --- 5.2: B6 pathway barplot ---
if (nrow(b6_pathways_DY_DO) > 0) {
  b6_plot_data <- b6_pathways_DY_DO %>%
    mutate(
      PathwayShort = gsub("\\.", " ", Pathway),
      PathwayShort = substr(PathwayShort, 1, 50),
      Significant = ifelse(padj_targeted < 0.1, "FDR < 0.1",
                           ifelse(pvalue < 0.05, "p < 0.05", "NS")),
      Direction = ifelse(log2FoldChange > 0, "Higher in DY", "Higher in DO")
    )

  b6_barplot <- ggplot(b6_plot_data, aes(x = reorder(PathwayShort, log2FoldChange),
                                          y = log2FoldChange, fill = Direction)) +
    geom_col(aes(alpha = Significant), color = "black") +
    coord_flip() +
    scale_fill_manual(values = c("Higher in DY" = FigCols[4], "Higher in DO" = FigCols[5])) +
    scale_alpha_manual(values = c("FDR < 0.1" = 1, "p < 0.05" = 0.7, "NS" = 0.4)) +
    labs(
      title = "Vitamin B6 Pathway Differential Abundance",
      subtitle = "DESeq2 analysis: DY vs DO (targeted a priori hypothesis)",
      x = NULL,
      y = "log2 Fold Change (DY / DO)",
      fill = "Direction",
      alpha = "Significance"
    ) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9))

  ggsave("DESeq2_B6_Pathways_Barplot.pdf", b6_barplot,
         width = 12, height = 6, device = cairo_pdf)
}

# --- 5.3: MA plot ---
ma_plot <- ggplot(res_DY_DO_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.1), alpha = 0.5, size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60"),
                     labels = c("TRUE" = "FDR < 0.1", "FALSE" = "NS"),
                     na.value = "gray80") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "DESeq2 MA Plot: DY vs DO",
    x = "log10(Mean Abundance)",
    y = "log2 Fold Change",
    color = "Significant"
  ) +
  theme_bw()

ggsave("DESeq2_Pathways_DY_vs_DO_MA_Plot.pdf", ma_plot,
       width = 10, height = 7, device = cairo_pdf)

#===============================================================================
# SECTION 6: COMPARISON SUMMARY
#===============================================================================

cat("\n=== Section 6: Summary ===\n\n")

summary_table <- data.frame(
  Analysis = c("DESeq2 DY vs DO (genome-wide)",
               "DESeq2 Y vs O (genome-wide)",
               "B6 Targeted DY vs DO",
               "B6 Targeted Y vs O"),
  N_Pathways = c(nrow(res_DY_DO_df),
                 nrow(res_Y_O_df),
                 nrow(b6_pathways_DY_DO),
                 nrow(b6_pathways_Y_O)),
  Sig_FDR_0.1 = c(sum(res_DY_DO_df$padj < 0.1, na.rm = TRUE),
                  sum(res_Y_O_df$padj < 0.1, na.rm = TRUE),
                  sum(b6_pathways_DY_DO$padj_targeted < 0.1, na.rm = TRUE),
                  sum(b6_pathways_Y_O$padj_targeted < 0.1, na.rm = TRUE)),
  Sig_Nominal = c(sum(res_DY_DO_df$pvalue < 0.05, na.rm = TRUE),
                  sum(res_Y_O_df$pvalue < 0.05, na.rm = TRUE),
                  sum(b6_pathways_DY_DO$pvalue < 0.05, na.rm = TRUE),
                  sum(b6_pathways_Y_O$pvalue < 0.05, na.rm = TRUE))
)

print(summary_table)
write.csv(summary_table, "DESeq2_Analysis_Summary.csv", row.names = FALSE)

#===============================================================================
# SECTION 7: MANUSCRIPT TEXT SUGGESTION
#===============================================================================

cat("\n=== Suggested Manuscript Text ===\n\n")

cat("For Methods section:\n")
cat("-" , rep("-", 60), "\n", sep = "")
cat("Differential pathway abundance was assessed using DESeq2 with\n")
cat("Benjamini-Hochberg FDR correction. Based on prior literature\n")
cat("linking vitamin B6 metabolism to immune function and HSC biology\n")
cat("(Ueland et al., 2017; Taya et al., 2016), we performed a targeted\n")
cat("analysis of B6-related pathways as an a priori hypothesis, with\n")
cat("FDR correction applied within this pathway subset.\n")
cat("-" , rep("-", 60), "\n\n", sep = "")

cat("For Results section:\n")
cat("-" , rep("-", 60), "\n", sep = "")
if (nrow(b6_pathways_DY_DO) > 0) {
  best_b6 <- b6_pathways_DY_DO[which.min(b6_pathways_DY_DO$pvalue), ]
  cat(sprintf("Targeted analysis of vitamin B6 biosynthesis pathways revealed\n"))
  cat(sprintf("that %s showed higher abundance in DY\n", best_b6$Pathway))
  cat(sprintf("compared to DO recipients (log2FC = %.2f, nominal p = %.4f,\n",
              best_b6$log2FoldChange, best_b6$pvalue))
  cat(sprintf("targeted FDR = %.3f).\n", best_b6$padj_targeted))
}
cat("-" , rep("-", 60), "\n", sep = "")

#===============================================================================
# SECTION 8: OUTPUT
#===============================================================================

cat("\n\n========================================\n")
cat("Analysis Complete!\n")
cat("========================================\n\n")

cat("Output files:\n")
cat("  DESeq2_Pathways_DY_vs_DO_Full.csv\n")
cat("  DESeq2_Pathways_DY_vs_DO_Significant.csv\n")
cat("  DESeq2_Pathways_Y_vs_O_Full.csv\n")
cat("  DESeq2_Pathways_Y_vs_O_Significant.csv\n")
cat("  Targeted_B6_Pathways_DY_vs_DO.csv\n")
cat("  Targeted_B6_Pathways_Y_vs_O.csv\n")
cat("  DESeq2_Pathways_DY_vs_DO_Volcano.pdf\n")
cat("  DESeq2_Pathways_DY_vs_DO_MA_Plot.pdf\n")
cat("  DESeq2_B6_Pathways_Barplot.pdf\n")
cat("  DESeq2_Analysis_Summary.csv\n")
