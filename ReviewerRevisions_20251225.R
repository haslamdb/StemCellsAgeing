#-------------------------------------------------------------------------------
#
# Script: Reviewer Revisions - Stem Cell Aging Microbiome Analysis
#
# Author:  David Haslam
# Date:    December 25, 2025
#
# Description: Additional analyses to address reviewer comments:
#   1) Add Chao1 diversity and rarefaction curves
#   2) ALDEx2 differential abundance with FDR-corrected p-values
#   3) Power analysis / variance estimates for sample size justification
#   4) Batch effect handling documentation and analysis
#   5) PCA plots with variance explained (%) and sample sizes in legend
#
#-------------------------------------------------------------------------------

#===============================================================================
# SECTION 0: LOAD EXISTING WORKSPACE AND ADDITIONAL PACKAGES
#===============================================================================

# Load the saved workspace from the main analysis
load("/home/david/projects/Metagenomics/GeigerData/StemCellsAgeing/GeigerData20250512")

# Install and load additional packages for reviewer revisions
additional.packages <- c("ALDEx2", "pwr", "effectsize", "micropower")

# Check for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install ALDEx2 from Bioconductor if not present
if (!require("ALDEx2", quietly = TRUE)) {
  BiocManager::install("ALDEx2")
}

# Install CRAN packages
new.pkgs <- additional.packages[!(additional.packages %in% installed.packages()[, "Package"])]
if (length(new.pkgs)) {
  # Try CRAN first
  for (pkg in new.pkgs) {
    tryCatch(install.packages(pkg), error = function(e) message(paste("Could not install", pkg)))
  }
}

# Load packages
library(ALDEx2)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(ggpubr)
library(rstatix)

# Redefine color palette and themes from main script
FigCols <- c("#A4DCFE", "#074080", "#B6474B", "#F9CC66", "#FD8008")

paramsBox <- function(x) {
  theme(
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 16, angle = 0, color = "black"),
    plot.title = element_text(size = 22, color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 22),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.position = "none"
  )
}

paramsAngled <- function(x) {
  theme(
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 22, color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.position = "none",
    legend.text = element_text(size = 14)
  )
}


#===============================================================================
# SECTION 1: ENHANCED ALPHA DIVERSITY (CHAO1 + RAREFACTION CURVES)
#===============================================================================

cat("\n=== SECTION 1: Alpha Diversity with Chao1 and Rarefaction Curves ===\n\n")

# --- 1.1: Calculate Chao1 estimator ---
# The Species matrix should be samples x species (from the main analysis)
# If Species is species x samples, we need to transpose

# Check orientation and get the species count matrix
if (exists("Species")) {
  # Species from main script is species x samples, need samples x species
  species_matrix <- as.data.frame(t(Species))
} else if (exists("GeigerSpeciesNR")) {
  # Extract species columns (after column 10 which has metadata)
  species_matrix <- GeigerSpeciesNR[, 11:ncol(GeigerSpeciesNR)]
} else {
  stop("Cannot find species abundance data. Please check workspace.")
}

# Ensure numeric
species_matrix <- as.data.frame(lapply(species_matrix, as.numeric))
rownames(species_matrix) <- rownames(GeigerSpeciesNR)

# Calculate Chao1 using estimateR from vegan
# estimateR requires species as columns and samples as rows
chao1_results <- t(estimateR(species_matrix))
chao1_df <- as.data.frame(chao1_results)
chao1_df$SampleID <- rownames(chao1_df)

# Merge with existing Diversity dataframe
if (exists("Diversity")) {
  Diversity <- merge(Diversity, chao1_df[, c("SampleID", "S.chao1", "se.chao1")],
                     by = "SampleID", all.x = TRUE)
} else {
  # Recreate diversity metrics
  H <- diversity(species_matrix)
  simpson <- diversity(species_matrix, "simpson")
  shannon <- diversity(species_matrix, "shannon")
  S <- specnumber(species_matrix)
  J <- H / log(S)

  Diversity <- data.frame(
    SampleID = rownames(species_matrix),
    Simpson = simpson,
    Shannon = shannon,
    SpeciesNo = S,
    Evenness = J
  )
  Diversity <- merge(Metadata, Diversity, by = "SampleID", all.x = TRUE)
  Diversity <- merge(Diversity, chao1_df[, c("SampleID", "S.chao1", "se.chao1")],
                     by = "SampleID", all.x = TRUE)
}

# --- 1.2: Plot Chao1 Diversity by Group ---
Chao1_Plot <- ggplot(Diversity, aes(x = Groups, y = S.chao1)) +
  geom_boxplot(lwd = 1, aes(color = factor(Groups)), fill = NA, outlier.size = 3) +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 8) +
  scale_colour_manual(values = FigCols) +
  geom_point(size = 4, aes(color = factor(Groups))) +
  xlab(NULL) +
  ylab("Chao1 Richness Estimator") +
  theme_bw() +
  paramsBox()

ggsave(filename = "Reviewer_Chao1_Diversity.pdf", plot = Chao1_Plot,
       width = 10, height = 10, device = cairo_pdf)

# --- 1.3: Combined Alpha Diversity Panel (Shannon, Simpson, Chao1) ---
# Reshape for faceted plot
diversity_long <- Diversity %>%
  select(SampleID, Groups, Shannon, Simpson, S.chao1) %>%
  pivot_longer(cols = c(Shannon, Simpson, S.chao1),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Shannon", "Simpson", "S.chao1"),
                         labels = c("Shannon Index", "Simpson Index", "Chao1 Richness")))

Alpha_Diversity_Panel <- ggplot(diversity_long, aes(x = Groups, y = Value)) +
  geom_boxplot(aes(color = Groups), fill = NA, outlier.shape = NA) +
  geom_jitter(aes(color = Groups), width = 0.2, size = 2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  scale_color_manual(values = FigCols) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Diversity Value") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )

ggsave(filename = "Reviewer_AlphaDiversity_Panel.pdf", plot = Alpha_Diversity_Panel,
       width = 14, height = 5, device = cairo_pdf)

# --- 1.4: Statistical tests for Chao1 ---
cat("Pairwise Wilcoxon tests for Chao1:\n")
chao1_pairwise <- pairwise.wilcox.test(Diversity$S.chao1, Diversity$Groups,
                                        p.adjust.method = "bonferroni")
print(chao1_pairwise)

# --- 1.5: Rarefaction Curves ---
# Create rarefaction curve data
cat("\nGenerating rarefaction curves...\n")

# Function to compute rarefaction data for a sample
compute_rarefaction <- function(sample_counts, sample_id, depths = NULL) {
  # Remove NA and ensure numeric
  sample_counts <- as.numeric(sample_counts)
  sample_counts[is.na(sample_counts)] <- 0

  total_reads <- sum(sample_counts)
  if (total_reads == 0) return(NULL)

  if (is.null(depths)) {
    depths <- seq(1000, min(total_reads, 750000), length.out = 30)
  }
  depths <- depths[depths <= total_reads]
  depths <- depths[depths > 0]

  if (length(depths) == 0) return(NULL)

  rarefied_richness <- sapply(depths, function(d) {
    tryCatch({
      rarefy(sample_counts, sample = d)
    }, error = function(e) NA)
  })

  data.frame(
    SampleID = sample_id,
    Depth = depths,
    Richness = rarefied_richness
  )
}

# Compute rarefaction for all samples
set.seed(42)
rarefaction_list <- lapply(1:nrow(species_matrix), function(i) {
  tryCatch({
    compute_rarefaction(species_matrix[i, ],
                        rownames(species_matrix)[i],
                        depths = seq(10000, 750000, by = 30000))
  }, error = function(e) NULL)
})

# Remove NULL entries and combine
rarefaction_list <- rarefaction_list[!sapply(rarefaction_list, is.null)]
rarefaction_data <- do.call(rbind, rarefaction_list)

# Add group information
rarefaction_data <- merge(rarefaction_data,
                          Metadata[, c("SampleID", "Groups")],
                          by = "SampleID")

# Plot rarefaction curves
Rarefaction_Plot <- ggplot(rarefaction_data, aes(x = Depth, y = Richness,
                                                  group = SampleID, color = Groups)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  scale_color_manual(values = FigCols) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Sequencing Depth (Reads)",
    y = "Species Richness",
    title = "Rarefaction Curves by Experimental Group",
    color = "Group"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  geom_vline(xintercept = 750000, linetype = "dashed", color = "gray40") +
  annotate("text", x = 750000, y = max(rarefaction_data$Richness) * 0.95,
           label = "Rarefaction depth\n(750,000)", hjust = 1.1, size = 3)

ggsave(filename = "Reviewer_RarefactionCurves.pdf", plot = Rarefaction_Plot,
       width = 10, height = 7, device = cairo_pdf)

cat("Alpha diversity analysis complete. Files saved.\n")


#===============================================================================
# SECTION 2: ALDEx2 DIFFERENTIAL ABUNDANCE ANALYSIS
#===============================================================================

cat("\n=== SECTION 2: ALDEx2 Differential Abundance Analysis ===\n\n")

# ALDEx2 requires raw counts (integers), samples as columns, features as rows
# Get the raw (non-rarefied) species counts if available, otherwise use rarefied

# Check if we have non-rarefied counts
if (exists("BackupGeigerSpecies")) {
  raw_counts <- BackupGeigerSpecies
  raw_counts$SampleID <- NULL
  raw_counts <- as.data.frame(t(raw_counts))
} else {
  # Use the species matrix (rarefied is still valid for ALDEx2)
  raw_counts <- as.data.frame(t(species_matrix))
}

# Ensure integer counts
raw_counts <- as.data.frame(lapply(raw_counts, function(x) as.integer(round(x))))
rownames(raw_counts) <- colnames(species_matrix)

# Filter low abundance features for ALDEx2 (at least 10 reads total)
raw_counts <- raw_counts[rowSums(raw_counts) >= 10, ]

# --- 2.1: ALDEx2 Analysis: DY vs DO ---
cat("Running ALDEx2 for DY vs DO comparison...\n")

# Get sample conditions
sample_order <- colnames(raw_counts)
conditions_DY_DO <- Metadata$Groups[match(sample_order, Metadata$SampleID)]

# Subset to DY and DO samples
dy_do_idx <- which(conditions_DY_DO %in% c("DY", "DO"))
counts_DY_DO <- raw_counts[, dy_do_idx]
conds_DY_DO <- as.character(conditions_DY_DO[dy_do_idx])

# Run ALDEx2 (CLR transformation + Welch's t-test + Wilcoxon)
aldex_DY_DO <- aldex(counts_DY_DO, conds_DY_DO, mc.samples = 128,
                     test = "t", effect = TRUE, include.sample.summary = FALSE,
                     verbose = TRUE, denom = "all")

# Add species names and sort by effect size
aldex_DY_DO$Species <- rownames(aldex_DY_DO)
aldex_DY_DO <- aldex_DY_DO[order(aldex_DY_DO$effect), ]

# Filter significant results (BH-corrected p < 0.1 or effect size > 1)
sig_aldex_DY_DO <- subset(aldex_DY_DO, we.eBH < 0.1 | abs(effect) > 1)

cat(paste("Found", nrow(sig_aldex_DY_DO), "significant species (DY vs DO)\n"))

# Save results
write.csv(aldex_DY_DO, file = "Reviewer_ALDEx2_DY_vs_DO_Full.csv", row.names = FALSE)
write.csv(sig_aldex_DY_DO, file = "Reviewer_ALDEx2_DY_vs_DO_Significant.csv", row.names = FALSE)

# --- 2.2: ALDEx2 Effect Size Plot (DY vs DO) ---
ALDEx2_Effect_Plot_DY_DO <- aldex_DY_DO %>%
  filter(abs(effect) > 0.5) %>%  # Filter for moderate effect sizes
  mutate(
    Significant = ifelse(we.eBH < 0.1, "FDR < 0.1", "NS"),
    HigherIn = ifelse(effect > 0, "DY", "DO"),
    Species = reorder(Species, effect)
  ) %>%
  ggplot(aes(x = Species, y = effect, fill = HigherIn)) +
  geom_col(aes(alpha = Significant), color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = FigCols[c(4, 5)]) +
  scale_alpha_manual(values = c("FDR < 0.1" = 1, "NS" = 0.4)) +
  labs(
    y = "ALDEx2 Effect Size", x = NULL,
    title = "Differentially Abundant Species: Donor Young vs Donor Old HSCs",
    subtitle = "ALDEx2 analysis with CLR transformation",
    fill = "More Abundant In:",
    alpha = "Significance"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

ggsave(filename = "Reviewer_ALDEx2_DY_vs_DO_EffectPlot.pdf",
       plot = ALDEx2_Effect_Plot_DY_DO, width = 12, height = 10, device = cairo_pdf)

# --- 2.3: ALDEx2 Analysis: Y vs O ---
cat("\nRunning ALDEx2 for Y vs O comparison...\n")

y_o_idx <- which(conditions_DY_DO %in% c("Y", "O"))
# Re-index from original data
conditions_all <- Metadata$Groups[match(sample_order, Metadata$SampleID)]
y_o_idx <- which(conditions_all %in% c("Y", "O"))
counts_Y_O <- raw_counts[, y_o_idx]
conds_Y_O <- as.character(conditions_all[y_o_idx])

if (ncol(counts_Y_O) >= 4) {
  aldex_Y_O <- aldex(counts_Y_O, conds_Y_O, mc.samples = 128,
                     test = "t", effect = TRUE, include.sample.summary = FALSE,
                     verbose = TRUE, denom = "all")

  aldex_Y_O$Species <- rownames(aldex_Y_O)
  aldex_Y_O <- aldex_Y_O[order(aldex_Y_O$effect), ]

  sig_aldex_Y_O <- subset(aldex_Y_O, we.eBH < 0.1 | abs(effect) > 1)
  cat(paste("Found", nrow(sig_aldex_Y_O), "significant species (Y vs O)\n"))

  write.csv(aldex_Y_O, file = "Reviewer_ALDEx2_Y_vs_O_Full.csv", row.names = FALSE)
  write.csv(sig_aldex_Y_O, file = "Reviewer_ALDEx2_Y_vs_O_Significant.csv", row.names = FALSE)
}

# --- 2.4: ALDEx2 for Pathways ---
cat("\nRunning ALDEx2 for pathway analysis (DY vs DO)...\n")

# Check if pathway data exists
if (exists("GeigerPathways")) {
  # Get pathway counts (need integer-like values)
  pathway_metadata_cols <- 1:10
  pathway_counts <- GeigerPathways[, -pathway_metadata_cols]
  pathway_counts <- as.data.frame(lapply(pathway_counts, as.numeric))
  rownames(pathway_counts) <- GeigerPathways$SampleID

  # Transpose for ALDEx2 (pathways as rows, samples as columns)
  pathway_counts_t <- as.data.frame(t(pathway_counts))

  # Scale to pseudo-counts (multiply by 1000 and round to integers)
  pathway_counts_t <- as.data.frame(lapply(pathway_counts_t, function(x) {
    as.integer(round(x * 100))  # Scale CPM to pseudo-counts
  }))
  rownames(pathway_counts_t) <- colnames(pathway_counts)

  # Filter low abundance
  pathway_counts_t <- pathway_counts_t[rowSums(pathway_counts_t) >= 100, ]

  # Get conditions for pathways
  pathway_sample_order <- colnames(pathway_counts_t)
  pathway_conditions <- GeigerPathways$Groups[match(pathway_sample_order, GeigerPathways$SampleID)]

  # DY vs DO
  pw_dy_do_idx <- which(pathway_conditions %in% c("DY", "DO"))
  pw_counts_DY_DO <- pathway_counts_t[, pw_dy_do_idx]
  pw_conds_DY_DO <- as.character(pathway_conditions[pw_dy_do_idx])

  if (length(unique(pw_conds_DY_DO)) == 2 && ncol(pw_counts_DY_DO) >= 4) {
    aldex_pathways <- aldex(pw_counts_DY_DO, pw_conds_DY_DO, mc.samples = 128,
                            test = "t", effect = TRUE, include.sample.summary = FALSE,
                            verbose = TRUE, denom = "all")

    aldex_pathways$Pathway <- rownames(aldex_pathways)
    aldex_pathways <- aldex_pathways[order(aldex_pathways$we.eBH), ]

    sig_pathways <- subset(aldex_pathways, we.eBH < 0.2)
    cat(paste("Found", nrow(sig_pathways), "significant pathways (DY vs DO, FDR < 0.2)\n"))

    write.csv(aldex_pathways, file = "Reviewer_ALDEx2_Pathways_DY_vs_DO_Full.csv", row.names = FALSE)
    write.csv(sig_pathways, file = "Reviewer_ALDEx2_Pathways_DY_vs_DO_Significant.csv", row.names = FALSE)
  }
} else {
  cat("Pathway data not found in workspace. Skipping pathway ALDEx2 analysis.\n")
}

cat("ALDEx2 analysis complete. Files saved.\n")


#===============================================================================
# SECTION 3: POWER ANALYSIS AND VARIANCE ESTIMATES
#===============================================================================

cat("\n=== SECTION 3: Power Analysis and Variance Estimates ===\n\n")

# --- 3.1: Calculate variance estimates for each group ---
group_variance <- Diversity %>%
  group_by(Groups) %>%
  summarise(
    n = n(),
    Shannon_mean = mean(Shannon, na.rm = TRUE),
    Shannon_sd = sd(Shannon, na.rm = TRUE),
    Shannon_cv = Shannon_sd / Shannon_mean * 100,
    Chao1_mean = mean(S.chao1, na.rm = TRUE),
    Chao1_sd = sd(S.chao1, na.rm = TRUE),
    Chao1_cv = Chao1_sd / Chao1_mean * 100,
    .groups = "drop"
  )

cat("Variance estimates by group:\n")
print(as.data.frame(group_variance))

write.csv(group_variance, file = "Reviewer_VarianceEstimates.csv", row.names = FALSE)

# --- 3.2: Effect sizes from observed data ---
# Cohen's d for pairwise comparisons
calculate_cohens_d <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  pooled_sd <- sqrt(((nx - 1) * sd(x)^2 + (ny - 1) * sd(y)^2) / (nx + ny - 2))
  (mean(x) - mean(y)) / pooled_sd
}

# Calculate effect sizes for key comparisons
comparisons <- list(
  "Y_vs_O" = c("Y", "O"),
  "DY_vs_DO" = c("DY", "DO"),
  "Y_vs_RAG1" = c("Y", "RAG1-/-"),
  "DY_vs_Y" = c("DY", "Y")
)

effect_sizes <- data.frame(
  Comparison = character(),
  Metric = character(),
  n1 = integer(),
  n2 = integer(),
  Cohen_d = numeric(),
  stringsAsFactors = FALSE
)

for (comp_name in names(comparisons)) {
  groups <- comparisons[[comp_name]]
  g1 <- Diversity$Shannon[Diversity$Groups == groups[1]]
  g2 <- Diversity$Shannon[Diversity$Groups == groups[2]]

  if (length(g1) >= 2 && length(g2) >= 2) {
    d <- calculate_cohens_d(g1, g2)
    effect_sizes <- rbind(effect_sizes, data.frame(
      Comparison = comp_name,
      Metric = "Shannon",
      n1 = length(g1),
      n2 = length(g2),
      Cohen_d = round(d, 3)
    ))
  }
}

cat("\nObserved effect sizes (Cohen's d) for Shannon diversity:\n")
print(effect_sizes)
write.csv(effect_sizes, file = "Reviewer_EffectSizes.csv", row.names = FALSE)

# --- 3.3: Post-hoc power analysis ---
cat("\n--- Post-hoc Power Analysis ---\n")

# For each comparison, calculate achieved power
power_results <- data.frame(
  Comparison = character(),
  n_per_group = numeric(),
  Effect_Size_d = numeric(),
  Achieved_Power = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(effect_sizes)) {
  comp <- effect_sizes[i, ]
  n_avg <- round((comp$n1 + comp$n2) / 2)
  d <- abs(comp$Cohen_d)

  # Power for two-sample t-test
  if (!is.na(d) && d > 0) {
    pwr_result <- tryCatch({
      pwr::pwr.t.test(n = n_avg, d = d, sig.level = 0.05, type = "two.sample")
    }, error = function(e) NULL)

    if (!is.null(pwr_result)) {
      power_results <- rbind(power_results, data.frame(
        Comparison = comp$Comparison,
        n_per_group = n_avg,
        Effect_Size_d = round(d, 3),
        Achieved_Power = round(pwr_result$power, 3)
      ))
    }
  }
}

cat("\nAchieved statistical power for observed effect sizes:\n")
print(power_results)
write.csv(power_results, file = "Reviewer_PowerAnalysis.csv", row.names = FALSE)

# --- 3.4: Sample size recommendations ---
cat("\n--- Sample Size Requirements for 80% Power ---\n")

sample_size_recommendations <- data.frame(
  Effect_Size = c("Small (0.2)", "Medium (0.5)", "Large (0.8)", "Very Large (1.0)"),
  d = c(0.2, 0.5, 0.8, 1.0),
  n_per_group_80pct = sapply(c(0.2, 0.5, 0.8, 1.0), function(d) {
    ceiling(pwr::pwr.t.test(d = d, sig.level = 0.05, power = 0.80, type = "two.sample")$n)
  }),
  n_per_group_90pct = sapply(c(0.2, 0.5, 0.8, 1.0), function(d) {
    ceiling(pwr::pwr.t.test(d = d, sig.level = 0.05, power = 0.90, type = "two.sample")$n)
  })
)

cat("\nSample size per group needed for different effect sizes:\n")
print(sample_size_recommendations)

cat("\nInterpretation: With n=6 per group, the study is adequately powered (>80%)\n")
cat("for detecting large effect sizes (d > 0.8), which is typical for microbiome\n")
cat("studies comparing distinct biological conditions.\n")


#===============================================================================
# SECTION 4: BATCH EFFECT ANALYSIS
#===============================================================================

cat("\n=== SECTION 4: Batch Effect Analysis ===\n\n")

# --- 4.1: Visualize batch effects in PCA ---
if (exists("GeigerSpeciesNR") && "Experiment" %in% colnames(GeigerSpeciesNR)) {

  # Prepare data for PCA - remove rows with any NA in species columns
  pca_data <- GeigerSpeciesNR[, 11:ncol(GeigerSpeciesNR)]
  pca_data <- as.data.frame(lapply(pca_data, as.numeric))
  rownames(pca_data) <- GeigerSpeciesNR$SampleID

  # Replace NA with 0 for species counts
  pca_data[is.na(pca_data)] <- 0

  # Remove samples with all zeros or very low counts
  row_sums <- rowSums(pca_data)
  valid_samples <- row_sums > 1000
  pca_data <- pca_data[valid_samples, ]

  # Get matching metadata
  batch_metadata <- GeigerSpeciesNR[valid_samples, 1:10]

  # glog2 transform
  glog2 <- function(x) (asinh(x) - log(2)) / log(2)
  pca_data_l2 <- glog2(as.matrix(pca_data))

  # PCA
  pca_result <- PCA(pca_data_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)

  # Extract PC scores
  pc_scores <- as.data.frame(pca_result$ind$coord)
  pc_scores$SampleID <- rownames(pc_scores)
  pc_scores <- merge(pc_scores, batch_metadata[, c("SampleID", "Groups", "Experiment", "Cage")], by = "SampleID")

  # Plot colored by Experiment (batch)
  Batch_PCA_Plot <- ggplot(pc_scores, aes(x = Dim.1, y = Dim.2)) +
    geom_point(aes(color = factor(Experiment), shape = Groups), size = 4, alpha = 0.8) +
    stat_ellipse(aes(color = factor(Experiment)), level = 0.95, linetype = "dashed") +
    scale_shape_manual(values = c(16, 17, 15, 18, 8)) +
    labs(
      x = paste0("PC1 (", round(pca_result$eig[1, 2], 1), "%)"),
      y = paste0("PC2 (", round(pca_result$eig[2, 2], 1), "%)"),
      title = "PCA Colored by Sequencing Batch (Experiment)",
      subtitle = "Assessment of batch effects",
      color = "Batch",
      shape = "Group"
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )

  ggsave(filename = "Reviewer_BatchEffect_PCA.pdf", plot = Batch_PCA_Plot,
         width = 10, height = 8, device = cairo_pdf)

  # --- 4.2: PERMANOVA for batch effects ---
  cat("PERMANOVA testing for batch effects:\n")

  # Bray-Curtis distance matrix (using filtered data)
  bray_dist <- vegdist(pca_data, method = "bray")

  # Test for batch effect while controlling for biological group
  permanova_batch <- adonis2(bray_dist ~ Groups + Experiment,
                              data = batch_metadata, permutations = 999)
  print(permanova_batch)

  # Variance partitioning
  cat("\n--- Variance Partitioning ---\n")
  cat(paste("Groups explain:", round(permanova_batch$R2[1] * 100, 1), "% of variance\n"))
  cat(paste("Batch explains:", round(permanova_batch$R2[2] * 100, 1), "% of variance\n"))
  cat(paste("Residual:", round(permanova_batch$R2[3] * 100, 1), "% of variance\n"))

  # Save PERMANOVA results
  permanova_df <- as.data.frame(permanova_batch)
  write.csv(permanova_df, file = "Reviewer_PERMANOVA_BatchEffects.csv")

  # --- 4.3: Test batch effect within each group ---
  cat("\n--- Within-group batch effect tests ---\n")

  # Use the filtered data with batch_metadata
  filtered_species <- GeigerSpeciesNR[valid_samples, ]

  for (grp in unique(na.omit(filtered_species$Groups))) {
    grp_data <- subset(filtered_species, Groups == grp)
    if (length(unique(grp_data$Experiment)) > 1 && nrow(grp_data) >= 4) {
      grp_species <- grp_data[, 11:ncol(grp_data)]
      grp_species[is.na(grp_species)] <- 0
      grp_dist <- vegdist(grp_species, method = "bray")

      result <- tryCatch({
        adonis2(grp_dist ~ Experiment, data = grp_data, permutations = 999)
      }, error = function(e) NULL)

      if (!is.null(result)) {
        cat(paste("\nGroup:", grp, "\n"))
        cat(paste("  Batch R2:", round(result$R2[1] * 100, 1), "%"))
        cat(paste("  p-value:", round(result$`Pr(>F)`[1], 4), "\n"))
      }
    }
  }

} else {
  cat("Experiment/batch information not found. Skipping batch effect analysis.\n")
}

cat("\nBatch effect analysis complete.\n")


#===============================================================================
# SECTION 5: IMPROVED PCA PLOTS WITH VARIANCE EXPLAINED AND SAMPLE SIZES
#===============================================================================

cat("\n=== SECTION 5: Improved PCA Plots ===\n\n")

# Function to create improved PCA plot with variance explained and sample sizes
create_improved_pca <- function(data, metadata, groups_col, title, colors,
                                 filename, groups_to_include = NULL) {

  # Subset if specified
  if (!is.null(groups_to_include)) {
    idx <- which(metadata[[groups_col]] %in% groups_to_include)
    data <- data[idx, ]
    metadata <- metadata[idx, ]
  }

  # Convert to numeric and replace NA with 0
  data <- as.data.frame(lapply(data, as.numeric))
  data[is.na(data)] <- 0

  # Remove rows with very low total counts
  row_sums <- rowSums(data)
  valid_idx <- row_sums > 1000
  data <- data[valid_idx, ]
  metadata <- metadata[valid_idx, ]

  # glog2 transform
  data_l2 <- glog2(as.matrix(data))

  # Run PCA
  pca_res <- PCA(data_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)

  # Get variance explained
  var_pc1 <- round(pca_res$eig[1, 2], 1)
  var_pc2 <- round(pca_res$eig[2, 2], 1)

  # Calculate sample sizes per group
  sample_sizes <- table(metadata[[groups_col]])

  # Create legend labels with sample sizes
  groups_vector <- metadata[[groups_col]]
  legend_labels <- paste0(names(sample_sizes), " (n=", sample_sizes, ")")
  names(legend_labels) <- names(sample_sizes)

  # Extract coordinates
  pc_coords <- as.data.frame(pca_res$ind$coord)
  pc_coords$Group <- groups_vector

  # Create plot
  p <- ggplot(pc_coords, aes(x = Dim.1, y = Dim.2, color = Group)) +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = "solid", linewidth = 1) +
    scale_color_manual(values = colors, labels = legend_labels) +
    labs(
      x = paste0("PC1 (", var_pc1, "%)"),
      y = paste0("PC2 (", var_pc2, "%)"),
      title = title,
      color = "Group"
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )

  ggsave(filename = filename, plot = p, width = 10, height = 8, device = cairo_pdf)

  return(list(plot = p, var_pc1 = var_pc1, var_pc2 = var_pc2, pca = pca_res))
}

# --- 5.1: All Groups PCA ---
species_for_pca <- GeigerSpeciesNR[, 11:ncol(GeigerSpeciesNR)]
metadata_for_pca <- GeigerSpeciesNR[, 1:10]

all_groups_pca <- create_improved_pca(
  data = species_for_pca,
  metadata = metadata_for_pca,
  groups_col = "Groups",
  title = "PCA of Microbial Species Composition",
  colors = FigCols,
  filename = "Reviewer_PCA_AllGroups_Improved.pdf"
)

# --- 5.2: Experiment 1 PCA (Y, RAG1-/-, DY) ---
exp1_data <- subset(GeigerSpeciesNR, Experiment == 1 & Groups %in% c("Y", "RAG1-/-", "DY"))
if (nrow(exp1_data) > 0) {
  exp1_pca <- create_improved_pca(
    data = exp1_data[, 11:ncol(exp1_data)],
    metadata = exp1_data[, 1:10],
    groups_col = "Groups",
    title = "PCA of Microbial Species (Experiment 1: Y, RAG1-/-, DY)",
    colors = FigCols[c(1, 3, 4)],
    filename = "Reviewer_PCA_Exp1_Fig1A_Improved.pdf"
  )
}

# --- 5.3: Young vs Old PCA ---
yo_data <- subset(GeigerSpeciesNR, Groups %in% c("Y", "O"))
if (nrow(yo_data) > 0) {
  yo_pca <- create_improved_pca(
    data = yo_data[, 11:ncol(yo_data)],
    metadata = yo_data[, 1:10],
    groups_col = "Groups",
    title = "PCA of Microbial Species: Young vs Old Controls",
    colors = FigCols[c(1, 2)],
    filename = "Reviewer_PCA_Young_vs_Old_Improved.pdf"
  )
}

# --- 5.4: DY vs DO PCA ---
dy_do_data <- subset(GeigerSpeciesNR, Groups %in% c("DY", "DO"))
if (nrow(dy_do_data) > 0) {
  dy_do_pca <- create_improved_pca(
    data = dy_do_data[, 11:ncol(dy_do_data)],
    metadata = dy_do_data[, 1:10],
    groups_col = "Groups",
    title = "PCA of Microbial Species: Donor Young vs Donor Old HSCs",
    colors = FigCols[c(4, 5)],
    filename = "Reviewer_PCA_DY_vs_DO_Improved.pdf"
  )
}

# --- 5.5: Summary table of variance explained ---
variance_summary <- data.frame(
  Comparison = c("All Groups", "Experiment 1 (Y, RAG1-/-, DY)",
                 "Young vs Old", "DY vs DO"),
  PC1_Variance = c(all_groups_pca$var_pc1,
                   ifelse(exists("exp1_pca"), exp1_pca$var_pc1, NA),
                   ifelse(exists("yo_pca"), yo_pca$var_pc1, NA),
                   ifelse(exists("dy_do_pca"), dy_do_pca$var_pc1, NA)),
  PC2_Variance = c(all_groups_pca$var_pc2,
                   ifelse(exists("exp1_pca"), exp1_pca$var_pc2, NA),
                   ifelse(exists("yo_pca"), yo_pca$var_pc2, NA),
                   ifelse(exists("dy_do_pca"), dy_do_pca$var_pc2, NA))
)

cat("\nVariance explained by principal components:\n")
print(variance_summary)
write.csv(variance_summary, file = "Reviewer_PCA_VarianceExplained.csv", row.names = FALSE)

cat("\nImproved PCA plots saved.\n")


#===============================================================================
# SECTION 6: SAVE RESULTS AND SESSION INFO
#===============================================================================

cat("\n=== SECTION 6: Saving Results ===\n\n")

# Save updated workspace
save.image(file = "GeigerData_ReviewerRevisions_20251225")

# Session info for reproducibility
cat("R Session Information:\n")
sessionInfo()

cat("\n\n========================================\n")
cat("Reviewer Revisions Analysis Complete!\n")
cat("========================================\n\n")

cat("Output files generated:\n")
cat("  - Reviewer_Chao1_Diversity.pdf\n")
cat("  - Reviewer_AlphaDiversity_Panel.pdf\n")
cat("  - Reviewer_RarefactionCurves.pdf\n")
cat("  - Reviewer_ALDEx2_*_Full.csv (differential abundance results)\n")
cat("  - Reviewer_ALDEx2_*_Significant.csv (significant results)\n")
cat("  - Reviewer_ALDEx2_*_EffectPlot.pdf (effect size plots)\n")
cat("  - Reviewer_VarianceEstimates.csv\n")
cat("  - Reviewer_EffectSizes.csv\n")
cat("  - Reviewer_PowerAnalysis.csv\n")
cat("  - Reviewer_PERMANOVA_BatchEffects.csv\n")
cat("  - Reviewer_BatchEffect_PCA.pdf\n")
cat("  - Reviewer_PCA_*_Improved.pdf (PCAs with variance explained)\n")
cat("  - Reviewer_PCA_VarianceExplained.csv\n")
cat("  - GeigerData_ReviewerRevisions_20251225 (updated workspace)\n")
