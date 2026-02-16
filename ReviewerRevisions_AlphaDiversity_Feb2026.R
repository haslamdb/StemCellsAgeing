#-------------------------------------------------------------------------------
#
# Script: Alpha Diversity Plots - Reviewer Revisions (February 2026)
#
# Author:  David Haslam / Selina Stahl / Hartmut Geiger
# Date:    February 2026
#
# Description: Generate α-diversity plots (Shannon, Simpson, Chao1) for
#   specific group comparisons with exact color codes and significance brackets.
#   Also generates rarefaction curves in TIF format.
#
# Required input: GeigerData20250512 workspace (or GeigerData_ReviewerRevisions_20251225)
#
#-------------------------------------------------------------------------------

#===============================================================================
# SECTION 0: SETUP
#===============================================================================

# Load workspace - use the reviewer revisions workspace if available
if (file.exists("GeigerData_ReviewerRevisions_20251225")) {
  load("GeigerData_ReviewerRevisions_20251225")
  cat("Loaded reviewer revisions workspace.\n")
} else {
  load("GeigerData20250512")
  cat("Loaded main workspace. Will recompute diversity metrics.\n")
}

# Load packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggpubr)
library(rstatix)
library(pwr)

#===============================================================================
# SECTION 1: PREPARE DIVERSITY DATA
#===============================================================================

cat("\n=== Section 1: Preparing Diversity Data ===\n\n")

# Get species matrix (samples x species)
if (exists("Species")) {
  species_matrix <- as.data.frame(t(Species))
} else if (exists("GeigerSpeciesNR")) {
  species_matrix <- GeigerSpeciesNR[, 11:ncol(GeigerSpeciesNR)]
} else {
  stop("Cannot find species abundance data.")
}

species_matrix <- as.data.frame(lapply(species_matrix, as.numeric))
rownames(species_matrix) <- rownames(GeigerSpeciesNR)
species_matrix[is.na(species_matrix)] <- 0

# Calculate all diversity metrics if not already in Diversity dataframe
if (!exists("Diversity") || !"S.chao1" %in% colnames(Diversity)) {
  H <- diversity(species_matrix, "shannon")
  simp <- diversity(species_matrix, "simpson")
  S <- specnumber(species_matrix)
  J <- H / log(S)
  chao1_results <- t(estimateR(species_matrix))
  chao1_df <- as.data.frame(chao1_results)

  Diversity <- data.frame(
    SampleID = rownames(species_matrix),
    Shannon = H,
    Simpson = simp,
    SpeciesNo = S,
    Evenness = J,
    S.chao1 = chao1_df$S.chao1,
    se.chao1 = chao1_df$se.chao1
  )
  Diversity <- merge(Metadata, Diversity, by = "SampleID", all.x = TRUE)
}

# Verify Groups column
cat("Groups in data:\n")
print(table(Diversity$Groups))

#===============================================================================
# SECTION 2: ALPHA DIVERSITY PLOTS - COMPARISON-SPECIFIC
#===============================================================================

cat("\n=== Section 2: Comparison-Specific Alpha Diversity Plots ===\n\n")

# --- Define color palettes for each comparison ---
# Comparison A: Y vs DY vs RAG1-/-
colors_A <- c("Y" = "#A4DCFE", "DY" = "#FECC66", "RAG1-/-" = "#A5333A")
# Comparison B: Y vs O
colors_B <- c("Y" = "#A4DCFE", "O" = "#074080")
# Comparison C: DY vs DO
colors_C <- c("DY" = "#FECC66", "DO" = "#FD8008")

# --- Helper function to create alpha diversity boxplot with significance ---
create_alpha_plot <- function(data, groups_to_include, metric_col, metric_label,
                               color_map, comparisons_list = NULL,
                               y_expand_mult = 0.15) {

  # Subset data
  plot_data <- data %>%
    filter(Groups %in% groups_to_include) %>%
    mutate(Groups = factor(Groups, levels = names(color_map)))

  # Create base plot
  p <- ggplot(plot_data, aes(x = Groups, y = .data[[metric_col]])) +
    geom_boxplot(aes(color = Groups), fill = NA, lwd = 1,
                 outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(color = Groups), width = 0.15, size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 5,
                 color = "black") +
    scale_color_manual(values = color_map) +
    labs(x = NULL, y = metric_label) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 16, color = "black", face = "bold"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 18),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, y_expand_mult)))

  # Add significance brackets if comparisons specified
  if (!is.null(comparisons_list)) {
    p <- p + stat_compare_means(
      comparisons = comparisons_list,
      method = "wilcox.test",
      label = "p.signif",       # Shows *, **, ***, ns
      bracket.size = 0.6,
      tip.length = 0.02,
      size = 5,
      step.increase = 0.08
    )
  }

  return(p)
}

# --- Define pairwise comparisons for each group set ---
comparisons_A <- list(c("Y", "DY"), c("Y", "RAG1-/-"), c("DY", "RAG1-/-"))
comparisons_B <- list(c("Y", "O"))
comparisons_C <- list(c("DY", "DO"))

# --- Metrics to plot ---
metrics <- list(
  list(col = "Shannon",  label = "Shannon Index"),
  list(col = "Simpson",  label = "Simpson Index"),
  list(col = "S.chao1",  label = "Chao1 Richness Estimator")
)

#===============================================================================
# COMPARISON A: Y vs DY vs RAG1-/-
#===============================================================================

cat("Generating Comparison A: Y vs DY vs RAG1-/-\n")

plots_A <- lapply(metrics, function(m) {
  create_alpha_plot(
    data = Diversity,
    groups_to_include = c("Y", "DY", "RAG1-/-"),
    metric_col = m$col,
    metric_label = m$label,
    color_map = colors_A,
    comparisons_list = comparisons_A,
    y_expand_mult = 0.25  # Extra space for brackets
  )
})

# Combined panel (3 metrics side by side)
panel_A <- plot_grid(plots_A[[1]], plots_A[[2]], plots_A[[3]],
                     ncol = 3, labels = c("A", "B", "C"),
                     label_size = 18, align = "h")

ggsave("Alpha_Diversity_Y_DY_RAG1.pdf", panel_A,
       width = 16, height = 6, device = cairo_pdf)
ggsave("Alpha_Diversity_Y_DY_RAG1.tif", panel_A,
       width = 16, height = 6, dpi = 300, compression = "lzw")

# Also save individual metric plots
for (i in seq_along(metrics)) {
  fname <- paste0("Alpha_", metrics[[i]]$col, "_Y_DY_RAG1")
  ggsave(paste0(fname, ".pdf"), plots_A[[i]],
         width = 6, height = 6, device = cairo_pdf)
  ggsave(paste0(fname, ".tif"), plots_A[[i]],
         width = 6, height = 6, dpi = 300, compression = "lzw")
}

#===============================================================================
# COMPARISON B: Y vs O
#===============================================================================

cat("Generating Comparison B: Y vs O\n")

plots_B <- lapply(metrics, function(m) {
  create_alpha_plot(
    data = Diversity,
    groups_to_include = c("Y", "O"),
    metric_col = m$col,
    metric_label = m$label,
    color_map = colors_B,
    comparisons_list = comparisons_B,
    y_expand_mult = 0.15
  )
})

panel_B <- plot_grid(plots_B[[1]], plots_B[[2]], plots_B[[3]],
                     ncol = 3, labels = c("A", "B", "C"),
                     label_size = 18, align = "h")

ggsave("Alpha_Diversity_Y_vs_O.pdf", panel_B,
       width = 16, height = 6, device = cairo_pdf)
ggsave("Alpha_Diversity_Y_vs_O.tif", panel_B,
       width = 16, height = 6, dpi = 300, compression = "lzw")

for (i in seq_along(metrics)) {
  fname <- paste0("Alpha_", metrics[[i]]$col, "_Y_vs_O")
  ggsave(paste0(fname, ".pdf"), plots_B[[i]],
         width = 6, height = 6, device = cairo_pdf)
  ggsave(paste0(fname, ".tif"), plots_B[[i]],
         width = 6, height = 6, dpi = 300, compression = "lzw")
}

#===============================================================================
# COMPARISON C: DY vs DO
#===============================================================================

cat("Generating Comparison C: DY vs DO\n")

plots_C <- lapply(metrics, function(m) {
  create_alpha_plot(
    data = Diversity,
    groups_to_include = c("DY", "DO"),
    metric_col = m$col,
    metric_label = m$label,
    color_map = colors_C,
    comparisons_list = comparisons_C,
    y_expand_mult = 0.15
  )
})

panel_C <- plot_grid(plots_C[[1]], plots_C[[2]], plots_C[[3]],
                     ncol = 3, labels = c("A", "B", "C"),
                     label_size = 18, align = "h")

ggsave("Alpha_Diversity_DY_vs_DO.pdf", panel_C,
       width = 16, height = 6, device = cairo_pdf)
ggsave("Alpha_Diversity_DY_vs_DO.tif", panel_C,
       width = 16, height = 6, dpi = 300, compression = "lzw")

for (i in seq_along(metrics)) {
  fname <- paste0("Alpha_", metrics[[i]]$col, "_DY_vs_DO")
  ggsave(paste0(fname, ".pdf"), plots_C[[i]],
         width = 6, height = 6, device = cairo_pdf)
  ggsave(paste0(fname, ".tif"), plots_C[[i]],
         width = 6, height = 6, dpi = 300, compression = "lzw")
}

#===============================================================================
# SECTION 3: STATISTICAL SUMMARY TABLE
#===============================================================================

cat("\n=== Section 3: Statistical Summary ===\n\n")

# Pairwise Wilcoxon tests for all metrics across all comparisons
all_comparisons <- list(
  "Y_vs_DY_vs_RAG1" = list(
    groups = c("Y", "DY", "RAG1-/-"),
    pairs = comparisons_A
  ),
  "Y_vs_O" = list(
    groups = c("Y", "O"),
    pairs = comparisons_B
  ),
  "DY_vs_DO" = list(
    groups = c("DY", "DO"),
    pairs = comparisons_C
  )
)

stat_results <- data.frame()

for (comp_name in names(all_comparisons)) {
  comp <- all_comparisons[[comp_name]]
  sub_data <- Diversity %>% filter(Groups %in% comp$groups)

  for (m in metrics) {
    for (pair in comp$pairs) {
      g1_vals <- sub_data[[m$col]][sub_data$Groups == pair[1]]
      g2_vals <- sub_data[[m$col]][sub_data$Groups == pair[2]]

      # Remove NAs
      g1_vals <- g1_vals[!is.na(g1_vals)]
      g2_vals <- g2_vals[!is.na(g2_vals)]

      if (length(g1_vals) >= 2 && length(g2_vals) >= 2) {
        wt <- wilcox.test(g1_vals, g2_vals)
        stat_results <- rbind(stat_results, data.frame(
          Comparison_Set = comp_name,
          Group1 = pair[1],
          Group2 = pair[2],
          Metric = m$col,
          n1 = length(g1_vals),
          n2 = length(g2_vals),
          Mean1 = round(mean(g1_vals), 3),
          Mean2 = round(mean(g2_vals), 3),
          p_value = signif(wt$p.value, 4),
          Significance = ifelse(wt$p.value < 0.001, "***",
                         ifelse(wt$p.value < 0.01, "**",
                         ifelse(wt$p.value < 0.05, "*", "ns"))),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

cat("Alpha Diversity Statistical Comparisons:\n")
print(stat_results)
write.csv(stat_results, "Alpha_Diversity_Statistics.csv", row.names = FALSE)


#===============================================================================
# SECTION 4: RAREFACTION CURVES (TIF FORMAT)
#===============================================================================

cat("\n=== Section 4: Rarefaction Curves (TIF Format) ===\n\n")

# Use the same color codes as comparison A (all groups use Y, O, RAG1-/-, DY, DO)
rarefaction_colors <- c("Y" = "#A4DCFE", "O" = "#074080",
                         "RAG1-/-" = "#A5333A", "DY" = "#FECC66", "DO" = "#FD8008")

# Compute rarefaction data
compute_rarefaction <- function(sample_counts, sample_id, depths = NULL) {
  sample_counts <- as.numeric(sample_counts)
  sample_counts[is.na(sample_counts)] <- 0
  total_reads <- sum(sample_counts)
  if (total_reads == 0) return(NULL)

  if (is.null(depths)) {
    depths <- seq(10000, min(total_reads, 750000), by = 30000)
  }
  depths <- depths[depths <= total_reads & depths > 0]
  if (length(depths) == 0) return(NULL)

  rarefied_richness <- sapply(depths, function(d) {
    tryCatch(rarefy(sample_counts, sample = d), error = function(e) NA)
  })

  data.frame(SampleID = sample_id, Depth = depths, Richness = rarefied_richness)
}

set.seed(42)
rarefaction_list <- lapply(1:nrow(species_matrix), function(i) {
  tryCatch(
    compute_rarefaction(species_matrix[i, ], rownames(species_matrix)[i],
                        depths = seq(10000, 750000, by = 30000)),
    error = function(e) NULL
  )
})

rarefaction_list <- rarefaction_list[!sapply(rarefaction_list, is.null)]
rarefaction_data <- do.call(rbind, rarefaction_list)

# Add group info
rarefaction_data <- merge(rarefaction_data,
                           Metadata[, c("SampleID", "Groups")],
                           by = "SampleID")

# Calculate sample sizes per group for legend
sample_sizes <- rarefaction_data %>%
  distinct(SampleID, Groups) %>%
  count(Groups) %>%
  mutate(label = paste0(Groups, " (n=", n, ")"))

# Create legend labels
legend_labels <- setNames(sample_sizes$label, sample_sizes$Groups)

# Plot
Rarefaction_Plot <- ggplot(rarefaction_data,
                            aes(x = Depth, y = Richness,
                                group = SampleID, color = Groups)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  scale_color_manual(values = rarefaction_colors, labels = legend_labels) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Sequencing Depth (Reads)",
    y = "Species Richness",
    color = "Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 16)
  ) +
  geom_vline(xintercept = 750000, linetype = "dashed", color = "gray40") +
  annotate("text", x = 720000,
           y = max(rarefaction_data$Richness, na.rm = TRUE) * 0.95,
           label = "Rarefaction\ndepth", hjust = 1, size = 3.5)

# Save as TIF (publication quality)
ggsave("Rarefaction_Curves.tif", Rarefaction_Plot,
       width = 10, height = 7, dpi = 300, compression = "lzw")
# Also save PDF for convenience
ggsave("Rarefaction_Curves.pdf", Rarefaction_Plot,
       width = 10, height = 7, device = cairo_pdf)

cat("Rarefaction curves saved (TIF and PDF).\n")


#===============================================================================
# SECTION 5: POWER ANALYSIS - INCLUDING DY vs DO
#===============================================================================

cat("\n=== Section 5: Power Analysis (Updated with DY vs DO) ===\n\n")

# Cohen's d function
calculate_cohens_d <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  nx <- length(x)
  ny <- length(y)
  if (nx < 2 || ny < 2) return(NA)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  if (pooled_sd == 0) return(NA)
  (mean(x) - mean(y)) / pooled_sd
}

# Comparisons for power analysis (as requested: Y vs O and DY vs DO)
power_comparisons <- list(
  "Y_vs_O" = c("Y", "O"),
  "DY_vs_DO" = c("DY", "DO")
)

# Calculate effect sizes and power for multiple metrics
power_metrics <- c("Shannon", "Simpson", "S.chao1")
power_results_all <- data.frame()

for (comp_name in names(power_comparisons)) {
  groups <- power_comparisons[[comp_name]]
  for (metric in power_metrics) {
    g1 <- Diversity[[metric]][Diversity$Groups == groups[1]]
    g2 <- Diversity[[metric]][Diversity$Groups == groups[2]]
    g1 <- g1[!is.na(g1)]
    g2 <- g2[!is.na(g2)]

    d <- calculate_cohens_d(g1, g2)
    n_avg <- round((length(g1) + length(g2)) / 2)

    achieved_power <- NA
    if (!is.na(d) && abs(d) > 0 && n_avg >= 2) {
      pwr_result <- tryCatch(
        pwr.t.test(n = n_avg, d = abs(d), sig.level = 0.05, type = "two.sample"),
        error = function(e) NULL
      )
      if (!is.null(pwr_result)) achieved_power <- round(pwr_result$power, 3)
    }

    power_results_all <- rbind(power_results_all, data.frame(
      Comparison = comp_name,
      Metric = metric,
      n1 = length(g1),
      n2 = length(g2),
      Mean1 = round(mean(g1), 3),
      Mean2 = round(mean(g2), 3),
      SD1 = round(sd(g1), 3),
      SD2 = round(sd(g2), 3),
      Cohen_d = round(d, 3),
      Achieved_Power = achieved_power,
      stringsAsFactors = FALSE
    ))
  }
}

cat("Power Analysis Results (Y vs O and DY vs DO):\n")
print(power_results_all)
write.csv(power_results_all, "Power_Analysis_YvsO_DYvsDO.csv", row.names = FALSE)

# Also compute: what sample size would be needed for 80% power?
cat("\n--- Sample Size Requirements for 80% Power ---\n")
for (i in 1:nrow(power_results_all)) {
  row <- power_results_all[i, ]
  if (!is.na(row$Cohen_d) && abs(row$Cohen_d) > 0) {
    n_needed <- tryCatch(
      ceiling(pwr.t.test(d = abs(row$Cohen_d), sig.level = 0.05,
                          power = 0.80, type = "two.sample")$n),
      error = function(e) NA
    )
    cat(sprintf("  %s (%s): d=%.3f, need n=%s per group for 80%% power\n",
                row$Comparison, row$Metric, row$Cohen_d,
                ifelse(is.na(n_needed), "NA", n_needed)))
  }
}

# PERMDISP / Bray-Curtis based effect sizes for beta diversity
cat("\n--- Beta Diversity Effect Sizes ---\n")
cat("(For PERMANOVA-based power, use observed R2 values from PERMANOVA output)\n")


#===============================================================================
# SECTION 6: SESSION INFO
#===============================================================================

cat("\n\n========================================\n")
cat("Reviewer Revisions (Feb 2026) Complete!\n")
cat("========================================\n\n")

cat("Output files generated:\n")
cat("  Alpha Diversity Panels:\n")
cat("    - Alpha_Diversity_Y_DY_RAG1.pdf / .tif\n")
cat("    - Alpha_Diversity_Y_vs_O.pdf / .tif\n")
cat("    - Alpha_Diversity_DY_vs_DO.pdf / .tif\n")
cat("  Individual Alpha Diversity Plots:\n")
cat("    - Alpha_Shannon_*.pdf / .tif\n")
cat("    - Alpha_Simpson_*.pdf / .tif\n")
cat("    - Alpha_S.chao1_*.pdf / .tif\n")
cat("  Rarefaction:\n")
cat("    - Rarefaction_Curves.tif / .pdf\n")
cat("  Statistics:\n")
cat("    - Alpha_Diversity_Statistics.csv\n")
cat("    - Power_Analysis_YvsO_DYvsDO.csv\n")

sessionInfo()
