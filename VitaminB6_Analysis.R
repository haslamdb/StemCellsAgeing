#-------------------------------------------------------------------------------
#
# Script: Vitamin B6 Pathway and Species Analysis
#
# Author:  David Haslam
# Date:    December 25, 2025
#
# Description: Analyze pathways and species related to vitamin B6 (pyridoxine)
#              metabolism in the context of HSC aging
#
#-------------------------------------------------------------------------------

#===============================================================================
# SECTION 0: LOAD WORKSPACE AND PACKAGES
#===============================================================================

# Load the workspace from reviewer revisions (or the main analysis)
if (file.exists("GeigerData_ReviewerRevisions_20251225")) {
  load("GeigerData_ReviewerRevisions_20251225")
} else {
  load("GeigerData20250512")
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(rstatix)

# Color palette
FigCols <- c("#A4DCFE", "#074080", "#B6474B", "#F9CC66", "#FD8008")

#===============================================================================
# SECTION 1: IDENTIFY VITAMIN B6 RELATED PATHWAYS
#===============================================================================

cat("\n=== SECTION 1: Vitamin B6 Related Pathways ===\n\n")

# Common vitamin B6 pathway identifiers in MetaCyc/HUMAnN3:
# - PYRIDOXSYN-PWY: Pyridoxal 5'-phosphate biosynthesis I
# - PWY-6892: Pyridoxal 5'-phosphate salvage II (plants and cyanobacteria)
# - PYRIDOX-PWY: Pyridoxine biosynthesis I
# - PWY-6467: Pyridoxine biosynthesis IV (Archaebacteria)
# - PWY-6466: Pyridoxine biosynthesis III (Actinobacteria)
# - PWY-6464: Pyridoxine biosynthesis I (archaea)
# - PWY-6465: Pyridoxine biosynthesis II
# - PLPSAL-PWY: Pyridoxal 5'-phosphate salvage I
# - PWY-6895: Vitamin B6 degradation
# - Also look for: pyridox, B6, PLP, pyridoxal, pyridoxamine

if (exists("GeigerPathways")) {
  cat("Searching for vitamin B6 related pathways...\n\n")

  # Get all pathway names
  pathway_cols <- colnames(GeigerPathways)[11:ncol(GeigerPathways)]

  # Search patterns for B6-related pathways
  b6_patterns <- c("PYRIDOX", "B6", "PLP", "pyridox", "pyridoxal",
                   "pyridoxamine", "pyridoxine", "vitamin.b6")

  # Find matching pathways
  b6_pathways <- pathway_cols[grepl(paste(b6_patterns, collapse = "|"),
                                     pathway_cols, ignore.case = TRUE)]

  cat("Found", length(b6_pathways), "vitamin B6-related pathways:\n")
  for (pw in b6_pathways) {
    cat("  -", pw, "\n")
  }

  if (length(b6_pathways) > 0) {
    # Extract B6 pathway data
    b6_data <- GeigerPathways[, c("SampleID", "Groups", "Experiment", b6_pathways)]

    # Standardize group names to match species data
    b6_data$Groups <- as.character(b6_data$Groups)
    b6_data$Groups <- gsub("Young_C57BL-Untransplanted", "Y", b6_data$Groups)
    b6_data$Groups <- gsub("Old_C57BL-Untransplanted", "O", b6_data$Groups)
    b6_data$Groups <- gsub("Young_RAG1-Untransplanted", "RAG1-/-", b6_data$Groups)
    b6_data$Groups <- gsub("RAG1-Young_HSCs", "DY", b6_data$Groups)
    b6_data$Groups <- gsub("RAG1-Old_HSCs", "DO", b6_data$Groups)
    b6_data$Groups <- gsub("RAG1-Rejuvenated_HSCs", "Rejuv", b6_data$Groups)

    # Convert to long format for plotting
    b6_long <- b6_data %>%
      pivot_longer(cols = all_of(b6_pathways),
                   names_to = "Pathway", values_to = "Abundance") %>%
      mutate(Abundance = as.numeric(Abundance))

    # Clean pathway names for plotting
    b6_long$PathwayShort <- gsub("\\.\\..*", "", b6_long$Pathway)
    b6_long$PathwayShort <- gsub("\\.", " ", b6_long$PathwayShort)

    # Filter to only include valid groups and remove NA
    b6_long <- b6_long %>%
      filter(!is.na(Groups) & Groups %in% c("Y", "O", "RAG1-/-", "DY", "DO"))

    # --- 1.1: Plot B6 pathway abundances by group ---
    B6_Pathway_Plot <- ggplot(b6_long, aes(x = Groups, y = Abundance, fill = Groups)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
      facet_wrap(~ PathwayShort, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = FigCols, drop = FALSE) +
      labs(
        title = "Vitamin B6 Pathway Abundance by Group",
        x = NULL, y = "Abundance (CPM)"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8),
        legend.position = "none"
      )

    ggsave("VitaminB6_Pathway_Abundance.pdf", B6_Pathway_Plot,
           width = 12, height = max(6, length(b6_pathways) * 1.5), device = cairo_pdf)

    # --- 1.2: Statistical tests for B6 pathways ---
    cat("\n--- Statistical Tests for B6 Pathways ---\n")

    b6_stats <- b6_long %>%
      group_by(Pathway) %>%
      kruskal_test(Abundance ~ Groups) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance()

    print(as.data.frame(b6_stats))
    write.csv(b6_stats, "VitaminB6_Pathway_KruskalWallis.csv", row.names = FALSE)

    # Pairwise comparisons for significant pathways
    cat("\n--- Pairwise Comparisons ---\n")

    b6_pairwise <- tryCatch({
      b6_long %>%
        filter(!is.na(Abundance) & is.finite(Abundance)) %>%
        group_by(Pathway) %>%
        wilcox_test(Abundance ~ Groups, p.adjust.method = "BH") %>%
        filter(p.adj < 0.1)
    }, error = function(e) {
      cat("Note: Some pairwise comparisons failed due to insufficient data\n")
      cat("Running individual pathway comparisons instead...\n\n")

      # Manual pairwise comparisons
      results <- data.frame()
      for (pw in unique(b6_long$Pathway)) {
        pw_data <- b6_long[b6_long$Pathway == pw & !is.na(b6_long$Abundance), ]
        for (g1 in c("Y", "O", "RAG1-/-")) {
          for (g2 in c("DY", "DO")) {
            d1 <- pw_data$Abundance[pw_data$Groups == g1]
            d2 <- pw_data$Abundance[pw_data$Groups == g2]
            if (length(d1) >= 2 && length(d2) >= 2) {
              test <- tryCatch(wilcox.test(d1, d2), error = function(e) NULL)
              if (!is.null(test)) {
                results <- rbind(results, data.frame(
                  Pathway = pw, group1 = g1, group2 = g2,
                  p = test$p.value
                ))
              }
            }
          }
        }
      }
      if (nrow(results) > 0) {
        results$p.adj <- p.adjust(results$p, method = "BH")
        results <- results[results$p.adj < 0.1, ]
      }
      return(results)
    })

    if (nrow(b6_pairwise) > 0) {
      print(as.data.frame(b6_pairwise))
      write.csv(b6_pairwise, "VitaminB6_Pathway_Pairwise.csv", row.names = FALSE)
    } else {
      cat("No significant pairwise comparisons (adjusted p < 0.1)\n")
    }

    # --- 1.3: Focus on DY vs DO comparison ---
    cat("\n--- DY vs DO Comparison for B6 Pathways ---\n")

    b6_DY_DO <- tryCatch({
      b6_long %>%
        filter(Groups %in% c("DY", "DO") & !is.na(Abundance) & is.finite(Abundance)) %>%
        group_by(Pathway) %>%
        wilcox_test(Abundance ~ Groups) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
    }, error = function(e) {
      cat("Running manual comparison for DY vs DO...\n")
      results <- data.frame()
      for (pw in unique(b6_long$Pathway)) {
        pw_data <- b6_long[b6_long$Pathway == pw & !is.na(b6_long$Abundance), ]
        dy <- pw_data$Abundance[pw_data$Groups == "DY"]
        do <- pw_data$Abundance[pw_data$Groups == "DO"]
        if (length(dy) >= 2 && length(do) >= 2) {
          test <- tryCatch(wilcox.test(dy, do), error = function(e) NULL)
          if (!is.null(test)) {
            results <- rbind(results, data.frame(
              Pathway = pw,
              group1 = "DY", group2 = "DO",
              n1 = length(dy), n2 = length(do),
              statistic = test$statistic,
              p = test$p.value
            ))
          }
        }
      }
      if (nrow(results) > 0) {
        results$p.adj <- p.adjust(results$p, method = "BH")
      }
      return(results)
    })

    if (nrow(b6_DY_DO) > 0) {
      print(as.data.frame(b6_DY_DO))
      write.csv(b6_DY_DO, "VitaminB6_DY_vs_DO.csv", row.names = FALSE)
    }

  } else {
    cat("\nNo vitamin B6 pathways found in the dataset.\n")
    cat("This may indicate:\n")
    cat("  1. Low abundance below detection threshold\n")
    cat("  2. Different pathway naming convention in your HUMAnN3 output\n")
  }

} else {
  cat("GeigerPathways not found. Loading pathway data...\n")
}


#===============================================================================
# SECTION 2: SPECIES KNOWN TO PRODUCE VITAMIN B6
#===============================================================================

cat("\n=== SECTION 2: Vitamin B6 Producing Species ===\n\n")

# Key species known to synthesize vitamin B6:
# - Bacteroides species (B. fragilis, B. thetaiotaomicron)
# - Prevotella species
# - Bifidobacterium species (B. longum, B. adolescentis)
# - Lactobacillus species
# - Enterococcus faecium
# - Streptococcus thermophilus

b6_producer_patterns <- c(
  "Bacteroides",
  "Prevotella",
  "Bifidobacterium",
  "Lactobacillus",
  "Enterococcus",
  "Streptococcus"
)

if (exists("GeigerSpeciesNR") || exists("Species")) {

  # Get species data
  if (exists("Species")) {
    species_cols <- rownames(Species)
    species_data <- as.data.frame(t(Species))
    species_data$SampleID <- rownames(species_data)
  } else {
    species_cols <- colnames(GeigerSpeciesNR)[11:ncol(GeigerSpeciesNR)]
    species_data <- GeigerSpeciesNR
  }

  # Find B6 producer species
  b6_species <- species_cols[grepl(paste(b6_producer_patterns, collapse = "|"),
                                    species_cols, ignore.case = TRUE)]

  cat("Found", length(b6_species), "potential vitamin B6 producing species\n")

  if (length(b6_species) > 0) {
    # Merge with metadata if needed
    if (!"Groups" %in% colnames(species_data)) {
      species_data <- merge(Metadata[, c("SampleID", "Groups", "Experiment")],
                            species_data, by = "SampleID")
    }

    # Calculate total B6 producer abundance
    b6_producer_cols <- intersect(b6_species, colnames(species_data))
    species_data$B6_Producers_Total <- rowSums(
      species_data[, b6_producer_cols, drop = FALSE], na.rm = TRUE
    )

    # Filter to valid groups
    species_data <- species_data %>%
      filter(!is.na(Groups) & Groups %in% c("Y", "O", "RAG1-/-", "DY", "DO"))

    # --- 2.1: Plot total B6 producer abundance ---
    B6_Producers_Plot <- ggplot(species_data, aes(x = Groups, y = B6_Producers_Total)) +
      geom_boxplot(aes(fill = Groups), outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
      scale_fill_manual(values = FigCols, drop = FALSE) +
      labs(
        title = "Total Abundance of Vitamin B6 Producing Species",
        subtitle = paste("Includes:", paste(b6_producer_patterns, collapse = ", ")),
        x = NULL,
        y = "Total Read Count"
      ) +
      theme_bw() +
      theme(legend.position = "none")

    ggsave("VitaminB6_Producers_Total.pdf", B6_Producers_Plot,
           width = 8, height = 6, device = cairo_pdf)

    # --- 2.2: Statistical test ---
    cat("\n--- Kruskal-Wallis Test for B6 Producers ---\n")
    kruskal_b6 <- kruskal.test(B6_Producers_Total ~ Groups, data = species_data)
    print(kruskal_b6)

    cat("\n--- Pairwise Wilcoxon Tests ---\n")
    pairwise_b6 <- pairwise.wilcox.test(species_data$B6_Producers_Total,
                                         species_data$Groups,
                                         p.adjust.method = "BH")
    print(pairwise_b6)

    # --- 2.3: Individual B6 producer species ---
    cat("\n--- Individual B6 Producer Species ---\n")

    # Get top 10 most abundant B6 producers
    b6_abundances <- colMeans(species_data[, b6_producer_cols, drop = FALSE], na.rm = TRUE)
    top_b6_species <- names(sort(b6_abundances, decreasing = TRUE))[1:min(15, length(b6_abundances))]

    cat("Top vitamin B6 producing species by mean abundance:\n")
    for (sp in top_b6_species) {
      cat(sprintf("  %s: %.1f reads\n", sp, b6_abundances[sp]))
    }

    # Create long format for top species
    b6_species_long <- species_data %>%
      select(SampleID, Groups, all_of(top_b6_species)) %>%
      pivot_longer(cols = all_of(top_b6_species),
                   names_to = "Species", values_to = "Abundance") %>%
      mutate(Abundance = as.numeric(Abundance))

    # Filter to valid groups
    b6_species_long <- b6_species_long %>%
      filter(!is.na(Groups) & Groups %in% c("Y", "O", "RAG1-/-", "DY", "DO"))

    # Plot individual species
    B6_Species_Individual <- ggplot(b6_species_long, aes(x = Groups, y = Abundance + 1)) +
      geom_boxplot(aes(fill = Groups), outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
      facet_wrap(~ Species, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = FigCols, drop = FALSE) +
      scale_y_log10() +
      labs(
        title = "Individual Vitamin B6 Producing Species",
        x = NULL, y = "Abundance (log10 + 1)"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 7),
        legend.position = "none"
      )

    ggsave("VitaminB6_Producers_Individual.pdf", B6_Species_Individual,
           width = 14, height = 12, device = cairo_pdf)

    # --- 2.4: DY vs DO comparison for B6 producers ---
    cat("\n--- DY vs DO Comparison for B6 Producers ---\n")

    dy_do_species <- subset(species_data, Groups %in% c("DY", "DO"))
    wilcox_dy_do <- wilcox.test(B6_Producers_Total ~ Groups, data = dy_do_species)
    cat(sprintf("Wilcoxon test (DY vs DO): W = %.1f, p = %.4f\n",
                wilcox_dy_do$statistic, wilcox_dy_do$p.value))

    # Effect size
    dy_mean <- mean(dy_do_species$B6_Producers_Total[dy_do_species$Groups == "DY"], na.rm = TRUE)
    do_mean <- mean(dy_do_species$B6_Producers_Total[dy_do_species$Groups == "DO"], na.rm = TRUE)
    cat(sprintf("Mean abundance - DY: %.1f, DO: %.1f (fold change: %.2f)\n",
                dy_mean, do_mean, dy_mean / do_mean))

  }
}


#===============================================================================
# SECTION 3: GENE FAMILIES RELATED TO VITAMIN B6
#===============================================================================

cat("\n=== SECTION 3: Vitamin B6 Related Gene Families ===\n\n")

# Key B6 biosynthesis genes:
# - pdxA, pdxB, pdxJ, pdxK, pdxH, pdxS, pdxT (E. coli pathway)
# - SNZ1, SNO1 (yeast)
# - PYRIDOX*, PLP*

if (exists("GeigerGenes")) {
  gene_cols <- colnames(GeigerGenes)[11:ncol(GeigerGenes)]

  b6_gene_patterns <- c("pdxA", "pdxB", "pdxJ", "pdxK", "pdxH", "pdxS", "pdxT",
                        "PYRIDOX", "PLP", "pyridox", "B6")

  b6_genes <- gene_cols[grepl(paste(b6_gene_patterns, collapse = "|"),
                               gene_cols, ignore.case = TRUE)]

  cat("Found", length(b6_genes), "vitamin B6-related gene families\n")

  if (length(b6_genes) > 0) {
    cat("\nB6-related genes found:\n")
    for (g in head(b6_genes, 20)) {
      cat("  -", g, "\n")
    }

    # Extract and analyze
    b6_gene_data <- GeigerGenes[, c("SampleID", "Groups", b6_genes)]

    # Standardize group names
    b6_gene_data$Groups <- as.character(b6_gene_data$Groups)
    b6_gene_data$Groups <- gsub("Young_C57BL-Untransplanted", "Y", b6_gene_data$Groups)
    b6_gene_data$Groups <- gsub("Old_C57BL-Untransplanted", "O", b6_gene_data$Groups)
    b6_gene_data$Groups <- gsub("Young_RAG1-Untransplanted", "RAG1-/-", b6_gene_data$Groups)
    b6_gene_data$Groups <- gsub("RAG1-Young_HSCs", "DY", b6_gene_data$Groups)
    b6_gene_data$Groups <- gsub("RAG1-Old_HSCs", "DO", b6_gene_data$Groups)
    b6_gene_data$Groups <- gsub("RAG1-Rejuvenated_HSCs", "Rejuv", b6_gene_data$Groups)

    # Calculate total B6 gene abundance
    b6_gene_data$B6_Genes_Total <- rowSums(
      sapply(b6_gene_data[, b6_genes], as.numeric), na.rm = TRUE
    )

    # Filter to valid groups
    b6_gene_data <- b6_gene_data %>%
      filter(!is.na(Groups) & Groups %in% c("Y", "O", "RAG1-/-", "DY", "DO"))

    # Plot
    B6_Genes_Plot <- ggplot(b6_gene_data, aes(x = Groups, y = B6_Genes_Total)) +
      geom_boxplot(aes(fill = Groups), outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
      scale_fill_manual(values = FigCols, drop = FALSE) +
      labs(
        title = "Total Abundance of Vitamin B6 Biosynthesis Genes",
        x = NULL, y = "Abundance (CPM)"
      ) +
      theme_bw() +
      theme(legend.position = "none")

    ggsave("VitaminB6_Genes_Total.pdf", B6_Genes_Plot,
           width = 8, height = 6, device = cairo_pdf)

    # Statistical test
    cat("\n--- Kruskal-Wallis Test for B6 Genes ---\n")
    print(kruskal.test(B6_Genes_Total ~ Groups, data = b6_gene_data))

    write.csv(b6_gene_data, "VitaminB6_GeneAbundances.csv", row.names = FALSE)

  }
} else {
  cat("Gene family data (GeigerGenes) not found in workspace.\n")
}


#===============================================================================
# SECTION 4: CORRELATION ANALYSIS
#===============================================================================

cat("\n=== SECTION 4: B6 Producers vs Alpha Diversity ===\n\n")

if (exists("Diversity") && exists("species_data") && "B6_Producers_Total" %in% colnames(species_data)) {

  # Merge B6 producer data with diversity metrics
  correlation_data <- merge(
    Diversity[, c("SampleID", "Groups", "Shannon", "S.chao1")],
    species_data[, c("SampleID", "B6_Producers_Total")],
    by = "SampleID"
  )

  # Correlation test
  cor_shannon <- cor.test(correlation_data$Shannon, correlation_data$B6_Producers_Total,
                          method = "spearman")
  cor_chao1 <- cor.test(correlation_data$S.chao1, correlation_data$B6_Producers_Total,
                        method = "spearman")

  cat(sprintf("Correlation with Shannon: rho = %.3f, p = %.4f\n",
              cor_shannon$estimate, cor_shannon$p.value))
  cat(sprintf("Correlation with Chao1: rho = %.3f, p = %.4f\n",
              cor_chao1$estimate, cor_chao1$p.value))

  # Scatter plot
  B6_Correlation_Plot <- ggplot(correlation_data, aes(x = B6_Producers_Total, y = Shannon)) +
    geom_point(aes(color = Groups), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    scale_color_manual(values = FigCols) +
    labs(
      title = "B6 Producer Abundance vs Shannon Diversity",
      subtitle = sprintf("Spearman rho = %.3f, p = %.4f", cor_shannon$estimate, cor_shannon$p.value),
      x = "Total B6 Producer Abundance",
      y = "Shannon Diversity Index"
    ) +
    theme_bw()

  ggsave("VitaminB6_Correlation_Diversity.pdf", B6_Correlation_Plot,
         width = 8, height = 6, device = cairo_pdf)
}


#===============================================================================
# SECTION 5: SUMMARY OUTPUT
#===============================================================================

cat("\n\n========================================\n")
cat("Vitamin B6 Analysis Complete!\n")
cat("========================================\n\n")

cat("Output files generated:\n")
cat("  - VitaminB6_Pathway_Abundance.pdf\n")
cat("  - VitaminB6_Pathway_KruskalWallis.csv\n")
cat("  - VitaminB6_Pathway_Pairwise.csv\n")
cat("  - VitaminB6_DY_vs_DO.csv\n")
cat("  - VitaminB6_Producers_Total.pdf\n")
cat("  - VitaminB6_Producers_Individual.pdf\n")
cat("  - VitaminB6_Genes_Total.pdf\n")
cat("  - VitaminB6_GeneAbundances.csv\n")
cat("  - VitaminB6_Correlation_Diversity.pdf\n")
