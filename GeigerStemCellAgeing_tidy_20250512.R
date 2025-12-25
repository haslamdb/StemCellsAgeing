#-------------------------------------------------------------------------------
#
# Script: Stem Cell Aging Microbiome Analysis
#
# Author:  David Haslam, Selina Stahl, Hartmut Geiger
# Date:    July 14, 2025
#
# Description: This script analyzes metagenomic data from a study on hematopoietic
#              stem cell (HSC) aging. It performs a comprehensive analysis
#              including data preprocessing, diversity calculations, principal
#              component analysis (PCA), differential abundance testing, and
#              functional pathway analysis.
#
#-------------------------------------------------------------------------------

#===============================================================================
# SECTION 1: SETUP - PACKAGE MANAGEMENT
#===============================================================================

# List of required packages for the analysis
list.of.packages <- c(
  "reshape2", "stringdist", "stringr", "plyr", "vegan", "labdsv", "pvclust",
  "ggplot2", "gplots", "RColorBrewer", "Heatplus", "fossil", "ade4", "scales",
  "extrafont", "ggbiplot", "MASS", "ggthemes", "tidyverse", "pheatmap",
  "magrittr", "readxl", "robustbase", "cowplot", "sda", "locfdr",
  "FactoMineR", "factoextra", "dunn.test", "FSA", "NBZIMM", "DataCombine",
  "gtools", "heatmap3", "glmmTMB", "DescTools", "CompQuadForm", "dirmult",
  "ecodist", "GUniFrac", "lme4", "Matrix", "permute", "GMPR", "bbmle",
  "VennDiagram", "Cairo", "gridExtra", "ggpubr", "rstatix"
)

# Identify and install any missing packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

# Load all required packages
lapply(list.of.packages, require, character.only = TRUE)

#===============================================================================
# SECTION 2: GLOBAL PARAMETERS & CUSTOM FUNCTIONS
#===============================================================================

# --- Set Working Directory ---
# NOTE: Update this path to your project's directory. Using R Projects is recommended.
# setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
# setwd("~/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")


# --- Color Palettes & Figure Themes ---

# Define custom color palettes for plots
metaphlan.colors <- colorRampPalette(c("#000033", "#007FFF", "cyan", "red", "yellow"))
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(50)
diverging.colors <- colorRampPalette(c("#000033", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#fee090", "#fdae61", "#f46d43", "#d73027", "red"))

# Define colors for experimental groups
# Y: Young C57BL/6 (#A4DCFE - light blue)
# O: Old C57BL/6 (#074080 - dark blue)
# RAG1-/-: RAG1-/- untransplanted (#B6474B - red-brown)
# DY: Donor Young HSCs (#F9CC66 - yellow)
# DO: Donor Old HSCs (#FD8008 - orange)
FigCols <- c("#A4DCFE", "#074080", "#B6474B", "#F9CC66", "#FD8008")


# --- Custom ggplot Theme Functions ---

# Standard theme for plots
params <- function(x) {
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 22, color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14)
  )
}

# Theme for plots with angled x-axis text
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

# Theme for boxplots
paramsBox <- function(x) {
  theme(
    axis.text.x = element_text(size = 22, color = "black", family = "Arial"),
    axis.text.y = element_text(size = 16, angle = 0, color = "black", family = "Arial"),
    plot.title = element_text(size = 22, color = "black", family = "Arial"),
    axis.title.x = element_text(size = 20, family = "Arial"),
    axis.title.y = element_text(size = 22, family = "Arial"),
    legend.title = element_text(size = 18, family = "Arial"),
    legend.text = element_text(size = 14, family = "Arial"),
    legend.position = "none"
  )
}


# --- Analysis & Utility Functions ---

# Function to plot a linear regression with key statistics in the title
ggplotRegression <- function(fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(
      title = paste(
        "Adj R2 = ", signif(summary(fit)$adj.r.squared, 5),
        "Intercept =", signif(fit$coef[[1]], 5),
        " Slope =", signif(fit$coef[[2]], 5),
        " P =", signif(summary(fit)$coef[2, 4], 5)
      )
    )
}

# Partitioning Around Medoids (PAM) clustering function
pam.clustering <- function(x, k) {
  require(cluster)
  cluster <- as.vector(pam(as.dist(x), k, diss = TRUE)$clustering)
  return(cluster)
}

# Jensen-Shannon Divergence (JSD) distance calculator
dist.JSD <- function(inMatrix, pseudocount = 0.000001, ...) {
  KLD <- function(x, y) sum(x * log(x / y))
  JSD <- function(x, y) sqrt(0.5 * KLD(x, (x + y) / 2) + 0.5 * KLD(y, (x + y) / 2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix <- apply(inMatrix, 1:2, function(x) ifelse(x == 0, pseudocount, x))
  
  for (i in 1:matrixColSize) {
    for (j in 1:matrixColSize) {
      resultsMatrix[i, j] <- JSD(as.vector(inMatrix[, i]), as.vector(inMatrix[, j]))
    }
  }
  colnames(resultsMatrix) <- rownames(resultsMatrix) <- colnames
  resultsMatrix <- as.dist(resultsMatrix)
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix)
}

# Function to remove low-abundance features (noise)
noise.removal <- function(dataframe, percent = 0.001, top = NULL) {
  Matrix <- dataframe
  bigones <- rowSums(Matrix) * 100 / (sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones, ]
  print(percent)
  return(Matrix_1)
}

# Generalized log base 2 transformation
glog2 <- function(x) {
  (asinh(x) - log(2)) / log(2)
}

# Function to compute a robust mean using Huber's M-estimator
robustMean <- function(x) {
  huberM(x)$mu
}

# --- Shrinkage Discriminant Analysis (SDA) Effect Size Functions ---

# Effect size calculation for Background (C57BL vs. RAG1)
compute_ef_Background <- function(d, g, min_shrink = 0.3) {
  raw_scores <- sda.ranking(d, g, diagonal = FALSE, verbose = TRUE, fdr = FALSE)
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs) / freqs / length(g))
  
  ef <- raw_scores[, "cat.C57BL"] * m["C57BL"] - raw_scores[, "cat.RAG1"] * m["RAG1"]
  
  n0 <- sum(g_summaries$idx[, "C57BL"])
  n1 <- sum(g_summaries$idx[, "RAG1"])
  m_stats <- 1 / sqrt(1 / n0 + 1 / n1)
  stats <- m_stats * ef
  
  lfdr <- locfdr(stats, nulltype = 1)$fdr
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - lfdr)
  
  res <- tibble(
    Species = names(ef),
    ef = ef,
    ef_shrunk = ef_shrunk,
    stat = stats,
    lfdr = lfdr
  )
  return(res)
}

# Effect size calculation for Donor HSCs (Young vs. Old)
compute_ef_YoungOldHSCs <- function(d, g, min_shrink = 0.3) {
  # Pre-computation checks and data cleaning
  if (any(is.na(d))) {
    warning("Input data contains NA values. Imputing with zeros.")
    d[is.na(d)] <- 0
  }
  g <- factor(g)
  if (any(table(g) < 2)) {
    stop("Each group needs at least 2 samples for sda.ranking.")
  }
  
  raw_scores <- sda.ranking(d, g, diagonal = FALSE, verbose = TRUE, fdr = FALSE)
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs) / freqs / length(g))
  
  ef <- raw_scores[, "cat.DY"] * m["DY"] - raw_scores[, "cat.DO"] * m["DO"]
  
  n0 <- sum(g_summaries$idx[, "DY"])
  n1 <- sum(g_summaries$idx[, "DO"])
  m_stats <- 1 / sqrt(1 / n0 + 1 / n1)
  stats <- m_stats * ef
  
  lfdr <- locfdr(stats, nulltype = 1)$fdr
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - lfdr)
  
  res <- tibble(
    Species = names(ef),
    ef = ef,
    ef_shrunk = ef_shrunk,
    stat = stats,
    lfdr = lfdr
  )
  return(res)
}

# Effect size calculation for Recipient Age (Young vs. Old)
compute_ef_Recipient <- function(d, g, min_shrink = 0.3) {
  raw_scores <- sda.ranking(d, g, diagonal = FALSE, verbose = TRUE, fdr = FALSE)
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs) / freqs / length(g))
  
  ef <- raw_scores[, "cat.Y"] * m["Y"] - raw_scores[, "cat.O"] * m["O"]
  
  n0 <- sum(g_summaries$idx[, "Y"])
  n1 <- sum(g_summaries$idx[, "O"])
  m_stats <- 1 / sqrt(1 / n0 + 1 / n1)
  stats <- m_stats * ef
  
  lfdr <- locfdr(stats, nulltype = 1)$fdr
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - lfdr)
  
  res <- tibble(
    Species = names(ef),
    ef = ef,
    ef_shrunk = ef_shrunk,
    stat = stats,
    lfdr = lfdr
  )
  return(res)
}


#===============================================================================
# SECTION 3: METADATA LOADING AND PREPROCESSING
#===============================================================================

# Load metadata from CSV file
Metadata <- read.csv("GeigerSampleKeyRevised20220910_Corrected.csv", header = TRUE, stringsAsFactors = FALSE)
names(Metadata)[1] <- "SampleID"
Metadata <- Metadata[, 1:7]

# Define factor levels for categorical variables
Metadata$Background <- factor(Metadata$Background, levels = c("C57BL", "RAG1"))
Metadata$Transplant <- factor(Metadata$Transplant, levels = c("None", "DY", "DO", "Casin_treated_DO", "RAG1"))
Metadata$Recipient <- factor(Metadata$Recipient, levels = c("Young", "Old"))
Metadata$ImmuneSystem <- factor(Metadata$ImmuneSystem, levels = c("Young", "Rejuvenated", "Old"))
Metadata$Cage <- factor(Metadata$Cage)
Metadata$Experiment <- factor(Metadata$Experiment)

# Create combined group identifiers
Metadata$Groups <- paste(Metadata$Recipient, Metadata$Background, Metadata$Transplant, Metadata$ImmuneSystem, sep = "-")
Metadata$CageGroups <- paste(Metadata$Groups, Metadata$Cage, sep = "-")

# Simplify and format group names for clarity in plots
Metadata$Groups <- gsub("Young-RAG1-Casin_treated_DO-Rejuvenated", "RAG1-Rejuvenated_HSCs", Metadata$Groups)
Metadata$Groups <- gsub("Young-RAG1-DY-Young", "RAG1-DY", Metadata$Groups)
Metadata$Groups <- gsub("Young-RAG1-DO-Old", "RAG1-DO", Metadata$Groups)
Metadata$Groups <- gsub("Young-RAG1-None-Young", "Young_RAG1-Untransplanted", Metadata$Groups)
Metadata$Groups <- gsub("Old-C57BL-None-Old", "Old_C57BL-Untransplanted", Metadata$Groups)
Metadata$Groups <- gsub("Young-C57BL-None-Young", "Young_C57BL-Untransplanted", Metadata$Groups)

# Further simplify group names for final figures
Metadata$Groups <- gsub("RAG1-Rejuvenated_HSCs", "Rejuv", Metadata$Groups)
Metadata$Groups <- gsub("RAG1-DY", "DY", Metadata$Groups)
Metadata$Groups <- gsub("RAG1-DO", "DO", Metadata$Groups)
Metadata$Groups <- gsub("Young_RAG1-Untransplanted", "RAG1-/-", Metadata$Groups)
Metadata$Groups <- gsub("Old_C57BL-Untransplanted", "O", Metadata$Groups)
Metadata$Groups <- gsub("Young_C57BL-Untransplanted", "Y", Metadata$Groups)

# Remove unused experimental groups
Metadata <- subset(Metadata, !Metadata$Groups %in% c("Young-RAG1-RAG1-Young", "Rejuv"))

# Set final factor levels for plotting order
Metadata$Groups <- factor(Metadata$Groups, levels = c("Y", "O", "RAG1-/-", "DY", "DO"))
Metadata$GroupsExp <- paste(Metadata$Groups, Metadata$Experiment, sep = "-")

# Store the list of sample IDs to be used in the analysis
GeigerSamples <- Metadata$SampleID


#===============================================================================
# SECTION 4: KRAKEN2 SPECIES DATA - LOADING AND MERGING
#===============================================================================

# NOTE: Update this path to where your Kraken2 output files are located.
# setwd("C:/Users/dbhas/OneDrive/Documents/Alignments/KrakenAlignments/Kraken2")
# setwd("~/Documents/Alignments/KrakenAlignments/Kraken2")

# Get a list of all species abundance files
AllKrakenFiles <- list.files()
SpeciesFileList <- grep("_species_abundance.txt", AllKrakenFiles, value = TRUE)
FileList <- gsub("_species_abundance.txt", "", SpeciesFileList)

# Filter for files corresponding to the samples in the metadata
NewSpeciesFileList <- subset(FileList, FileList %in% GeigerSamples)

# --- Loop to load and merge all species files into a single dataframe ---

# Initialize with the first file
first_file_path <- paste0(NewSpeciesFileList[1], "_species_abundance.txt")
NewSpeciesNR <- read.csv(first_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(NewSpeciesNR) <- c("Species", "TaxonomyID", "TaxonomyRank", "Count")
NewSpeciesNR <- NewSpeciesNR[, c("Species", "Count")]

# Clean species names
NewSpeciesNR$Species <- gsub("[\\s_\\[\\]\\-/]|X,", ".", NewSpeciesNR$Species)
NewSpeciesNR$Species <- gsub("\\.\\.", ".", NewSpeciesNR$Species)
NewSpeciesNR <- subset(NewSpeciesNR, !duplicated(Species))
names(NewSpeciesNR)[2] <- NewSpeciesFileList[1]

# Loop through the rest of the files and merge them
for (i in 2:length(NewSpeciesFileList)) {
  file_path <- paste0(NewSpeciesFileList[i], "_species_abundance.txt")
  x <- read.csv(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  x <- x[, c(1, 4)]
  names(x) <- c("Species", NewSpeciesFileList[i])
  
  # Clean species names
  x$Species <- gsub("[\\s_\\[\\]\\-/]|X,", ".", x$Species)
  x$Species <- gsub("\\.\\.", ".", x$Species)
  x <- subset(x, !duplicated(Species))
  
  NewSpeciesNR <- merge(x, NewSpeciesNR, by = "Species", all = TRUE)
}

# Convert to a matrix with species as rows and samples as columns
row.names(NewSpeciesNR) <- NewSpeciesNR$Species
NewSpeciesNR$Species <- NULL
NewSpeciesNR[is.na(NewSpeciesNR)] <- 0

# Clean up intermediate objects
rm(list = FileList)

# Navigate back to the main project directory
# setwd("~/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")


#===============================================================================
# SECTION 5: SPECIES DATA - CLEANING AND NORMALIZATION
#===============================================================================

# Transpose so samples are rows and species are columns
GeigerSpecies <- as.data.frame(t(NewSpeciesNR))

# Remove potential contaminants (if necessary, currently commented out)
# SaliniCol <- grep("Salinibacter", colnames(GeigerSpecies))
# GeigerSpecies <- GeigerSpecies[, -SaliniCol]
# HumanRow <- grep("Homo.sapiens", colnames(GeigerSpecies))
# GeigerSpecies <- GeigerSpecies[, -HumanRow]
# MouseRow <- grep("Mus.musculus", colnames(GeigerSpecies))
# GeigerSpecies <- GeigerSpecies[, -MouseRow]

# Filter species present in less than 5% of samples
TenPercentCutoff <- floor(nrow(GeigerSpecies) / 20)
NonZeroCounts <- sapply(GeigerSpecies, function(x) sum(x > 0))
TenPercentNotZero <- which(NonZeroCounts >= TenPercentCutoff)
GeigerSpeciesNR <- GeigerSpecies[, TenPercentNotZero]

# Remove low-abundance species (less than 0.001% of total reads)
GeigerSpeciesNR <- as.data.frame(t(noise.removal(t(GeigerSpeciesNR), 0.001)))

# Remove samples with very low read counts (less than 750,000)
LowSamples <- which(rowSums(GeigerSpeciesNR) <= 750000)
if (length(LowSamples) > 0) {
  GeigerSpeciesNR <- GeigerSpeciesNR[-LowSamples, ]
}

# Rarefy samples to an even depth (750,000 reads) for normalization
GeigerSpeciesNR <- data.frame(rrarefy(GeigerSpeciesNR, 750000))
GeigerSpeciesNR$SampleID <- row.names(GeigerSpeciesNR)

# Backup the processed data
BackupGeigerSpecies <- GeigerSpeciesNR

# Merge with metadata to create the final analysis-ready dataframe
GeigerSpeciesNR <- merge(Metadata, GeigerSpeciesNR, by = "SampleID", all.x = TRUE)
row.names(GeigerSpeciesNR) <- GeigerSpeciesNR$SampleID

# Separate species matrix for downstream analysis
Species <- as.data.frame(t(GeigerSpeciesNR[, c(11:ncol(GeigerSpeciesNR))]))


#===============================================================================
# SECTION 6: ALPHA DIVERSITY ANALYSIS
#===============================================================================

# Calculate various alpha diversity metrics
H <- diversity(t(Species))
simpson <- diversity(t(Species), "simpson")
shannon <- diversity(t(Species), "shannon")
invsimp <- diversity(t(Species), "inv")
alpha <- fisher.alpha(t(Species))
S <- specnumber(t(Species)) # Species richness
J <- H / log(S) # Pielou's evenness

# Combine metrics into a single dataframe
Diversity <- data.frame(
  Simpson = simpson,
  Shannon = shannon,
  InvSimpson = invsimp,
  SpeciesNo = S,
  Evenness = J,
  SampleID = row.names(GeigerSpeciesNR)
)

# Merge diversity data with metadata
Diversity <- merge(Metadata, Diversity, by = "SampleID", all.x = TRUE)

# --- Plot Shannon Diversity ---
Diversity_Groups <- ggplot(Diversity, aes(x = Groups, y = Shannon)) +
  geom_boxplot(lwd = 1, aes(color = factor(Groups)), fill = NA, outlier.size = 3) +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 8) +
  scale_colour_manual(values = FigCols) +
  geom_point(size = 4, aes(color = factor(Groups))) +
  xlab(NULL) +
  ylab("Shannon Diversity Index") +
  theme_bw() +
  paramsBox()

ggsave(filename = "Diversity_NewGroups.pdf", plot = Diversity_Groups, width = 10, height = 10)

# Facet by experiment
Diversity_NewGroups_Experiment <- ggplot(Diversity, aes(x = Groups, y = Shannon)) +
  geom_boxplot(lwd = 1, aes(color = factor(Groups)), fill = NA, outlier.size = 3) +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 8) +
  scale_colour_manual(values = FigCols) +
  geom_point(size = 4, aes(color = factor(Groups))) +
  xlab(NULL) +
  ylab("Shannon Diversity Index") +
  facet_grid(. ~ Experiment) +
  theme_bw() +
  paramsAngled()

ggsave(filename = "Diversity_NewGroups_Experiment.pdf", plot = Diversity_NewGroups_Experiment, width = 10, height = 10)

# Statistical test for diversity
pairwise.wilcox.test(Diversity$Shannon, Diversity$Groups, p.adjust.method = "bonferroni")


#===============================================================================
# SECTION 7: BETA DIVERSITY & PRINCIPAL COMPONENT ANALYSIS (PCA)
#===============================================================================

# --- Multi-Response Permutation Procedure (MRPP) ---
# Tests for significant differences between predefined groups
clean_data <- na.omit(GeigerSpeciesNR[, 11:ncol(GeigerSpeciesNR)])
clean_groups <- GeigerSpeciesNR[!is.na(rowSums(GeigerSpeciesNR[, 11:ncol(GeigerSpeciesNR)])), "Groups"]
Groups.mrpp <- mrpp(clean_data, clean_groups, distance = "bray")
print(Groups.mrpp) # p = 0.001 indicates significant group separation

# --- PCA of Experimental Groups (Figure 1A) ---
# Subset data for Experiment 1: Y, RAG1-/-, DY
PCASpecies_Fig1A <- subset(GeigerSpeciesNR, Experiment == 1 & Groups %in% c("DY", "RAG1-/-", "Y"))

metadata_Fig1A <- PCASpecies_Fig1A[, 1:10]
cts_Fig1A <- as.matrix(PCASpecies_Fig1A[, -(1:10)])
cts_l2_Fig1A <- glog2(cts_Fig1A)
grps_Fig1A <- metadata_Fig1A$Groups

# Perform PCA
pca_Fig1A <- PCA(cts_l2_Fig1A, scale.unit = FALSE, ncp = 5, graph = FALSE)

# Plot PCA
Fig1A_PCA_Plot <- fviz_pca_ind(pca_Fig1A,
                               geom.ind = "point", pointsize = 2,
                               title = "PCA of Microbial Species (Experiment 1)",
                               subtitle = "Groups: Young, RAG1-/-, Donor Young HSC",
                               col.ind = grps_Fig1A,
                               addEllipses = TRUE, ellipse.type = "confidence",
                               legend.title = "Mouse Group",
                               palette = FigCols[c(1, 3, 4)],
                               mean.point = TRUE
)

pdf("Figure1A_PCA.pdf", width = 8, height = 7)
print(Fig1A_PCA_Plot)
dev.off()


# --- PCA of Young vs. Old Control Mice (Figure 1D) ---
# Subset data for Experiment 1: Y, O
PCASpecies_Fig1D <- subset(GeigerSpeciesNR, Experiment == 1 & Groups %in% c("Y", "O"))

metadata_Fig1D <- PCASpecies_Fig1D[, 1:10]
cts_Fig1D <- as.matrix(PCASpecies_Fig1D[, -(1:10)])
cts_l2_Fig1D <- glog2(cts_Fig1D)
grps_Fig1D <- metadata_Fig1D$Groups

# Perform PCA
pca_Fig1D <- PCA(cts_l2_Fig1D, scale.unit = FALSE, ncp = 5, graph = FALSE)

# Plot PCA
Fig1D_PCA_Plot <- fviz_pca_ind(pca_Fig1D,
                               geom.ind = "point", pointsize = 2,
                               title = "PCA of Microbial Species (Experiment 1)",
                               subtitle = "Groups: Young vs. Old Controls",
                               col.ind = grps_Fig1D,
                               addEllipses = TRUE, ellipse.type = "confidence",
                               legend.title = "Mouse Group",
                               palette = FigCols[c(1, 2)],
                               mean.point = TRUE
)

pdf("Figure1D_PCA.pdf", width = 8, height = 7)
print(Fig1D_PCA_Plot)
dev.off()

# --- PCA of Transplanted Mice by Experiment (Supp. Figure 4) ---
# Goal: Plot each experiment's PCA on a consistent set of axes for comparison.

# 1. Run a single PCA on ALL transplanted samples from all experiments
all_transplant_indices <- which(metadata$Transplant != "None")
all_transplant_cts <- cts_l2[all_transplant_indices, ]
all_transplant_metadata <- metadata[all_transplant_indices, ]
all_transplant_pca <- PCA(all_transplant_cts, scale.unit = FALSE, ncp = 5, graph = FALSE)

# 2. Create a function to plot a single experiment using the master PCA coordinates
plot_experiment <- function(exp_num, master_pca, all_meta, cols) {
  p <- fviz_pca_ind(master_pca,
                    geom.ind = "point", pointsize = 3,
                    title = paste("Experiment", exp_num),
                    habillage = all_meta$Transplant,
                    addEllipses = TRUE, ellipse.type = "confidence",
                    legend.title = "Transplant",
                    palette = cols,
                    select.ind = list(name = all_meta$SampleID[all_meta$Experiment == exp_num])
  ) + theme_minimal()
  return(p)
}

# 3. Generate plots for each experiment
p1 <- plot_experiment(1, all_transplant_pca, all_transplant_metadata, FigCols[c(4, 5)])
p2 <- plot_experiment(2, all_transplant_pca, all_transplant_metadata, FigCols[c(4, 5)])
p3 <- plot_experiment(3, all_transplant_pca, all_transplant_metadata, FigCols[c(4, 5)])

# 4. Arrange the plots in a grid
combined_grid <- grid.arrange(p1, p2, p3, ncol = 3,
                              top = "PCA of Transplanted Mice by Experiment (Consistent Axes)")

pdf("SuppFig4_Combined_Transplant_PCA_Grid.pdf", width = 15, height = 6)
print(combined_grid)
dev.off()


#===============================================================================
# SECTION 8: BRAY-CURTIS DISTANCE ANALYSIS (FIGURE 1B & 1C)
#===============================================================================

# Calculate Bray-Curtis distance matrix for Experiment 1 samples
cts_Fig1A_dist <- vegdist(cts_l2_Fig1A, method = "bray")
bray_matrix <- as.matrix(cts_Fig1A_dist)

# Melt the matrix into a long format for plotting
bray_long <- bray_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "Distance") %>%
  # Remove self-comparisons and duplicates
  filter(Sample1 < Sample2) %>%
  # Add group information
  left_join(metadata_Fig1A %>% select(SampleID, Group1 = Groups), by = c("Sample1" = "SampleID")) %>%
  left_join(metadata_Fig1A %>% select(SampleID, Group2 = Groups), by = c("Sample2" = "SampleID"))

# Categorize pairwise comparisons
categorize <- function(g1, g2) {
  if (g1 == g2) return(NA) # Skip intra-group
  pair <- sort(c(as.character(g1), as.character(g2)))
  return(paste(pair, collapse = "-"))
}
bray_long$Category <- mapply(categorize, bray_long$Group1, bray_long$Group2)

# Filter for comparisons of interest and set factor order for plotting
filtered_data <- bray_long %>%
  filter(Category %in% c("DY-Y", "RAG1-/- -Y", "DY-RAG1-/-")) %>%
  mutate(Category = factor(Category, levels = c("DY-Y", "RAG1-/- -Y", "DY-RAG1-/-")))

# --- Create Boxplot (Figure 1B) ---
p_bray <- ggplot(filtered_data, aes(x = Category, y = Distance)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, aes(color = Category)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "red") +
  scale_color_manual(values = c("DY-Y" = "#7CAE00", "RAG1-/- -Y" = "#00BFC4", "DY-RAG1-/-" = "#F8766D")) +
  labs(
    title = "Pairwise Bray-Curtis Distances Between Groups",
    x = "Comparison",
    y = "Bray-Curtis Distance"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# --- Statistical Analysis (Figure 1C) ---
# Kruskal-Wallis test (non-parametric ANOVA)
kruskal_result <- kruskal.test(Distance ~ Category, data = filtered_data)
print(kruskal_result)

# Dunn's post-hoc test for pairwise comparisons
dunn_result <- filtered_data %>%
  dunn_test(Distance ~ Category, p.adjust.method = "bonferroni")
print(dunn_result)
write.csv(dunn_result, file = "Figure1C_Dunn_Test_Results.csv", row.names = FALSE)

# Add stats to the plot
p_bray_with_stats <- p_bray +
  stat_pvalue_manual(dunn_result, label = "p.adj.signif", y.position = 1.05)

pdf("Figure1B_C_BrayCurtisDistances.pdf", width = 8, height = 7)
print(p_bray_with_stats)
dev.off()


#===============================================================================
# SECTION 9: DIFFERENTIAL ABUNDANCE ANALYSIS
#===============================================================================

# --- Comparison 1: Young (Y) vs. Old (O) Untransplanted Mice ---
YoungOldSamples <- subset(GeigerSpeciesNR, Groups %in% c("Y", "O"))
metadata_YO <- YoungOldSamples[, 1:10]
cts_YO <- as.matrix(YoungOldSamples[, -(1:10)])
cts_l2_YO <- glog2(cts_YO)

# Calculate effect sizes using SDA
res_Recipient <- compute_ef_Recipient(d = cts_l2_YO, g = metadata_YO$Groups, min_shrink = 0.3)

# Filter for significant species and plot effect sizes (Supp. Figure 1B)
Recipient_Effect_Plot <- res_Recipient %>%
  filter(abs(ef_shrunk) > 1.25) %>%
  mutate(
    HigherIn = if_else(ef_shrunk > 0, "Young", "Old"),
    Species = reorder(Species, ef_shrunk)
  ) %>%
  ggplot(aes(x = Species, y = ef_shrunk, fill = HigherIn)) +
  geom_col(color = "black") +
  coord_flip() +
  scale_fill_manual(values = FigCols[c(1, 2)]) +
  labs(
    y = "Shrunken Effect Size", x = NULL,
    title = "Differentially Abundant Species: Young vs. Old Controls",
    fill = "More Abundant In:"
  ) +
  theme_bw()

ggsave(filename = "SuppFigure1B_Species_Effect_Young_vs_Old.pdf", plot = Recipient_Effect_Plot, width = 10, height = 8, device = cairo_pdf)

# --- Comparison 2: Donor Young (DY) vs. Donor Old (DO) HSCs ---
YoungOldHSC <- subset(GeigerSpeciesNR, Groups %in% c("DY", "DO"))
metadata_HSC <- YoungOldHSC[, 1:10]
cts_HSC <- as.matrix(YoungOldHSC[, -(1:10)])
cts_l2_HSC <- glog2(cts_HSC)

# Calculate effect sizes using SDA
res_YoungOldHSC <- compute_ef_YoungOldHSCs(d = cts_l2_HSC, g = metadata_HSC$Groups, min_shrink = 0.3)

# Filter for significant species and plot effect sizes (Figure 4)
YoungOldHSC_Effect_Plot <- res_YoungOldHSC %>%
  filter(abs(ef_shrunk) > 1) %>%
  mutate(
    HigherIn = if_else(ef_shrunk > 0, "DY", "DO"),
    Species = reorder(Species, ef_shrunk)
  ) %>%
  ggplot(aes(x = Species, y = ef_shrunk, fill = HigherIn)) +
  geom_col(color = "black") +
  coord_flip() +
  scale_fill_manual(values = FigCols[c(4, 5)]) +
  labs(
    y = "Shrunken Effect Size", x = NULL,
    title = "Differentially Abundant Species: Donor Young vs. Donor Old HSCs",
    fill = "More Abundant In:"
  ) +
  theme_bw() +
  paramsBox()

ggsave(filename = "Figure4_Species_Effect_DY_vs_DO.pdf", plot = YoungOldHSC_Effect_Plot, width = 16, height = 10, device = cairo_pdf)


#===============================================================================
# SECTION 10: PHYLUM-LEVEL ANALYSIS (F/B RATIO)
#===============================================================================

# --- Load Phylum-level Kraken2 Data ---
# NOTE: Update this path to where your Kraken2 output files are located.
# setwd("C:/Users/dbhas/OneDrive/Documents/Alignments/KrakenAlignments/Kraken2")
# setwd("~/Documents/Alignments/KrakenAlignments/Kraken2")

AllKrakenFiles_phylum <- list.files()
phylumFileList <- grep("_phylum_abundance.txt", AllKrakenFiles_phylum, value = TRUE)
FileList_phylum <- gsub("_phylum_abundance.txt", "", phylumFileList)

NewphylumFileList <- subset(FileList_phylum, FileList_phylum %in% GeigerSamples)

# Initialize with the first file
first_phylum_path <- paste0(NewphylumFileList[1], "_phylum_abundance.txt")
NewphylumNR <- read.csv(first_phylum_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
NewphylumNR <- NewphylumNR[, c(1, 4)]
names(NewphylumNR) <- c("phylum", NewphylumFileList[1])

# Clean phylum names
NewphylumNR$phylum <- gsub("[\\s_\\[\\]\\-/]|X,", ".", NewphylumNR$phylum)
NewphylumNR$phylum <- gsub("\\.\\.", ".", NewphylumNR$phylum)
NewphylumNR <- subset(NewphylumNR, !duplicated(phylum))

# Loop through the rest of the files and merge
for (i in 2:length(NewphylumFileList)) {
  file_path <- paste0(NewphylumFileList[i], "_phylum_abundance.txt")
  x <- read.csv(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  x <- x[, c(1, 4)]
  names(x) <- c("phylum", NewphylumFileList[i])
  x$phylum <- gsub("[\\s_\\[\\]\\-/]|X,", ".", x$phylum)
  x$phylum <- gsub("\\.\\.", ".", x$phylum)
  x <- subset(x, !duplicated(phylum))
  NewphylumNR <- merge(x, NewphylumNR, by = "phylum", all = TRUE)
}

# Convert to a matrix
row.names(NewphylumNR) <- NewphylumNR$phylum
NewphylumNR$phylum <- NULL
NewphylumNR[is.na(NewphylumNR)] <- 0

# --- Calculate and Plot Firmicutes/Bacteroidetes Ratio (Supp. Figure 1A) ---
PhylumTableT <- as.data.frame(t(NewphylumNR))
PhylumTableT$SampleID <- rownames(PhylumTableT)

# Subset for Young vs. Old untransplanted mice
PhylumYOTable <- PhylumTableT[PhylumTableT$SampleID %in% subset(Metadata, Groups %in% c("Y", "O"))$SampleID, ]
PhylumYOTable <- merge(Metadata[, c("SampleID", "Groups")], PhylumYOTable, by = "SampleID")
PhylumYOTable$Groups <- factor(PhylumYOTable$Groups, levels = c("Y", "O"))

# Calculate F/B Ratio
FirmicutesCol <- grep("Firmicutes", colnames(PhylumYOTable), ignore.case = TRUE)
BacteroidetesCol <- grep("Bacteroidetes", colnames(PhylumYOTable), ignore.case = TRUE)
PhylumYOTable$Firmicutes_Bacteroidetes <- as.numeric(PhylumYOTable[[FirmicutesCol]] / PhylumYOTable[[BacteroidetesCol]])

# Plot the F/B Ratio
FB_Ratio_Boxplot <- ggplot(PhylumYOTable, aes(x = Groups, y = Firmicutes_Bacteroidetes)) +
  geom_boxplot(lwd = 1, aes(color = Groups), fill = NA, outlier.size = 3) +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 8) +
  geom_point(size = 4, aes(color = Groups)) +
  scale_color_manual(values = FigCols[c(1, 2)]) +
  labs(y = "Firmicutes : Bacteroidetes Ratio", x = NULL) +
  theme_bw() +
  paramsBox()

# Save the plot
ggsave(filename = "SuppFigure1A_FB_Ratio.pdf", plot = FB_Ratio_Boxplot, width = 8, height = 12, device = cairo_pdf)

# Perform Wilcoxon test
stat_test_fb <- PhylumYOTable %>%
  wilcox_test(Firmicutes_Bacteroidetes ~ Groups) %>%
  add_significance()
print(stat_test_fb) # p = 0.15 (ns)


#===============================================================================
# SECTION 11: FUNCTIONAL ANALYSIS (HUMAnN3)
#===============================================================================

# --- 11.1: Pathway Analysis - Data Loading and Preprocessing ---
# NOTE: Update this path to your HUMAnN3 output directory
# setwd("~/path/to/humann3_output/")

# Read in non-stratified Pathway table
PathwayTable <- read.csv("GeigerFiles_pathabundance20211018-cpm_unstratified.tsv", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
names(PathwayTable)[1] <- "Pathway"
colnames(PathwayTable) <- gsub(".paired_Abundance.CPM", "", colnames(PathwayTable))
PathwayTable <- PathwayTable[-c(1:2), ] # Remove header lines from HUMAnN3 table

# Clean pathway names and sample names
row.names(PathwayTable) <- PathwayTable$Pathway
PathwayTable$Pathway <- NULL
colnames(PathwayTable) <- gsub("X", "", colnames(PathwayTable))

# Transpose and merge with metadata
GeigerPathwayTable_t <- t(PathwayTable)
GeigerPathwayTabledf <- as.data.frame(GeigerPathwayTable_t)
GeigerPathwayTabledf$SampleID <- row.names(GeigerPathwayTabledf)

GeigerPathways <- merge(Metadata, GeigerPathwayTabledf, by = "SampleID", all.y = TRUE)
GeigerPathways <- subset(GeigerPathways, !is.na(Groups))
row.names(GeigerPathways) <- GeigerPathways$SampleID

# Clean column names for R compatibility
valid_colnames <- make.names(colnames(GeigerPathways), unique = TRUE)
colnames(GeigerPathways) <- valid_colnames


# --- 11.2: Pathway Analysis - PCA ---
pathway_metadata <- GeigerPathways[, 1:10]
pathway_cts <- as.matrix(GeigerPathways[, -(1:10)])
pathway_cts[is.na(pathway_cts)] <- 0 # Impute NA with 0
pathway_cts_l2 <- glog2(pathway_cts)
pathway_grps <- factor(pathway_metadata$Groups)

# Perform PCA
pathway_pca <- PCA(pathway_cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)

# Plot PCA
PathwaysPCA_Plot <- fviz_pca_ind(pathway_pca,
                                 geom.ind = "point",
                                 title = "PCA of Functional Pathways",
                                 col.ind = pathway_grps,
                                 addEllipses = TRUE, ellipse.type = "confidence",
                                 legend.title = "Mouse Group",
                                 palette = FigCols,
                                 mean.point = TRUE
)

pdf("PathwaysPCA.pdf", width = 8, height = 7)
print(PathwaysPCA_Plot)
dev.off()


# --- 11.3: Pathway Analysis - Differential Abundance ---

# Function to run pairwise Wilcoxon tests across all pathways
run_wilcox_test <- function(data, groups_col, start_col) {
  results <- list()
  for (i in start_col:ncol(data)) {
    formula <- as.formula(paste0("`", colnames(data)[i], "` ~ ", groups_col))
    if (length(unique(data[[groups_col]])) == 2) {
      test_res <- wilcox.test(formula, data = data)
      p_val <- test_res$p.value
    } else {
      # Fallback for more than 2 groups, though pairwise is intended
      test_res <- kruskal.test(formula, data = data)
      p_val <- test_res$p.value
    }
    results[[colnames(data)[i]]] <- p_val
  }
  
  WilcoxTable <- data.frame(
    Pathway = names(results),
    Unadjusted_p = unlist(results)
  )
  WilcoxTable$FDR <- p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
  WilcoxTable <- WilcoxTable[order(WilcoxTable$Unadjusted_p), ]
  return(WilcoxTable)
}

# Comparison 1: Donor Young (DY) vs. Donor Old (DO) HSCs
GeigerPathwaysYO <- subset(GeigerPathways, Groups %in% c("DY", "DO"))
WilcoxTable_DY_DO <- run_wilcox_test(GeigerPathwaysYO, "Groups", 11)
write.csv(WilcoxTable_DY_DO, file = "SignificantGeigerPathways_DY_vs_DO.csv", row.names = FALSE)

# Comparison 2: Young RAG1 vs. Young C57BL
GeigerPathwaysRC <- subset(GeigerPathways, Groups %in% c("RAG1-/-", "Y"))
WilcoxTable_RAG_C57 <- run_wilcox_test(GeigerPathwaysRC, "Groups", 11)
write.csv(WilcoxTable_RAG_C57, file = "SignificantGeigerPathways_RAG_vs_C57BL.csv", row.names = FALSE)

# Comparison 3: Young (Y) vs. Old (O) C57BL
GeigerPathwaysYOWT <- subset(GeigerPathways, Groups %in% c("Y", "O"))
WilcoxTable_Y_O <- run_wilcox_test(GeigerPathwaysYOWT, "Groups", 11)
write.csv(WilcoxTable_Y_O, file = "SignificantGeigerPathways_Y_vs_O.csv", row.names = FALSE)


# --- 11.4: Pathway Analysis - Overlap and Visualization ---
sigPathways_Y_O <- subset(WilcoxTable_Y_O, Unadjusted_p < 0.05)
sigPathways_DY_DO <- subset(WilcoxTable_DY_DO, Unadjusted_p < 0.05)

OverlappingYOPathways <- intersect(sigPathways_Y_O$Pathway, sigPathways_DY_DO$Pathway)

# Venn Diagram
PathwayVennDiagram <- draw.pairwise.venn(
  area1 = nrow(sigPathways_Y_O),
  area2 = nrow(sigPathways_DY_DO),
  cross.area = length(OverlappingYOPathways),
  category = c("Young vs Old Mice", "Young vs Old HSCs"),
  lty = "blank",
  fill = c("lightblue", "pink"),
  cat.pos = c(0, 0),
  fontfamily = "Helvetica",
  cat.fontfamily = "Helvetica"
)

ggsave(filename = "PathwayVennDiagram.pdf", plot = PathwayVennDiagram, width = 6, height = 6)

# --- 11.5: Gene Family Analysis - Data Loading and Preprocessing ---
# NOTE: This follows the same logic as pathway analysis
GeneTable <- read.csv("GeigerFiles_genefamilies20211018-cpm_unstratified.tsv", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
names(GeneTable)[1] <- "Gene"
colnames(GeneTable) <- gsub(".paired_Abundance.CPM", "", colnames(GeneTable))
GeneTable <- GeneTable[-c(1:2), ]

# Filter low-abundance genes to reduce complexity
GeneTable[, -1] <- lapply(GeneTable[, -1], as.numeric)
GeneTable$Total <- rowSums(GeneTable[, -1])
GeneTableSubset <- subset(GeneTable, Total > 500) # Keep genes with total CPM > 500
row.names(GeneTableSubset) <- GeneTableSubset$Gene
GeneTableSubset$Gene <- NULL
GeneTableSubset$Total <- NULL

# Transpose and merge with metadata
GeigerGeneTable_t <- t(GeneTableSubset)
GeigerGeneTabledf <- as.data.frame(GeigerGeneTable_t)
GeigerGeneTabledf$SampleID <- row.names(GeigerGeneTabledf)

GeigerGenes <- merge(Metadata, GeigerGeneTabledf, by = "SampleID", all.y = TRUE)
GeigerGenes <- subset(GeigerGenes, !is.na(Groups))
row.names(GeigerGenes) <- GeigerGenes$SampleID
colnames(GeigerGenes) <- make.names(colnames(GeigerGenes), unique = TRUE)

# --- 11.6: Gene Family Analysis - PCA ---
gene_metadata <- GeigerGenes[, 1:10]
gene_cts <- as.matrix(GeigerGenes[, -(1:10)])
gene_cts[is.na(gene_cts)] <- 0
gene_cts_l2 <- glog2(gene_cts)
gene_grps <- factor(gene_metadata$Groups)

# Perform PCA
gene_pca <- PCA(gene_cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)

# Plot PCA
GenesPCA_Plot <- fviz_pca_ind(gene_pca,
                              geom.ind = "point",
                              title = "PCA of Gene Families",
                              col.ind = gene_grps,
                              addEllipses = TRUE, ellipse.type = "confidence",
                              legend.title = "Mouse Group",
                              palette = FigCols,
                              mean.point = TRUE
)

pdf("GenesPCA.pdf", width = 8, height = 7)
print(GenesPCA_Plot)
dev.off()

# --- 11.7: Gene Family Analysis - Differential Abundance & Heatmaps ---
# The same differential abundance testing strategy (e.g., run_wilcox_test)
# and heatmap generation can be applied to the `GeigerGenes` dataframe as was
# done for pathways. For brevity, this is left as an exercise for the user,
# following the patterns established in section 11.3.

# Example: calculate log-fold change for a comparison
calculate_lfc <- function(data, group_col, groups_to_compare, start_col) {
  data_subset <- subset(data, get(group_col) %in% groups_to_compare)
  mean_abund <- aggregate(.~get(group_col), data = data_subset[, c(group_col, colnames(data)[start_col:ncol(data)])], FUN = mean)
  rownames(mean_abund) <- mean_abund[,1]
  mean_abund <- mean_abund[,-1]
  
  lfc_table <- as.data.frame(t(mean_abund))
  lfc_table$logratio <- foldchange2logratio(foldchange(lfc_table[,1], lfc_table[,2]))
  lfc_table$Gene <- rownames(lfc_table)
  return(lfc_table)
}

# Generate heatmap for significant genes in Y vs O
lfc_genes_Y_O <- calculate_lfc(GeigerGenes, "Groups", c("Y", "O"), 11)
sig_genes_Y_O <- run_wilcox_test(subset(GeigerGenes, Groups %in% c("Y", "O")), "Groups", 11)
lfc_genes_Y_O_sig <- merge(lfc_genes_Y_O, sig_genes_Y_O, by.x="Gene", by.y="Pathway") # Note: col name is "Pathway" from function
lfc_genes_Y_O_sig <- subset(lfc_genes_Y_O_sig, FDR < 0.2 & abs(logratio) > 1)
lfc_genes_Y_O_sig <- lfc_genes_Y_O_sig[order(lfc_genes_Y_O_sig$logratio),]

if(nrow(lfc_genes_Y_O_sig) > 1){
  matrixtable <- lfc_genes_Y_O_sig[, c("logratio"), drop=FALSE]
  rownames(matrixtable) <- lfc_genes_Y_O_sig$Gene
  
  pdf("Heatmap_SigGenes_Y_vs_O.pdf", width=10, height=15)
  heatmap.2(as.matrix(matrixtable), Rowv = FALSE, Colv = FALSE, col = diverging.colors(100), trace = "none",
            margins=c(5,25), main = "Significant Genes (Y vs O)", xlab="Log-Ratio")
  dev.off()
}


#===============================================================================
# SECTION 12: SAVE WORKSPACE
#===============================================================================

# Save the entire R workspace environment for later use
save.image(file = "Geiger_StemCell_Analysis_Workspace.RData")

cat("Analysis complete. Workspace saved to Geiger_StemCell_Analysis_Workspace.RData\n")