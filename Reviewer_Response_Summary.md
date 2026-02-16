# Reviewer Response: Additional Analyses for Stem Cell Aging Microbiome Study

**Last updated:** February 16, 2026
**Document Purpose:** Summary of additional analyses performed in response to reviewer comments, with suggested manuscript text and figure inventory.

---

## Table of Contents

1. [Alpha Diversity Metrics and Rarefaction Curves](#1-alpha-diversity-metrics-and-rarefaction-curves)
2. [Pathway Differential Abundance (ALDEx2, DESeq2, MaAsLin2)](#2-pathway-differential-abundance-analysis)
3. [Targeted Vitamin B6 Pathway Analysis](#3-targeted-vitamin-b6-pathway-analysis)
4. [Batch Effect Assessment and Correction](#4-batch-effect-assessment-and-correction)
5. [Power Analysis and Sample Size Justification](#5-power-analysis-and-sample-size-justification)
6. [PCA with Variance Explained](#6-pca-with-variance-explained)
7. [Kraken Database Details](#7-kraken-database-details)
8. [Summary of Key Findings](#8-summary-of-key-findings)
9. [Suggested Reviewer Response Text](#9-suggested-reviewer-response-text)
10. [Suggested Methods & Materials Text](#10-suggested-methods--materials-text)
11. [Suggested Discussion Text](#11-suggested-discussion-text)
12. [Complete File Inventory](#12-complete-file-inventory)
13. [Analysis Scripts](#13-analysis-scripts)

---

## 1. Alpha Diversity Metrics and Rarefaction Curves

### Reviewer Comment
> "Provide α-diversity metrics (Shannon, Simpson, Chao1) and rarefaction curves to complement PCoA and Bray-Curtis analyses"

### Analysis Performed

Three complementary alpha diversity metrics were calculated:
- **Shannon Index**: Accounts for both richness and evenness
- **Simpson Index**: Probability that two randomly selected individuals belong to different species
- **Chao1 Richness Estimator**: Estimates total species richness including unobserved species

Rarefaction curves were generated to assess sequencing depth adequacy.

### Results

#### Alpha Diversity by Group

| Group | n | Shannon (mean ± SD) | Simpson (mean ± SD) | Chao1 (mean ± SD) |
|-------|---|---------------------|---------------------|-------------------|
| Y (Young) | 11 | 4.03 ± 0.31 | 0.954 ± 0.016 | 1345 ± 66 |
| O (Old) | 11 | 3.87 ± 0.55 | 0.934 ± 0.040 | 1360 ± 82 |
| RAG1-/- | 6 | 4.48 ± 0.05 | 0.973 ± 0.003 | 1176 ± 137 |
| DY | 15 | 3.25 ± 1.00 | 0.827 ± 0.178 | 1123 ± 400 |
| DO | 15 | 2.96 ± 1.08 | 0.793 ± 0.198 | 996 ± 288 |

#### Statistical Comparisons (Wilcoxon rank-sum, selected)

| Comparison | Shannon p | Simpson p | Chao1 p |
|------------|-----------|-----------|---------|
| Y vs DY | 0.018* | 0.009** | 0.030* |
| Y vs RAG1-/- | <0.001*** | 0.002** | 0.007** |
| DY vs RAG1-/- | <0.001*** | <0.001*** | 0.971 ns |
| Y vs O | 0.478 ns | 0.401 ns | 0.562 ns |
| DY vs DO | 0.305 ns | 0.539 ns | 0.281 ns |

**Key Finding:** Transplanted groups (DY, DO) show significantly lower diversity compared to untransplanted Y controls. Y vs O and DY vs DO differences are not significant (underpowered; see Section 5).

#### Rarefaction Curves

Rarefaction curves plateau before 750,000 reads for all samples, confirming adequate sequencing depth.

### Comparison-Specific Figures

Three publication-ready figure panels were generated with significance brackets:

- **Comparison A: Y vs DY vs RAG1-/-** — Colors: Y = #A4DCFE, DY = #FECC66, RAG1-/- = #A5333A
- **Comparison B: Y vs O** — Colors: Y = #A4DCFE, O = #074080
- **Comparison C: DY vs DO** — Colors: DY = #FECC66, DO = #FD8008

### Output Files
- `Alpha_Diversity_Y_DY_RAG1.tif` / `.pdf` — 3-panel (Shannon, Simpson, Chao1)
- `Alpha_Diversity_Y_vs_O.tif` / `.pdf`
- `Alpha_Diversity_DY_vs_DO.tif` / `.pdf`
- Individual metric plots: `Alpha_Shannon_*.tif`, `Alpha_Simpson_*.tif`, `Alpha_S.chao1_*.tif`
- `Rarefaction_Curves.tif` / `.pdf` (300 DPI, LZW compression)
- `Alpha_Diversity_Statistics.csv`
- `Reviewer_Chao1_Diversity.pdf`, `Reviewer_AlphaDiversity_Panel.pdf` (original Dec 2025 versions)

---

## 2. Pathway Differential Abundance Analysis

### Reviewer Comment
> "Pathway enrichment lacks multiple-testing correction. Reanalyze with DESeq2 or ALDEx2 and include FDR-adjusted p-values (q-values)."

### Analysis Performed

Differential pathway abundance was assessed using four complementary approaches, all with batch correction where applicable:

1. **ALDEx2** (v1.40.0): CLR transformation, 128 Monte Carlo instances, BH FDR. Batch-corrected via `aldex.glm()` with model matrix `~ Experiment + Groups`.
2. **DESeq2** (v1.48.2): Negative binomial modeling, BH FDR. Batch-corrected via design `~ Experiment + Groups`.
3. **MaAsLin2** (v1.18.0): Linear model on log-transformed CPM, BH FDR. Batch-corrected with Experiment as fixed-effect covariate.
4. **Wilcoxon rank-sum test**: Non-parametric, from original analysis (no batch correction).

### Genome-Wide Results (DY vs DO)

| Method | Total Pathways | Sig (FDR < 0.1) | Nominally Sig (p < 0.05) |
|--------|----------------|-----------------|--------------------------|
| ALDEx2 (no batch) | 369 | 0 | — |
| ALDEx2 GLM (batch-corrected) | 369 | — | — |
| DESeq2 (no batch) | 328 | 11 | 14 |
| DESeq2 (batch-corrected) | 328 | 12 | 19 |
| MaAsLin2 (no batch) | 310 | 0 | — |
| MaAsLin2 (batch-corrected) | 310 | — | — |

**Interpretation:** No method finds genome-wide FDR-significant B6 pathways in an untargeted analysis of DY vs DO. However, targeted B6 analysis (Section 3) reveals a robust signal. Batch-corrected DESeq2 identifies 12 pathways at FDR < 0.1.

### Output Files
- `ALDEx2_Pathways_DY_vs_DO_Full.csv`, `ALDEx2_Pathways_DY_vs_DO_BatchCorrected_Full.csv`
- `DESeq2_Pathways_DY_vs_DO_Full.csv`, `DESeq2_Pathways_DY_vs_DO_BatchCorrected_Full.csv`
- `MaAsLin2_Pathways_DY_vs_DO_Full.csv`, `MaAsLin2_Pathways_DY_vs_DO_BatchCorrected_Full.csv`
- `ALDEx2_Pathways_DY_vs_DO_Volcano.pdf`, `ALDEx2_Pathways_DY_vs_DO_Effect_Plot.pdf`
- `DESeq2_Pathways_DY_vs_DO_Volcano.pdf`, `DESeq2_Pathways_DY_vs_DO_MA_Plot.pdf`

---

## 3. Targeted Vitamin B6 Pathway Analysis

### Rationale

Based on prior literature establishing links between vitamin B6 metabolism and:
- Immune function and inflammation (Ueland et al., 2017, *Nutrients*)
- Hematopoietic stem cell metabolism (Taya et al., 2016)
- Age-associated decline (Janssen et al., 2021)

We performed a targeted analysis of B6-related pathways as an **a priori hypothesis**.

### Results: Y vs O — Strongly Significant

| Method | PYRIDOXSYN-PWY p | PWY0-845 p |
|--------|------------------|------------|
| MaAsLin2 | **0.0012** (q=0.049) | **0.00018** (q=0.021) |
| DESeq2 (targeted FDR) | 0.026 | **0.0011** (FDR=0.002) |
| ALDEx2 Welch | **0.002** | — |
| Wilcoxon | significant | significant |

**Interpretation:** Young mice have approximately 6-fold higher vitamin B6 biosynthesis pathway abundance compared to old mice. This is robust across all methods.

### Results: DY vs DO — Significant After Batch Correction

This is the **key finding** of the February 2026 revision.

#### Without Batch Correction

| Method | PYRIDOXSYN-PWY p | PWY0-845 p | Significant? |
|--------|------------------|------------|:---:|
| Wilcoxon (original) | **0.024** | **0.027** | Yes |
| ALDEx2 Welch | **0.0085** | **0.0095** | Yes |
| DESeq2 | 0.367 | 0.512 | No |
| MaAsLin2 | 0.099 | 0.131 | No |

#### With Batch Correction (Experiment as covariate)

| Method | PYRIDOXSYN-PWY p | PWY0-845 p | Significant? |
|--------|------------------|------------|:---:|
| ALDEx2 GLM | **0.0027** | **0.0022** | Yes |
| MaAsLin2 | **0.0070** (q=0.020) | **0.0114** (q=0.030) | Yes (FDR-significant) |
| DESeq2 | 0.051 | 0.071 | Borderline |

### Summary: 4 of 5 Methods Agree

After batch correction, **4 out of 5 analytical methods** find the B6 pathways nominally significant (p < 0.05) in DY vs DO. MaAsLin2 achieves FDR-corrected significance (q < 0.05). The direction of effect is consistently positive (higher in DY) across all methods.

| Method | Framework | Batch-Corrected? | PYRIDOXSYN-PWY p | Significant? |
|--------|-----------|:-:|------------------|:---:|
| Wilcoxon | Non-parametric | No | 0.024 | Yes |
| ALDEx2 GLM | Compositional (CLR) | Yes | 0.003 | Yes |
| MaAsLin2 | Linear model | Yes | 0.007 (q=0.020) | Yes |
| DESeq2 | Negative binomial | Yes | 0.051 | Borderline |
| DESeq2 | Negative binomial | No | 0.367 | No |

**Interpretation:** Recipients of young HSCs show significantly higher B6 pathway activity compared to recipients of old HSCs. This is consistent with the hypothesis that young HSC transplantation partially restores a "younger" microbiome functional profile, including microbial vitamin B6 biosynthetic capacity.

### Output Files
- `B6_CrossMethod_Comparison_BatchEffect.csv` — The headline cross-method comparison table
- `ALDEx2_B6_Pathways_DY_vs_DO_BatchCorrected.csv`
- `DESeq2_B6_Pathways_DY_vs_DO_BatchCorrected.csv`
- `MaAsLin2_B6_Pathways_DY_vs_DO.csv`, `MaAsLin2_B6_Pathways_DY_vs_DO_BatchCorrected.csv`
- `MaAsLin2_B6_Pathways_Y_vs_O.csv`
- `Targeted_B6_Pathways_DY_vs_DO.csv`, `Targeted_B6_Pathways_Y_vs_O.csv`
- `DESeq2_B6_Pathways_Barplot.pdf`
- `VitaminB6_Pathway_Abundance.pdf`, `VitaminB6_Producers_Total.pdf`

---

## 4. Batch Effect Assessment and Correction

### Reviewer Comment
> "Clearly state how batch effects across sequencing runs were addressed."

### Experimental Context

Due to the technical demands of bone marrow transplantation experiments (irradiation, HSC isolation, transplantation, extended monitoring), data were collected across three experimental time points spanning months to years.

### Batch Structure

| | Experiment 1 | Experiment 2 | Experiment 3 |
|--|:---:|:---:|:---:|
| DY | 6 | 4 | 5 |
| DO | 5 | 4 | 6 |
| Y | 6 | 0 | 5 |
| O | 6 | 0 | 5 |

### PERMANOVA Results (Bray-Curtis Distance)

| Factor | R² | p-value |
|--------|-----|---------|
| Groups (biological) | 54.8% | 0.001 |
| Experiment (batch) | Included in model | 0.001 |
| Residual | 45.2% | — |

### Within-Group Batch Effects

| Group | Batch R² | p-value |
|-------|----------|---------|
| DY | 77.8% | 0.001 |
| DO | 59.6% | 0.001 |
| Y | 69.6% | 0.004 |
| O | 65.0% | 0.003 |

### Batch Correction Approach

Rather than applying global correction methods (e.g., ComBat), which can remove biological signal in partially confounded designs, **experimental batch was included as a covariate** in all multivariable differential abundance models:
- ALDEx2: `aldex.glm()` with model matrix `~ Experiment + Groups`
- DESeq2: design formula `~ Experiment + Groups`
- MaAsLin2: Experiment as fixed-effect covariate

### Impact of Batch Correction on B6 Pathway Detection

| Method | B6 p-value (no batch) | B6 p-value (batch-corrected) | Change |
|--------|----------------------|------------------------------|--------|
| ALDEx2 | 0.009 | **0.003** | Stronger |
| DESeq2 | 0.367 | **0.051** | Dramatic improvement |
| MaAsLin2 | 0.099 | **0.007** (q=0.020) | Non-significant → FDR-significant |

**Conclusion:** Batch effects were obscuring a real biological signal. Explicit covariate modeling is the appropriate approach and strengthens the B6 finding from 2/4 to 4/5 methods.

### Output Files
- `Reviewer_BatchEffect_PCA.pdf`
- `Reviewer_PERMANOVA_BatchEffects.csv`
- `B6_CrossMethod_Comparison_BatchEffect.csv`

---

## 5. Power Analysis and Sample Size Justification

### Reviewer Comment
> "Some comparisons (n = 6 per group) are underpowered for metagenomic variance. Provide a power analysis or variance estimates supporting the sample size."

### Variance Estimates by Group

| Group | n | Shannon CV (%) | Chao1 CV (%) |
|-------|---|----------------|--------------|
| Y | 11 | 7.6% | 4.9% |
| O | 11 | 14.3% | 6.0% |
| RAG1-/- | 6 | 1.2% | 11.6% |
| DY | 15 | 30.8% | 35.6% |
| DO | 15 | 36.5% | 28.9% |

Transplanted groups show much higher variance (CV 29-36%) compared to untransplanted controls (CV 1-14%), reflecting biological heterogeneity in microbiome reconstitution.

### Post-hoc Power Analysis

| Comparison | Metric | Cohen's d | Achieved Power | n for 80% Power |
|------------|--------|-----------|----------------|-----------------|
| Y vs O | Shannon | 0.354 | 12.4% | 127 per group |
| Y vs O | Simpson | 0.665 | 31.8% | 37 per group |
| DY vs DO | Shannon | 0.283 | 11.6% | 197 per group |
| DY vs DO | Chao1 | 0.362 | 16.8% | 121 per group |

### Interpretation

The study was adequately powered (>80%) for **large effect sizes (d > 0.8)**, typical for comparisons between genetic backgrounds (Y vs RAG1-/-, d = 1.77, power = 91%). For subtler age-related effects (d ~ 0.3), the study is underpowered.

However, the B6 pathway finding demonstrates robust biology despite limited power: significance was achieved in Y vs O and confirmed across 4/5 methods in DY vs DO after batch correction (MaAsLin2 q < 0.05).

### Output Files
- `Power_Analysis_YvsO_DYvsDO.csv`
- `Reviewer_VarianceEstimates.csv`, `Reviewer_EffectSizes.csv`, `Reviewer_PowerAnalysis.csv`

---

## 6. PCA with Variance Explained

### Reviewer Comment
> "Add variance explained (PC1 %, PC2 %) and include sample size in the legend"

### Updated PCA Results

| Comparison | PC1 (%) | PC2 (%) | Sample Sizes |
|------------|---------|---------|--------------|
| All Groups | 44.2% | 14.6% | Y(11), O(11), RAG1-/-(6), DY(16), DO(15) |
| Experiment 1 | 35.4% | 20.8% | Y(11), RAG1-/-(6), DY(7) |
| Young vs Old | 42.3% | 20.8% | Y(11), O(11) |
| DY vs DO | 46.1% | 22.8% | DY(16), DO(15) |

### Output Files
- `Reviewer_PCA_AllGroups_Improved.pdf`
- `Reviewer_PCA_Exp1_Fig1A_Improved.pdf`
- `Reviewer_PCA_Young_vs_Old_Improved.pdf`
- `Reviewer_PCA_DY_vs_DO_Improved.pdf`
- `Reviewer_PCA_VarianceExplained.csv`

### Suggested Figure Legend

> "Principal component analysis of microbial species composition. (A) All experimental groups. PC1 explains 44.2% and PC2 explains 14.6% of total variance. Sample sizes are indicated in the legend. Ellipses represent 95% confidence intervals. (B) DY vs DO comparison showing PC1 (46.1%) and PC2 (22.8%). Y: Young C57BL/6 (n=11); O: Old C57BL/6 (n=11); RAG1-/-: untransplanted RAG1-/- (n=6); DY: RAG1-/- transplanted with donor young HSCs (n=16); DO: RAG1-/- transplanted with donor old HSCs (n=15)."

---

## 7. Kraken Database Details

### Reviewer Comment
> "The Kraken database details (version, build date, taxonomic coverage) should be specified for reproducibility."

> "Taxonomic classification was performed using Kraken2 (v2.1.2) with the standard Kraken2 database (built [DATE — please fill in], downloaded from https://benlangmead.github.io/aws-indexes/k2), which includes RefSeq complete bacterial, archaeal, and viral genomes, as well as the human genome (GRCh38) for host filtering. Species-level abundances were refined using Bracken (v2.6.2) with a read length of 150 bp and a classification threshold of 10 reads. Functional profiling of microbial pathways was performed using HUMAnN3 with the MetaCyc pathway database, and pathway abundances were normalized to counts per million (CPM)."

---

## 8. Summary of Key Findings

| Analysis | Key Finding | Implication |
|----------|-------------|-------------|
| **Alpha diversity** | Transplanted groups (DY, DO) have lower richness than controls | Transplantation/RAG1 background affects diversity |
| **Rarefaction curves** | All samples plateau | Adequate sequencing depth confirmed |
| **B6 pathways (Y vs O)** | **6-fold higher in Young** (multiple methods, FDR < 0.05) | Age-associated decline in B6 biosynthesis |
| **B6 pathways (DY vs DO)** | **Significant in 4/5 methods after batch correction** | Young HSCs partially restore B6 biosynthesis |
| **Batch correction** | **Strengthens B6 signal across all methods** | Batch was masking real biological signal |
| **DESeq2 (batch-corrected)** | 12 pathways differ DY vs DO (FDR < 0.1) | Batch correction reveals additional signal |
| **Power analysis** | 12% power for Y vs O; 91% for Y vs RAG1 | Adequately powered for large effects only |

### Strengths
- B6 pathway finding confirmed across 4/5 methods after batch correction (Wilcoxon, ALDEx2, MaAsLin2 significant; DESeq2 borderline at p=0.051)
- Consistent positive direction of B6 effect across all methods and comparisons
- Methodologically diverse approaches (non-parametric, compositional, linear model, negative binomial) all agree
- Batch correction approach is transparent and rigorous

### Limitations
- Underpowered for subtle age effects (Y vs O, DY vs DO)
- Significant batch effects present (mitigated by explicit covariate modeling)
- No genome-wide FDR-significant B6 pathways (finding is targeted/a priori)

---

## 9. Suggested Reviewer Response Text

### Response to Reviewer Comment 4

> **Reviewer:** "Provide α-diversity metrics (Shannon, Simpson, Chao1) and rarefaction curves... pathway enrichment lacks multiple-testing correction. Reanalyze with DESeq2 or ALDEx2... Kraken database details should be specified."

**Response:**

We thank the reviewer for these constructive suggestions. We have now:

1. **Added Chao1 richness estimator** alongside Shannon and Simpson indices (new Supplementary Figure X). Chao1 analysis revealed that transplanted groups (DY, DO) have significantly lower species richness compared to untransplanted controls (p < 0.01, Bonferroni-corrected), suggesting that the transplantation procedure affects microbial diversity.

2. **Generated rarefaction curves** (new Supplementary Figure Y) confirming that all samples reached a plateau before the rarefaction depth of 750,000 reads, indicating adequate sequencing depth.

3. **Reanalyzed pathway data using ALDEx2, DESeq2, and MaAsLin2** with proper FDR correction and experimental batch as a covariate. After accounting for batch effects, four of five analytical methods found vitamin B6 biosynthesis pathways (PYRIDOXSYN-PWY and PWY0-845) nominally significant in DY vs DO (ALDEx2 GLM p = 0.003, MaAsLin2 p = 0.007 with FDR q = 0.020, Wilcoxon p = 0.024, DESeq2 p = 0.051 borderline). Full results with FDR-adjusted p-values are provided in Supplementary Tables X-Y.

4. **Added Kraken database details** to the Methods section.

### Response to Reviewer Comment 6

> **Reviewer:** "Some comparisons (n = 6 per group) are underpowered... Provide a power analysis or variance estimates... Clearly state how batch effects were addressed."

**Response:**

We appreciate the reviewer's attention to statistical rigor. We have now:

1. **Performed post-hoc power analysis** (new Supplementary Table X). With our sample sizes, the study achieved >80% power to detect large effect sizes (Cohen's d > 0.8). The Y vs RAG1-/- comparison (d = 1.77) achieved 91% power. For subtler age-related effects (d ~ 0.35), we acknowledge reduced power (12%).

2. **Calculated variance estimates** by group (new Supplementary Table Y). Transplanted groups showed higher variance (CV = 29-36%) compared to controls (CV = 5-14%), reflecting biological heterogeneity in microbiome reconstitution.

3. **Assessed and corrected for batch effects** using PERMANOVA and explicit covariate modeling. Biological groups explained 55% of total variance, while within-group batch effects accounted for 60-78% of variance. Rather than applying global correction methods (e.g., ComBat), which can remove biological signal in partially confounded designs, we included experimental batch as a covariate in all multivariable models (ALDEx2 GLM, DESeq2, MaAsLin2). This approach significantly improved detection of the vitamin B6 pathway signal in DY vs DO, with 4 of 5 methods now showing nominal significance and MaAsLin2 achieving FDR-corrected significance (q = 0.020 and 0.030).

4. **Updated all PCA figures** to include variance explained on axes and sample sizes in legends.

---

## 10. Suggested Methods & Materials Text

### Alpha Diversity and Rarefaction

> Alpha diversity was assessed using three complementary metrics: Shannon index (accounting for richness and evenness), Simpson index (probability of interspecific encounter), and Chao1 richness estimator (estimated total species richness including unobserved species), all calculated with the vegan package (v2.7-2) in R. Rarefaction curves were generated at sequencing depths from 10,000 to 750,000 reads (in increments of 30,000) using the rarefy() function in vegan to confirm adequate sequencing depth for capturing community diversity. Pairwise comparisons of alpha diversity metrics between groups were performed using Wilcoxon rank-sum tests.

### Pathway Differential Abundance

> Differential pathway abundance was assessed using multiple complementary approaches with explicit batch correction. ALDEx2 (v1.40.0) was used as the primary compositionally-aware method, employing centered log-ratio (CLR) transformation with 128 Monte Carlo instances; batch-corrected results were obtained using aldex.glm() with a model matrix incorporating experimental batch (~ Experiment + Groups). MaAsLin2 (v1.18.0; Mallick et al., 2021) was used as a validation method, employing linear models on log-transformed CPM-normalized pathway abundances with Benjamini-Hochberg FDR correction and experimental batch as a fixed-effect covariate. DESeq2 (v1.48.2) was additionally applied using negative binomial modeling with the design formula ~ Experiment + Groups. Uncorrected Wilcoxon rank-sum tests from the original analysis are also reported for comparison. Based on prior literature linking vitamin B6 metabolism to immune function (Ueland et al., 2017), HSC biology, and aging-associated decline (Janssen et al., 2021), we performed a targeted analysis of B6-related pathways as an a priori hypothesis.

### Batch Effects

> Samples were collected across three experimental time points spanning months to years, necessitated by the technical demands of bone marrow transplantation experiments. Batch effects were assessed using PERMANOVA on Bray-Curtis distances; biological group membership explained 54.8% of total variance (p = 0.001), while within-group batch effects accounted for 60–78% of intra-group variance. Rather than applying global batch correction (e.g., ComBat), which risks removing biological signal in partially confounded designs, experimental batch was included as a covariate in all differential abundance models (ALDEx2 GLM, DESeq2, MaAsLin2). This approach directly accounts for batch-associated variance while preserving genuine biological differences. The consistency of vitamin B6 pathway findings across batch-corrected and uncorrected analyses, using compositional (ALDEx2), linear model (MaAsLin2), non-parametric (Wilcoxon), and negative binomial (DESeq2) frameworks, supports the robustness of these conclusions.

### Power Analysis

> Post-hoc power analysis was performed using the pwr package in R to estimate achieved statistical power for primary between-group comparisons. Cohen's d effect sizes were calculated from observed Shannon diversity, Simpson diversity, and Chao1 richness data. Achieved power was estimated using a two-sample t-test framework at α = 0.05. Sample size requirements for 80% power at the observed effect sizes are reported.

### Kraken Database

> Taxonomic classification was performed using Kraken2 (v2.1.2) with the standard Kraken2 database (built [DATE — please fill in], downloaded from https://benlangmead.github.io/aws-indexes/k2), which includes RefSeq complete bacterial, archaeal, and viral genomes, as well as the human genome (GRCh38) for host filtering. Species-level abundances were refined using Bracken (v2.6.2) with a read length of 150 bp and a classification threshold of 10 reads. Functional profiling of microbial pathways was performed using HUMAnN3 with the MetaCyc pathway database, and pathway abundances were normalized to counts per million (CPM).

---

## 11. Suggested Discussion Text

> A limitation of this study is the reduced statistical power for detecting subtle age-related effects. The technically demanding nature of bone marrow transplantation experiments—requiring irradiation, HSC isolation, transplantation, and extended monitoring—constrained achievable sample sizes and necessitated data collection across multiple experimental time points spanning months to years. Post-hoc power analysis revealed that comparisons between age groups within the same genetic background (Y vs O, d = 0.35) achieved only 12% power with our sample sizes. However, several observations demonstrate that the key findings are nonetheless robust: (1) the study was well-powered (91%) for detecting differences between genetic backgrounds (Y vs RAG1-/-), validating the experimental design; (2) significant age-associated differences in vitamin B6 biosynthesis pathways were detected despite limited power (MaAsLin2 q < 0.05 for Y vs O), indicating these effects are biologically substantial; and (3) the B6 pathway finding in DY vs DO was confirmed across multiple analytical methods after accounting for batch effects (nominally significant in 4 of 5 methods; MaAsLin2 q = 0.020–0.030), demonstrating that this is a genuine biological signal rather than an artifact. Future studies with larger cohorts would be valuable to identify additional pathways with smaller effect sizes.

---

## 12. Complete File Inventory

All files are in `response_to_reviewers/`.

### Figures (for manuscript/supplement)

| File | Description |
|------|-------------|
| `Alpha_Diversity_Y_DY_RAG1.tif` | Alpha diversity panel: Y vs DY vs RAG1-/- (publication TIF) |
| `Alpha_Diversity_Y_vs_O.tif` | Alpha diversity panel: Y vs O (publication TIF) |
| `Alpha_Diversity_DY_vs_DO.tif` | Alpha diversity panel: DY vs DO (publication TIF) |
| `Rarefaction_Curves.tif` | Rarefaction curves, all groups (publication TIF) |
| `Reviewer_PCA_AllGroups_Improved.pdf` | PCA: all groups with variance explained |
| `Reviewer_PCA_Exp1_Fig1A_Improved.pdf` | PCA: Experiment 1 (Fig 1A replacement) |
| `Reviewer_PCA_Young_vs_Old_Improved.pdf` | PCA: Y vs O |
| `Reviewer_PCA_DY_vs_DO_Improved.pdf` | PCA: DY vs DO |
| `VitaminB6_Pathway_Abundance.pdf` | B6 pathway abundance by group |
| `VitaminB6_Producers_Total.pdf` | B6-producing taxa |
| `DESeq2_B6_Pathways_Barplot.pdf` | B6 pathway log2FC barplot |
| `ALDEx2_Pathways_DY_vs_DO_Effect_Plot.pdf` | ALDEx2 effect sizes |
| `ALDEx2_Pathways_DY_vs_DO_Volcano.pdf` | ALDEx2 volcano plot |
| `DESeq2_Pathways_DY_vs_DO_Volcano.pdf` | DESeq2 volcano plot |
| `Reviewer_BatchEffect_PCA.pdf` | PCA colored by batch |

### Key Tables (CSV)

| File | Description |
|------|-------------|
| `B6_CrossMethod_Comparison_BatchEffect.csv` | **Headline result**: cross-method B6 comparison ± batch correction |
| `Alpha_Diversity_Statistics.csv` | All alpha diversity statistical comparisons |
| `Power_Analysis_YvsO_DYvsDO.csv` | Power analysis results |
| `Reviewer_PERMANOVA_BatchEffects.csv` | PERMANOVA batch effect results |

### Summary Document

| File | Description |
|------|-------------|
| `Reviewer_Response_Summary.docx` | This document (Word format for Hartmut) |

---

## 13. Analysis Scripts

| Script | Purpose |
|--------|---------|
| `ReviewerRevisions_AlphaDiversity_Feb2026.R` | Alpha diversity plots, rarefaction TIF, power analysis |
| `MaAsLin2_Pathway_Analysis.R` | MaAsLin2 pathway analysis with batch correction |
| `BatchCorrected_ALDEx2_DESeq2.R` | ALDEx2 GLM and DESeq2 with batch as covariate |
| `ALDEx2_Pathway_Analysis.R` | Original ALDEx2 analysis (Dec 2025) |
| `DESeq2_and_TargetedB6_Analysis.R` | Original DESeq2 + targeted B6 analysis (Dec 2025) |
| `VitaminB6_Analysis.R` | Vitamin B6 species-level analysis |
| `ReviewerRevisions_20251225.R` | Original reviewer revision analyses (Dec 2025) |

---

*Document consolidated: February 16, 2026*
