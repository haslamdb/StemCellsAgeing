# Stem Cell Aging Microbiome Analysis

Metagenomic analysis of gut microbiome changes in a hematopoietic stem cell (HSC) transplantation aging model.

## Study Design

Five experimental groups:
- **Y**: Young C57BL/6 mice (untransplanted controls)
- **O**: Old C57BL/6 mice (untransplanted controls)
- **RAG1-/-**: Untransplanted RAG1-deficient mice
- **DY**: RAG1-/- mice transplanted with young donor HSCs
- **DO**: RAG1-/- mice transplanted with old donor HSCs

## Analysis Scripts

| Script | Description |
|--------|-------------|
| `GeigerStemCellAgeing_tidy_20250512.R` | Main analysis: species diversity, PCA, differential abundance (SDA) |
| `ReviewerRevisions_20251225.R` | Reviewer-requested: Chao1, rarefaction, power analysis, batch effects |
| `ALDEx2_Pathway_Analysis.R` | ALDEx2 differential pathway abundance with FDR correction |
| `DESeq2_and_TargetedB6_Analysis.R` | DESeq2 pathway analysis + targeted vitamin B6 hypothesis testing |
| `VitaminB6_Analysis.R` | Targeted analysis of B6 biosynthesis pathways and producing species |

## Key Outputs

### Reviewer Response Materials
- `Reviewer_Response_Summary.md` - Comprehensive summary with suggested manuscript/reviewer text

### Alpha Diversity
- `Reviewer_Chao1_Diversity.pdf`
- `Reviewer_AlphaDiversity_Panel.pdf`
- `Reviewer_RarefactionCurves.pdf`

### Pathway Analysis (ALDEx2 & DESeq2)
- `ALDEx2_Pathways_*_Full.csv` - Complete ALDEx2 results with FDR
- `DESeq2_Pathways_*_Full.csv` - Complete DESeq2 results with FDR
- `*_Volcano.pdf` - Volcano plots

### Vitamin B6 Analysis
- `Targeted_B6_Pathways_*.csv` - B6 pathway results with targeted FDR
- `VitaminB6_Pathway_Abundance.pdf`
- `VitaminB6_Producers_Total.pdf`

### Power & Batch Analysis
- `Reviewer_PowerAnalysis.csv`
- `Reviewer_VarianceEstimates.csv`
- `Reviewer_PERMANOVA_BatchEffects.csv`
- `Reviewer_BatchEffect_PCA.pdf`

## Data Files (Not Tracked)

Large data files are excluded from this repository:
- `GeigerData*` - Processed Kraken2/Bracken/HUMAnN3 outputs
- `DataBefore*/` - Archived intermediate data

## Methods Summary

- **Taxonomic classification**: Kraken2 + Bracken
- **Functional profiling**: HUMAnN3 (MetaCyc pathways)
- **Differential abundance**: SDA, ALDEx2 (CLR + Monte Carlo), DESeq2 (negative binomial)
- **Multiple testing correction**: Benjamini-Hochberg FDR
- **Diversity metrics**: Shannon, Simpson, Chao1 (vegan)

## Key Findings

1. **B6 biosynthesis declines with age**: Young mice show 6-fold higher vitamin B6 pathway abundance vs old (FDR = 0.002)
2. **Transplantation reduces diversity**: DY and DO groups have significantly lower Chao1 richness than untransplanted controls
3. **Genetic background dominates**: RAG1-/- vs C57BL/6 shows very large effect sizes (Cohen's d = 1.77)
