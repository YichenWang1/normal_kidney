# Systemic mutagen exposures reported by normal kidney cell genomes

This repository contains the customized analysis code for the manuscript Wang  *et al.* *Systemic mutagen exposures reported by normal kidney cell genomes*. 

---

## Repository Structure

```
.
├── data/                        # Input data files
├── Fig2/                        # Figure 2 – Signature overview
├── Fig3/                        # Figure 3 – Mutation rate analyses
├── Fig4/                        # Figure 4 – Example patient case studies
├── Supplementary/               # Supplementary figures
└── signature_variation/         # Hypothesis testing & enrichment analyses in the main text
```

---

## Scripts and Figures

### Figure 2 – Mutational Signature Overview
| Script | Description |
|--------|-------------|
| `Fig2/Fig2ab_signature_overview.Rmd` | SBS and indel signature compositions in normal and cancer bulk kidney, stratified by country |

### Figure 3 – Mutation Rate Analyses
| Script | Description |
|--------|-------------|
| `Fig3/Fig3ab_mutation_rate.Rmd` | Age-dependent SBS and indel mutation rates across nephron structures from donors without kidney cancer |
| `Fig3/Fig3c_mutation_rate_signature.R` | Per-signature age-dependent mutation rates across LCM structures, aggregating all donors from all country |
| `Fig3/Fig3d_LCM_SBS_overview.Rmd` | SBS signature composition across nephron structures |

### Figure 4 – Example Patient Case Studies
| Script | Description |
|--------|-------------|
| `Fig4/Fig4_example_signature_overview.Rmd` | Signature profiles across tissue types for four representative patients |
| `Fig4/plot_snv_spectrums.R` | Plotting SBS96 mutation spectra for individual samples |

### Supplementary Figures
| Script | Description |
|--------|-------------|
| `Supplementary/S1.Rmd` | Bulk normal kidney SBS and indel mutation burden by country |
| `Supplementary/S2.Rmd` | HDP *de novo* signature decomposition and cosine similarity to reference signatures |
| `Supplementary/S3.Rmd` | SigProfiler *de novo* signature decomposition and cosine similarity to reference signatures|
| `Supplementary/S4.Rmd` | Indel signature composition across nephron structures |
| `Supplementary/S6.Rmd` | Geographic variation in SBS12, SBS22a/b/c and ID23 burden |
| `Supplementary/S7.Rmd` | Country-level mutation rate estimates for SBS40b in bulk kidney, proximal tubules, and kidney cancer |

### Variation of Signatures across Different Nephron Structures and Epidemiology Analysis
| Script | Description |
|--------|-------------|
| `signature_variation/siganture_hypotheis_tesing_bulk.Rmd` | Hypothesis testing for age, sex, and geographic effects on signature burden in bulk samples |
| `signature_variation/signatures_enrichment_lcm.Rmd` | Signature enrichment in LCM biopsies across nephron structures and countries |

<!-- --- -->

<!-- ## Data -->
<!-- 
Input files are in `data/`. Key files:

| File | Description |
|------|-------------|
| `bulk_normal_kidney_metadata.csv` | Sample-level metadata for bulk normal kidney (country, age, sex, alcohol, tobacco, mutation burden) |
| `LCM_normal_kidney_metadata.csv` | Sample-level metadata for LCM samples (structure, patient, demographics, mutation burden) |
| `SBS96_kidney_signature_weights_table.csv` | Pre-computed SBS signature weights per sample |
| `ID83_kidney_signature_ID_weights_table.csv` | Pre-computed indel signature weights per sample |
| `bulk_SBS96_matrix_corrected.txt` | SBS96 mutation count matrix for bulk samples |
| `LCM_SBS96_matrix_corrected_merged.txt` | SBS96 mutation count matrix for LCM samples |
| `RCC_SBS96_matrix.txt` | SBS96 mutation count matrix for renal cell carcinoma |
| `COSMIC_v3.4_SBS_GRCh37_with_liver_platinum_sig.txt` | COSMIC v3.4 reference SBS signatures | -->

---

## Dependencies

All analyses and visualization are implemented in **R**. Key packages:
 `nlme`, `tidyverse`, `dplyr`, `tidyr`, `ggplot2`, `ggpubr`, `cowplot`, `RColorBrewer`, `lsa`, `emmeans`.

<!-- ---

## Citation

> Wang Y, *et al.* (2025). *Title TBC.* Journal TBC. -->

---

## Contact

Yichen Wang — [yw2@sanger.ac.uk](mailto:yw2@sanger.ac.uk)
Wellcome Sanger Institute

---

## License

MIT License — see [LICENSE](LICENSE) for details.
