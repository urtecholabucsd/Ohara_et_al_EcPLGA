# MO_EcoliVax: supplementary microbiome figure

Reproduces the three-panel supplementary microbiomics figure for the E. coli OMV-PLGA vaccine study.

- **A.** Family-level relative abundance (group means, top 12 families + Other), T0 to T5
- **B.** Within-subject Bray-Curtis distance from each mouse's own T0 sample (mean +/- SEM)
- **C.** Area under each mouse's BC-from-T0 curve, by group (Mann-Whitney vs RBC-PLGA)

Panels are combined with `cowplot::plot_grid`.

## Layout

| Folder | Contents |
|---|---|
| `input/` | `20260527_MO_EcoliVax.csv`, preprocessed long-format zOTU table (seqid, OTU, reads, relabd, Mouse, Day, Vaccine Administered) |
| `ref/` | `tax.csv`, zOTU taxonomy |
| `scripts/` | `Supplementary_Microbiomics.R` |
| `figs/` | `Supplementary_microbiomics.pdf` and `.png` |
| `output/` | Statistics (Kruskal-Wallis, Mann-Whitney, PERMANOVA) and per-mouse BC/AUC tables |

## Running

From the repo root:

```bash
Rscript scripts/Supplementary_Microbiomics.R
```

or open the script in RStudio and click **Source**.

Required R packages: `tidyverse`, `cowplot`, `vegan`, `ggsignif`.

## Note on interpretation

Cage is fully confounded with treatment (2 cages per group), and groups already differ by PERMANOVA at Day 0, before any injection. The within-subject measures in panels B and C use each mouse as its own control to sidestep that baseline difference.
