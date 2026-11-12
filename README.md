# MO_EcoliVax: Q3 baseline-corrected microbiome analysis

Within-subject analyses that use each mouse's Day 0 sample as its own control, to separate vaccine effects from cage/baseline confounding.

1. Bray-Curtis distance from Day 0 (line plot, per-day box plots, LMM, Mann-Whitney U)
2. PCoA on CLR-delta from Day 0 (Aitchison distance, per-day panels in one plot grid, pairwise PERMANOVA)

## Layout

| Folder | Contents |
|---|---|
| `input/` | `zotu.csv` (zOTU count table), `20260527_MO_EcoliVax.csv` (processed sample table from 16S preprocessing) |
| `ref/` | `tax.csv` (zOTU taxonomy), `MO_Metadata.xlsx` (sample/mouse/cage metadata) |
| `scripts/` | `Q3_BaselineCorrected.Rmd` |
| `figs/` | Figures written by the script |
| `output/` | Statistics tables written by the script |

## Running

Open `scripts/Q3_BaselineCorrected.Rmd` in RStudio and click **Knit**, or from the repo root:

```r
rmarkdown::render("scripts/Q3_BaselineCorrected.Rmd")
```

Required R packages: `tidyverse`, `cowplot`, `lme4`, `lmerTest`, `vegan`, `ggsignif`, `readxl`, `rmarkdown`.

## Outputs

- `figs/Q3_bc_from_baseline_line.png`
- `figs/Q3_bc_from_baseline_boxplots.png`
- `figs/Q3_clr_delta_pcoa_per_day.png` (per-day PCoA panels arranged with `cowplot::plot_grid`, shared legend)
- `output/Q3_bc_from_baseline_mannwhitney.csv`
- `output/Q3_clr_delta_permanova.csv`
