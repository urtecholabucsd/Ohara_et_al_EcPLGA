# ============================================================
# Supplementary microbiomics figure
#   Panel A: Family-level relative abundance (stacked bar), group mean, T0-T5
#   Panel B: Within-subject Bray-Curtis distance from baseline (T0), T0-T5
#   Panel C: AUC of each mouse's BC-from-T0 trajectory, by group
# Combined with cowplot::plot_grid and saved to figs/.
# Statistics are written to output/.
#
# Input data: input/20260527_MO_EcoliVax.csv  (preprocessed, long format:
#   seqid, OTU, reads, relabd, log10relabd, Mouse, Day, Vaccine Administered)
# Reference:  ref/tax.csv  (OTU -> family/genus taxonomy)
#
# Run from the repo root (Rscript scripts/Supplementary_Microbiomics.R)
# or source it from RStudio; paths resolve relative to scripts/.
# ============================================================

library(tidyverse)
library(cowplot)
library(vegan)
library(ggsignif)

# --- Resolve working directory to scripts/ ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
    nzchar(rstudioapi::getActiveDocumentContext()$path)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else if (basename(getwd()) != "scripts" && dir.exists("scripts")) {
  setwd("scripts")
}

fig_dir <- "../figs"
out_dir <- "../output"
dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)

# --- Load preprocessed data ---
dat <- read_csv("../input/20260527_MO_EcoliVax.csv")
tax <- read_csv("../ref/tax.csv") %>% select(OTU, family, genus)

# Display labels use a true Delta for the lpxM mutant
lpxm_label <- "ΔlpxM EC-PLGA"
vax_order  <- c("RBC-PLGA", "WT EC-PLGA", lpxm_label)
vax_colors <- setNames(c("#636363", "#e6550d", "#3182bd"), vax_order)

# Day -> timepoint label (T0..T5)
week_map <- c("0" = "T0", "7" = "T1", "14" = "T2", "21" = "T3", "28" = "T4", "35" = "T5")

dat <- dat %>%
  left_join(tax, by = "OTU") %>%
  mutate(
    Day    = as.integer(Day),
    Mouse  = as.factor(Mouse),
    Week   = factor(week_map[as.character(Day)], levels = c("T0","T1","T2","T3","T4","T5")),
    family = coalesce(family, "Unclassified"),
    genus  = coalesce(genus,  "Unclassified"),
    `Vaccine Administered` = recode(`Vaccine Administered`, "-lpxM EC-PLGA" = lpxm_label),
    `Vaccine Administered` = factor(`Vaccine Administered`, levels = vax_order)
  )

# Inside-panel legend placement that works on ggplot2 < 3.5 and >= 3.5
legend_inside <- function(x, y) {
  if (packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside", legend.position.inside = c(x, y))
  } else {
    theme(legend.position = c(x, y))
  }
}

# ============================================================
# PANEL A: Family-level relative abundance stacked bar (top 12 + Other)
# ============================================================
rel_fam <- dat %>%
  group_by(seqid, Mouse, Day, Week, `Vaccine Administered`, family) %>%
  summarise(relabd = sum(relabd), .groups = "drop")

top_fam <- rel_fam %>%
  group_by(family) %>%
  summarise(mean_abd = mean(relabd), .groups = "drop") %>%
  arrange(desc(mean_abd)) %>%
  head(12) %>%
  pull(family)

fam_avg <- rel_fam %>%
  mutate(taxon = if_else(family %in% top_fam, family, "Other")) %>%
  group_by(`Vaccine Administered`, Week, taxon) %>%
  summarise(relabd = sum(relabd), .groups = "drop") %>%
  group_by(`Vaccine Administered`, Week) %>%
  mutate(relabd = relabd / sum(relabd)) %>%   # renormalize after collapsing
  ungroup() %>%
  mutate(taxon = factor(taxon, levels = c(top_fam, "Other")))

# Italicize family names in the legend (dropping SILVA-style "_2" suffixes),
# keep "Other" upright
fam_labels <- lapply(levels(fam_avg$taxon), function(l) {
  if (l == "Other") l else bquote(italic(.(sub("_[0-9]+$", "", l))))
})

pA <- ggplot(fam_avg, aes(x = Week, y = relabd, fill = taxon)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2, width = 0.95) +
  facet_wrap(~`Vaccine Administered`, nrow = 1) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_fill_discrete(labels = fam_labels) +
  labs(x = NULL, y = "Relative abundance", fill = "Family") +
  theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 14),
    legend.title     = element_text(face = "bold", size = 16),
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(0.4, "cm"),
    plot.margin      = margin(5.5, 25, 5.5, 5.5)
  )

# ============================================================
# PANEL B: Within-subject Bray-Curtis distance from baseline (T0)
#   Genus-level relative abundance from relabd
# ============================================================
genus_rel <- dat %>%
  group_by(seqid, genus) %>%
  summarise(relabd = sum(relabd), .groups = "drop") %>%
  pivot_wider(names_from = genus, values_from = relabd, values_fill = 0) %>%
  column_to_rownames("seqid")

sample_ids <- dat %>%
  select(seqid, Mouse, Day, Week, `Vaccine Administered`) %>% distinct()

d0_map <- sample_ids %>% filter(Day == 0) %>% select(Mouse, seqid) %>% deframe()

bc_from_baseline <- sample_ids %>%
  filter(Mouse %in% names(d0_map)) %>%
  rowwise() %>%
  mutate(
    bc_from_d0 = if (Day == 0) 0 else
      vegan::vegdist(
        rbind(genus_rel[seqid, ], genus_rel[d0_map[as.character(Mouse)], ]),
        method = "bray"
      )[1]
  ) %>%
  ungroup()

summary_bc <- bc_from_baseline %>%
  group_by(`Vaccine Administered`, Week) %>%
  summarise(mean_bc = mean(bc_from_d0),
            sem_bc  = sd(bc_from_d0) / sqrt(n()),
            .groups = "drop")

pB <- ggplot(summary_bc, aes(x = Week, y = mean_bc, color = `Vaccine Administered`,
                             group = `Vaccine Administered`)) +
  geom_errorbar(aes(ymin = mean_bc - sem_bc, ymax = mean_bc + sem_bc), width = 0.15) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = vax_colors) +
  labs(x = NULL, y = "Bray-Curtis distance from T0", color = "Vaccine Group") +
  theme_cowplot() +
  legend_inside(0.55, 0.22) +
  theme(legend.title = element_text(face = "bold", size = 16),
        legend.text  = element_text(size = 13),
        plot.margin  = margin(10, 5.5, 5.5, 18))

# ============================================================
# PANEL C: AUC of BC-from-baseline trajectory (per mouse), by group
#   Trapezoidal integration over Day (0..35); small minimal barplot + jitter
# ============================================================
auc_trapz <- function(x, y) {
  o <- order(x); x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

auc_df <- bc_from_baseline %>%
  group_by(Mouse, `Vaccine Administered`) %>%
  summarise(auc = auc_trapz(Day, bc_from_d0), .groups = "drop")

auc_summ <- auc_df %>%
  group_by(`Vaccine Administered`) %>%
  summarise(mean_auc = mean(auc), sem_auc = sd(auc) / sqrt(n()), .groups = "drop")

pC <- ggplot(auc_df, aes(x = `Vaccine Administered`, y = auc)) +
  geom_errorbar(data = auc_summ,
                aes(x = `Vaccine Administered`, y = mean_auc,
                    ymin = mean_auc - sem_auc, ymax = mean_auc + sem_auc),
                width = 0.25, inherit.aes = FALSE) +
  geom_col(data = auc_summ, aes(x = `Vaccine Administered`, y = mean_auc,
                                fill = `Vaccine Administered`),
           width = 0.7, alpha = 0.7, inherit.aes = FALSE) +
  geom_jitter(aes(color = `Vaccine Administered`), height = 0, width = 0.12,
              size = 1.3, alpha = 0.8) +
  geom_signif(comparisons = list(c("RBC-PLGA", "WT EC-PLGA"),
                                 c("RBC-PLGA", lpxm_label)),
              test = "wilcox.test", map_signif_level = FALSE,
              textsize = 3, step_increase = 0.22) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.18))) +
  scale_fill_manual(values = vax_colors) +
  scale_color_manual(values = vax_colors) +
  labs(x = NULL, y = "BC-from-T0 (AUC)") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.margin = margin(25, 5.5, 5.5, 18))

# ============================================================
# STATISTICS (printed to console and saved to output/)
# ============================================================
stats_log <- capture.output({
  # (1) Cage-robust test of vaccine effect: AUC of within-subject BC-from-baseline
  cat("\n--- AUC of BC-from-T0, Kruskal-Wallis across 3 groups ---\n")
  print(kruskal.test(auc ~ `Vaccine Administered`, data = auc_df))
  cat("Group AUC mean +/- sd:\n")
  print(auc_df %>% group_by(`Vaccine Administered`) %>%
          summarise(mean = mean(auc), sd = sd(auc), n = n(), .groups = "drop"))
  cat("Pairwise Mann-Whitney (RBC vs each treatment):\n")
  print(pairwise.wilcox.test(auc_df$auc, auc_df$`Vaccine Administered`,
                             p.adjust.method = "none"))

  # (2) PERMANOVA on genus-level Bray-Curtis: treatment, time, interaction
  #     NOTE: cage is confounded with treatment; Day 0 (pre-injection) groups
  #     already separate, so the "treatment" term reflects baseline cage effects,
  #     not vaccination. Time is the dominant driver.
  bc_all <- vegdist(genus_rel[sample_ids$seqid, ], method = "bray")
  perm_meta <- sample_ids
  set.seed(42)
  cat("\n--- PERMANOVA (genus Bray-Curtis) ---\n")
  print(adonis2(bc_all ~ `Vaccine Administered` * Day, data = perm_meta,
                permutations = 4999, by = "terms"))
  cat("\n--- Per-day PERMANOVA (treatment), incl. Day 0 pre-injection ---\n")
  for (d in sort(unique(perm_meta$Day))) {
    ids <- perm_meta$seqid[perm_meta$Day == d]
    set.seed(42)
    a <- adonis2(vegdist(genus_rel[ids, ], "bray") ~
                   perm_meta$`Vaccine Administered`[match(ids, perm_meta$seqid)],
                 permutations = 4999)
    cat(sprintf("Day %2d: R2=%.3f, p=%.4f\n", d, a$R2[1], a$`Pr(>F)`[1]))
  }
})
cat(stats_log, sep = "\n")
writeLines(stats_log, file.path(out_dir, "Supplementary_microbiomics_stats.txt"))
write_csv(auc_df, file.path(out_dir, "bc_from_T0_auc_per_mouse.csv"))
write_csv(bc_from_baseline, file.path(out_dir, "bc_from_T0_per_sample.csv"))

# ============================================================
# COMBINE and save to figs/
# ============================================================
row1 <- plot_grid(pA, labels = "A")
row2 <- plot_grid(pB, pC, NULL, labels = c("B", "C", ""), ncol = 3, rel_widths = c(1, 0.3, 0.5))
combined <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 0.9))

ggsave(file.path(fig_dir, "Supplementary_microbiomics.pdf"),
       combined, width = 11, height = 6)
ggsave(file.path(fig_dir, "Supplementary_microbiomics.png"),
       combined, width = 11, height = 6, dpi = 300, bg = "white")
