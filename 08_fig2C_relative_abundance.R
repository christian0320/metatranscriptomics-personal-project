#!/usr/bin/env Rscript

library(tidyverse)

PROJECT <- "/home/christianarturo/project1"

# === Input files ===
REL_FILE <- file.path(PROJECT, "data/processed/relative_abundance_species.csv")
TAX_FILE <- file.path(PROJECT, "data/taxonomy/taxonomic_information_bacteria.txt")
OUTDIR   <- file.path(PROJECT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ======================================================
# 1. Load relative abundance (species-level)
# ======================================================
rel <- read_csv(REL_FILE)

rel_long <- rel %>%
    pivot_longer(-species, names_to = "sample", values_to = "abundance")

# ======================================================
# 2. Load taxonomy (Species → Genus)
# ======================================================
taxonomy <- read_tsv(TAX_FILE) %>%
    select(Strain, Genus) %>%
    rename(species = Strain)

taxonomy$Genus[is.na(taxonomy$Genus)] <- "NA"

# ======================================================
# 3. Merge abundances with taxonomy
# ======================================================
merged <- rel_long %>%
    left_join(taxonomy, by = "species")

# ======================================================
# 4. Aggregate ABUNDANCE by genus
# ======================================================
rel_genus <- merged %>%
    group_by(Genus, sample) %>%
    summarise(total_abundance = sum(abundance), .groups = "drop")

# ======================================================
# 5. Convert to % relative abundance
# ======================================================
rel_genus <- rel_genus %>%
    group_by(sample) %>%
    mutate(percent = 100 * total_abundance / sum(total_abundance)) %>%
    ungroup()

# ======================================================
# 6. REMOVE genera that never appear (0% everywhere)
# ======================================================
rel_genus <- rel_genus %>%
    group_by(Genus) %>%
    mutate(genus_total = sum(percent)) %>%
    ungroup() %>%
    filter(genus_total > 0) %>%
    select(-genus_total)

# ======================================================
# 7. Rename sample IDs
# ======================================================
sample_names <- c(
  "SRR24611159" = "Matrix – Soil sample 1",
  "SRR24611160" = "Matrix – Soil sample 2",
  "SRR24611161" = "Matrix – Soil sample 3",
  "SRR24611162" = "Root – Root sample 1",
  "SRR24611163" = "Root – Root sample 2",
  "SRR24611164" = "Root – Root sample 3"
)

rel_genus$sample <- recode(rel_genus$sample, !!!sample_names)

# ======================================================
# 8. Reorder genera by total %
# ======================================================
genus_order <- rel_genus %>%
    group_by(Genus) %>%
    summarise(total = sum(percent)) %>%
    arrange(desc(total)) %>%
    pull(Genus)

rel_genus$Genus <- factor(rel_genus$Genus, levels = genus_order)

# ======================================================
# 9. PLOT (Fig 2C style)
# ======================================================
p <- rel_genus %>%
    ggplot(aes(x = sample, y = percent, fill = Genus)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    labs(
        title = "Relative abundance of bacterial strains (Fig 2C)",
        x = "Sample",
        y = "Relative Abundance (%)"
    )

ggsave(
    filename = file.path(OUTDIR, "Fig2C_relative_abundance.png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
)

cat("Fig 2C generated at:", file.path(OUTDIR, "Fig2C_relative_abundance.png"), "\n")
