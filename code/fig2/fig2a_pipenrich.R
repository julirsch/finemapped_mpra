# Load libraries
library(tidyverse)
library(vroom)
library(GenomicRanges)
library(rtracklayer)
library(plyranges)
library(matrixStats)
library(annotables)
library(BuenColors)
library(PNWColors)
library(cowplot)
library(patchwork)
source("code/utils.R")
source("code/theme.R")
options(stringsAsFactors = FALSE)

# Read in processed MPRA results
mpra_df <- vroom::vroom("data/preprocess/core_mpra.txt.gz")
mpra_meta_df <- vroom::vroom("data/preprocess/mpra_meta.txt.gz")
mpra_df <- mpra_df %>%
  left_join(mpra_meta_df,
            by = c("variant", "cohort")) %>%
  ungroup() 

# Filter to core trait and control variants without overlaps
mpra_df <- mpra_df %>% 
  dplyr::filter(cohort == "control" | (pip > 0.1 | cs_id > 0)) %>%
  group_by(variant) %>%
  dplyr::filter(n_distinct(cohort) == 1 | all(cohort %in% c("UKBB", "BBJ"))) %>% 
  ungroup()

# Filter to non-coding trait variants
mpra_pre <- mpra_df %>% 
  dplyr::filter(cohort != "control") %>%
  dplyr::filter(consequence %ni% c("LoF", "missense", "synonymous")) %>%
  group_by(variant, cohort) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  ungroup() %>%
  distinct(variant, cohort, pip, cell_type, emVar, emVar_any, log2Skew, Skew_logPadj, consequence) %>%
  dplyr::mutate(pip_bin = cut(pip, pip_bin_breaks, include.lowest = T, ordered_result = T),
                name = pip_bin) %>%
  dplyr::mutate(cohort = case_when(cohort %in% c("UKBB", "BBJ") ~ "Complex traits",
                                   cohort %in% c("GTEx") ~ "eQTLs",
                                   T ~ NA))

# Filter to non-coding control variants
mpra_pre_ctrls <- mpra_df %>% 
  dplyr::filter(cohort == "control") %>%
  dplyr::filter(consequence %ni% c("LoF", "missense", "synonymous")) %>%
  distinct(variant, cohort, cell_type, emVar, emVar_any, log2Skew, Skew_logPadj, type) %>% 
  dplyr::mutate(cohort = case_when(type %in% c("loc_CS", "loc_PIP10", "annot_PIP10", "annot_CS", "null_PIP10", "null_CS") ~ "Complex traits",
                                   type %in% c("3tissue_locctrl", "49tissue_locctrl", "3tissue_annotctrl", "49tissue_annotctrl") ~ "eQTLs")) %>% 
  dplyr::mutate(name = type)

# Re-combine and keep only variants measured in all cell-types for that group
mpra_comb_df <- bind_rows(mpra_pre, mpra_pre_ctrls) %>%
  dplyr::select(-pip, -type) %>%
  distinct() %>%
  group_by(variant, cohort, name) %>%
  dplyr::filter(n_distinct(cell_type) == 4) %>%
  arrange(-emVar, -Skew_logPadj) %>% 
  filter(row_number() == 1) %>%
  ungroup() %>%
  distinct()

# Get proportion emVar by group
mpra_prop_df <- mpra_comb_df %>%
  dplyr::count(cohort, name, emVar) %>%
  group_by(cohort, name) %>%
  dplyr::mutate(tot = sum(n)) %>%
  dplyr::filter(emVar == T) %>%
  dplyr::select(-emVar) %>% 
  dplyr::mutate(prop = n / tot,
                lower = binom.confint(n, tot, method = "wilson")$lower,
                upper = binom.confint(n, tot, method = "wilson")$upper) %>%
  arrange(-prop) %>%
  ungroup()

# Rename
mpra_prop_df <- mpra_prop_df %>% 
  dplyr::mutate(name = ifelse(name == "null_PIP10", "Null", name),
                name = ifelse(name %in% c("annot_PIP10", "49tissue_annotctrl"), "annot-matched", name),
                name = ifelse(name %in% c("loc_PIP10", "49tissue_locctrl"), "loc-matched", name)) %>%
  dplyr::filter(name %ni% c("null_CS", "annot_CS", "loc_CS", "3tissue_annotctrl", "3tissue_locctrl"))

# Hypothesis tests for emVar enrichment
# High-PIP complex traits vs Null
mpra_prop_df %>%
  dplyr::filter(cohort == "Complex traits" & name == "(0.9,1]" | cohort == "Complex traits" & name == "Null") %>%
  dplyr::summarize(pval = prop.test(n, tot)$p.value)
# High-PIP eQTLs vs Null
mpra_prop_df %>%
  dplyr::filter(cohort == "eQTLs" & name == "(0.9,1]" | cohort == "Complex traits" & name == "Null") %>%
  dplyr::summarize(pval = prop.test(n, tot)$p.value)

# Plot proportion of emVars by group (Fig. 2a)
p1 <- mpra_prop_df %>%
  dplyr::mutate(name = factor(name, levels = c("Null", "annot-matched", "loc-matched", "[0,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]", "(0.9,1]"))) %>%
  ggplot(aes(x = name, y = prop, group = cohort, fill = cohort)) + 
  geom_col(position = position_dodge(0.9, preserve = "single")) + 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                linewidth = 0.5, width = 0, col = "black",
                position = position_dodge(0.9, preserve = "single")) +
  pretty_plot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")) +
  ylab("Proportion of variants") +
  xlab("")
p1

# Save plot (Fig. 2a)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot(
  paste0("figures/fig2/fig2a_propemVar.pdf"),
  p1,
  base_height = 3,
  base_width = 3,
  device = cairo_pdf
)


