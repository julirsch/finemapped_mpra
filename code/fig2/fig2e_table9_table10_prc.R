# Load libraries
library(tidyverse)
library(vroom)
library(BuenColors)
library(binom)
library(cowplot)
library(patchwork)
source("code/utils.R")
source("code/github_ready/prc.functions.R")
options(stringsAsFactors = FALSE)

# Read in processed MPRA results
mpra_df <- vroom("/mnt/sdb/gwas_eqtl_mpra/data/preprocess/core_mpra.txt.gz") %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0))

# Add additional annotation columns
eqtl_mpra_df <- mpra_df %>%
  dplyr::filter(cohort == 'GTEx' | grepl('tissue',type)) %>%
  add_annots("CRE_yes_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 1, CRE == 0), .) %>%
  add_annots("CRE_no_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 0, CRE == 1), .) %>%
  add_annots("CRE_yes_active_any_yes", rlang::exprs(active_any == 1, active_any == 0, CRE == 1, CRE == 0), .)
traits_mpra_df <- mpra_df %>% 
  dplyr::filter(cohort %in% c('UKBB',"BBJ") | type %in% c('loc_CS','loc_PIP10','annot_CS','annot_PIP10','null_CS','null_PIP10')) %>% 
  add_annots("CRE_yes_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 1, CRE == 0), .) %>%
  add_annots("CRE_no_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 0, CRE == 1), .) %>%
  add_annots("CRE_yes_active_any_yes", rlang::exprs(active_any == 1, active_any == 0, CRE == 1, CRE == 0), .)

# Make df of causal variants for PrC
eqtl_prset_df <- make_prset(eqtl_mpra_df, "eQTLs")
traits_prset_df <-  make_prset(traits_mpra_df, "complex_traits")

# Count number of causal variants in data frame
bind_rows(eqtl_prset_df %>% 
            dplyr::filter(causal == T) %>% 
            dplyr::select(variant),
          traits_prset_df %>% 
            dplyr::filter(causal == T) %>% 
            dplyr::select(variant)) %>%
  distinct(variant)


# Save off supplementary table 9
mpra_mc_save <- eqtl_prset_df %>% dplyr::select(c("variant",
                                                  "pip",
                                                  "causal",
                                                  "emVar_any", 
                                                  "CRE",
                                                  "CRE_yes_emVar_any_yes",
                                                  "CRE_no_emVar_any_yes",
                                                  'active_any', 
                                                  "CRE_yes_active_any_yes")) %>%
  distinct() %>%
  dplyr::mutate(var_type = ifelse(causal==T, "eQTL test","eQTL control"))

pip9_save <- traits_prset_df %>% dplyr::select(c("variant",
                                                 "pip",
                                                 "causal",
                                                 "emVar_any", 
                                                 "CRE",
                                                 "CRE_yes_emVar_any_yes",
                                                 "CRE_no_emVar_any_yes",
                                                 'active_any', 
                                                 "CRE_yes_active_any_yes")) %>%
  distinct() %>%
  dplyr::mutate(var_type = ifelse(causal==T, "Complex traits test","Complex traits control")) 

test <- bind_rows(mpra_mc_save,pip9_save)
test %>% vroom::vroom_write(file = "tables/stable9.txt", delim = "\t", col_names = T)


# Calculate precision and recall
annotlist <- c("emVar_any", "CRE", "active_any", "CRE_yes_emVar_any_yes", "CRE_no_emVar_any_yes")
eqtl_out <- lapply(annotlist, function(x) {prec_rec(x, eqtl_prset_df)}) %>%
  bind_rows() %>%
  dplyr::mutate(cohort = "eQTL", .before = annot)
traits_out <- lapply(annotlist, function(x) {prec_rec(x, traits_prset_df)}) %>%
  bind_rows() %>%
  dplyr::mutate(cohort = "Complex traits", .before = annot)
out <- bind_rows(eqtl_out, traits_out)

# Plot basic PrC (Fig. 2e)
p1 <- ggplot(data = out, aes(x = recall, y = precision, color = annot, shape = cohort)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = prec_lower, ymax = prec_upper), width = 0) +
  geom_errorbarh(aes(xmin = rec_lower, xmax = rec_upper), height = 0) + 
  BuenColors::pretty_plot(fontsize = 20) +
  theme(legend.position = "none", 
        aspect.ratio = 1,
        panel.border = element_blank(),
        axis.line = element_line()) + 
  scale_color_manual(values = BuenColors::jdb_palette("corona")[c(4, 2, 5, 3, 1)]) +
  geom_hline(yintercept = 0.5, color = "grey50", linetype = "dashed")

# Save plot (Fig. 2e)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(6))

plt_combined
cowplot::save_plot('figures/fig2e.pdf',
                   plt_combined,
                   base_height = 6,
                   base_width = 6,
                   device = cairo_pdf
)

# save table (part of Supp Table 10)
out %>%
  dplyr::select(annot, cohort, precision, recall, TP, FP, FN, TN, total_pos, total, FPR, N, prec_upper, prec_lower, rec_upper,rec_lower) %>%
  dplyr::mutate(panel = 'Fig. 2e',
                description = 'Precision and recall') %>%
  vroom::vroom_write('tables/stable10_2e.txt')
