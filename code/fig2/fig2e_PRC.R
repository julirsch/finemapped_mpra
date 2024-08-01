# Load libraries
library(tidyverse)
library(vroom)
library(BuenColors)
library(binom)
library(cowplot)
library(patchwork)
source("code/utils.R")
source("code/fig2/PRC_functions.R")
options(stringsAsFactors = FALSE)

# Read in processed MPRA results
eqtl_mpra_df <- vroom("data/preprocess/mpra_gtex.txt.gz")
traits_mpra_df <- vroom("data/preprocess/mpra_complex_traits.txt.gz") 

# Add additional annotation columns
eqtl_mpra_df <- eqtl_mpra_df %>% 
  add_annots("CRE_yes_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 1, CRE == 0), .) %>%
  add_annots("CRE_no_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 0, CRE == 1), .) %>%
  add_annots("CRE_yes_active_any_yes", rlang::exprs(active_any == 1, active_any == 0, CRE == 1, CRE == 0), .)
traits_mpra_df <- traits_mpra_df %>% 
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
cowplot::save_plot('figures/fig2/fig2e_prc.pdf',
                   plt_combined,
                   base_height = 6,
                   base_width = 6,
                   device = cairo_pdf
)

