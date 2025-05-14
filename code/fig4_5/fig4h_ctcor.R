# Load libraries
library(tidyverse)
library(vroom)
library(corrplot)

options(stringsAsFactors = FALSE)

code_path<-""
source(paste0(code_path,"code/Preprocess/utils.R"))
source(paste0(code_path,"code/Preprocess/theme.R"))

##Preprocessing in Preprocess/satmut.R
# Compute nucleotide level correlation across cell-types
out_df <- satmut_df %>%
  distinct(ID, sat_ref, cell_type, pos, mut_base, log2Skew, log2FC_baseline) %>%
  na.omit() %>%
  pivot_wider(names_from = c("cell_type"), values_from = c("log2Skew", "log2FC_baseline")) %>%
  na.omit() %>%
  dplyr::group_by(sat_ref) %>%
  dplyr::summarize(cor = Hmisc::rcorr(log2Skew_HEPG2, log2Skew_K562, type = "pearson")$r[2],
                   pval = cor.test(log2Skew_HEPG2, log2Skew_K562, type = "pearson")$p.value,
                   log2FC_baseline_HEPG2 = max(log2FC_baseline_HEPG2),
                   log2FC_baseline_K562 = max(log2FC_baseline_K562)) %>%
  arrange(-cor)

# Plot correlation across cell-types (Fig. 4h)
p1 <- out_df %>%
  ggplot(aes(x = log2FC_baseline_HEPG2, y = log2FC_baseline_K562, fill = cor)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradientn(colors = COL1("YlOrBr", n = 200)) +
  pretty_plot() +
  plot_theme

# Plot plot (Fig. 4h)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(2))
cowplot::save_plot(
  paste0("figures/fig4/fig4h_satmut_ctscor.pdf"),
  plt_combined,
  base_height = 2,
  base_width = 2.2,
  device = cairo_pdf
)

# Proportion with high correlation
out_df %>%
  dplyr::count(cor > 0.6) %>% 
  dplyr::mutate(denom = sum(n),
                perc = n / denom)
