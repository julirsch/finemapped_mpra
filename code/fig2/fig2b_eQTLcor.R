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
mpra_df <- vroom("data/preprocess/core_mpra.txt.gz")
mpra_meta_df <- vroom("data/preprocess/mpra_meta.txt.gz")
mpra_df <- mpra_df %>%
  left_join(mpra_meta_df,
            by = c("variant", "cohort")) %>%
  ungroup()

# Read in TSS distance
gtex_tss <- vroom::vroom("data/tss/GTEx_tss.txt.gz")

# Filter to GTEx 95% CS emVars
# Select highest PIP, eQTL beta for most significant tissue across all genes
# Carefully call distal CRE and proximal CRE
out_df <- mpra_df %>%
  dplyr::filter(cohort %in% c("GTEx"), cs_id > 0) %>%
  dplyr::filter(emVar_any == T) %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(pip == max(pip), beta_marginal = mean(beta_marginal[abs(z) == max(abs(z))])) %>%
  dplyr::filter(row_number() == 1) %>%
  distinct() %>% 
  inner_join(gtex_tss,
             by = c("gene" = "gene_id", "variant_hg38" = "variant_id")) %>%
  dplyr::mutate(pip_bin = cut(pip, c(0, 0.1, 0.5, 1))) %>%
  dplyr::mutate(tss_bin = cut(abs(tss_distance), c(0, 1000, 5000, 20000, 50000, 1000000), include.lowest = T),
                promoter_dist = ifelse(tss_distance >  5000, 1, 0)) %>%
  dplyr::select(variant, pip_bin, tss_bin, tss_distance, log2Skew_meta, beta_marginal, promoter, promoter_dist, CRE) %>%
  ungroup() %>% 
  dplyr::mutate(CRE_type = case_when(CRE > 0 & abs(tss_distance) > 10000 & promoter == 0 ~ "distal CRE",
                                     CRE > 0 & abs(tss_distance) < 2000 & promoter == 1 ~ "promoter",
                                     CRE == 0 & abs(tss_distance) > 10000 & promoter == 0 ~ "distal non-CRE")) %>%
  drop_na() 

# Compute (Spearman) correlations and statistical summaries
out_df <- out_df %>%
  group_by(CRE_type, pip_bin) %>%
  dplyr::mutate(log2Skew_meta_sag = ifelse(sign(log2Skew_meta) == sign(beta_marginal), log2Skew_meta, NA)) %>%
  dplyr::summarize(cor = cor(rank(log2Skew_meta), rank(beta_marginal), method = "pearson"),
                   lower = cor.test(rank(log2Skew_meta), rank(beta_marginal))$conf.int[1],
                   upper = cor.test(rank(log2Skew_meta), rank(beta_marginal))$conf.int[2],
                   pval = cor.test(rank(log2Skew_meta), rank(beta_marginal))$p.value[1],
                   cor_abs = cor(rank(abs(log2Skew_meta)), rank(abs(beta_marginal)), method = "pearson"),
                   lower_abs = cor.test(rank(abs(log2Skew_meta)), rank(abs(beta_marginal)), na.rm = T)$conf.int[1],
                   upper_abs = cor.test(rank(abs(log2Skew_meta)), rank(abs(beta_marginal)), na.rm = T)$conf.int[2],
                   pval_abs = cor.test(rank(abs(log2Skew_meta)), rank(abs(beta_marginal)))$p.value[1],
                   binom_pval = binom.test(c(sum(!is.na(log2Skew_meta_sag)), length(log2Skew_meta_sag)))$p.value) %>%
  pivot_longer(!c(CRE_type, pip_bin)) %>%
  dplyr::mutate(type = ifelse(grepl("abs", name), "abs", "raw"),
                name = gsub("_abs", "", name)) %>%
  dplyr::filter(!(CRE_type == "distal non-CRE" & type == "abs")) %>%
  dplyr::mutate(type = ifelse(CRE_type == "distal non-CRE", "adj", type),
                CRE_type = ifelse(CRE_type == "distal non-CRE", "distal CRE", CRE_type)) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  ungroup()
 
# Create "adjusted" correlation and prepare to plot 
out_df <- out_df %>%
  group_by(pip_bin) %>%
  dplyr::mutate(cor_cumsum = case_when(type == "adj" ~ cor[type == "raw" & CRE_type == "distal CRE"] - cor[type == "adj"],
                                       T ~ cor),
                lower_cumsum = case_when(type == "adj" ~ lower[type == "raw" & CRE_type == "distal CRE"] - lower[type == "adj"],
                                         T ~ lower),
                upper_cumsum = case_when(type == "adj" ~ upper[type == "raw" & CRE_type == "distal CRE"] - upper[type == "adj"],
                                         T ~ upper)) %>%
  ungroup() %>%
  group_by(CRE_type, pip_bin) %>%
  dplyr::mutate(cor = case_when(type == "raw" ~ cor[type == "raw"] - max(cor_cumsum[type != "raw"]),
                                T ~ cor)) %>%
  ungroup() %>%
  group_by(pip_bin) %>%
  dplyr::mutate(cor = case_when(type == "adj" ~ cor_cumsum[type == "adj"] - cor[type == "abs" & CRE_type == "distal CRE"],
                                T ~ cor)) %>%
  ungroup() %>% 
  dplyr::mutate(type = ordered(type, c("raw", "adj", "abs"))) %>%
  distinct(pip_bin, CRE_type, type, cor, cor_cumsum, lower_cumsum, upper_cumsum, pval, binom_pval)

# Combine p-values (distal CRE + promoter) with Fisher to report
out_df %>%
  dplyr::filter(pip_bin %in% c("(0.5,1]"))
# Pearson correlation absolute effect size p-value
metap::sumlog(c(2.00e-13, 9.53e-4))$p
# Binomial sign p-value
metap::sumlog(c(7.11e-19, 8.24e-14))$p

# Plot correlations for high-PIP variants (Fig. 2b)
p1 <- out_df %>%
  dplyr::filter(pip_bin %in% c("(0.5,1]")) %>%
  ggplot(aes(x = CRE_type, y = cor)) +
  geom_bar(aes(alpha = type), stat = "identity")+
  geom_linerange(aes(ymin = lower_cumsum, ymax = upper_cumsum), position = position_dodge(width = 0)) +
  pretty_plot() +
  #scale_fill_manual(values = jdb_palette("FantasticFox")) + 
  xlab("") + 
  ylab("Correlation between eQTL beta and MPRA allelic effect") + 
  plot_theme +
  theme(legend.position = "none")

# Save plot (Fig. 2b)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(2))
cowplot::save_plot(
  paste0("figures/fig2/fig2b_eQTLcorr.pdf"),
  plt_combined,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf
)

