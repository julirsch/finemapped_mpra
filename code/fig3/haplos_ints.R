# produces 'data/preprocess/haplos/haplos_ints.txt.gz'
# fig 3d heatmap, 3b point/line plot
# supplementary figure for fig3 num_int_emVars_type.pdf
unixtools::set.tempdir("/mnt/sdb/rstudio/sirajl/")
setwd('code/figures/fig3')
source("code/utils.R")
# Load libraries
library(tidyverse)
library(BuenColors)
library(binom)
library(ggplot2)
library(cowplot)
library(patchwork)
library(vroom)
library("PNWColors")
library(ggpointdensity)
library(ggpubr)
library(ggrastr)
library(circlize)
library(ComplexHeatmap)
library(metafor)
library(GenomicRanges)
options(stringsAsFactors = FALSE)

# read in mpra
mpra <- vroom('ata/preprocess/core_mpra.txt.gz')
haplos_new_df <- vroom('data/preprocess/haplos/haplos_master_table.txt.gz')
out_meta_info <- vroom('data/preprocess/haplos/haplos_CS_annots.txt.gz')

ints<- haplos_new_df %>% 
  group_by(v1,v2) %>% 
  dplyr::mutate(int_emVar_meta_any = any(int_emVar_meta_new)) %>% 
  ungroup() %>% 
  dplyr::filter(int_emVar_meta_new == T) %>%
  dplyr::mutate(indiv_effects = log2Skew_v1incv2dec + log2Skew_v1decv2inc) %>%
  dplyr::mutate(diff = log2Skew_v1incv2inc - indiv_effects,
                way = ifelse(diff > 0 ,TRUE, FALSE)) %>%
  #distinct(v1,v2,cell_type,log2Skew_v1decv2dec, log2Skew_v1incv2inc,int_log2Skew_meta_padj, log2Skew_v1incv2dec,log2Skew_v1decv2inc, indiv_effects, diff, pval) %>%
  group_by(v1,v2) %>% 
  dplyr::filter(pval == max(pval)) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup()

ints <- ints %>%
  dplyr::mutate(int_type = ifelse(diff < 0, "dampening", "amplifying"))
ints <- ints %>%
  dplyr::mutate(v1v2 = paste(v1,v2, sep = ';'))


ints %>% dplyr::filter(grepl('6:169720239', v2))

ints %>% vroom::vroom_write('data/preprocess/haplos/haplos_ints.txt.gz', delim = '\t')

ints %>% distinct(v1,v2)

col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
a <- Heatmap(ints %>%
          dplyr::select(log2Skew_v1incv2dec, log2Skew_v1decv2inc, indiv_effects, log2Skew_v1incv2inc, diff),
        cluster_columns = FALSE,
        col = col_fun) 
a
### save manually!!! as /mnt/sdb/gwas_eqtl_mpra/figures/fig3/fig3_heatmap.pdf

p <- left_join(ints, out_meta_info %>% distinct(v1, v2, cohort, any_same_CS_cohort), by = c('v1','v2')) %>%
  dplyr::count(cohort, any_same_CS_cohort, int_type) %>%
  na.omit() %>% 
  dplyr::mutate(group = paste(cohort, any_same_CS_cohort, sep = ' ')) %>% 
  ggplot(., aes(x=group, y=n, fill = int_type)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c('red','blue'))+
  xlab('same CS') + 
  pretty_plot()+
  theme(legend.position.inside = c(0.9,0.9)) 

p

cowplot::save_plot(
  'figures/fig3/s12l_num_int_emVars_type.pdf',
  p,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf)

# Set up plotting theme
my_theme <-
  BuenColors::pretty_plot(fontsize = 12) +
  BuenColors::L_border() +
  theme(
    plot.background = element_blank(),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"),
    plot.title = element_text(hjust = 4e-3, margin = margin(b = -12)),
    # legend.position = "none",
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.title = element_text(margin = margin(0.0, 0.0, 0.0, 0.0)),
    legend.background = element_blank(),
    legend.key.size = unit(0.2, "cm")
  )

# identify which interacting emVars are chosen through cell-type-window-meta-analysis
ct_meta_ints <- haplos_new_df %>%
  distinct(v1,v2, cell_type, center_var, window,
           log2Skew_v1decv2inc,log2Skew_v1incv2dec, log2Skew_v1incv2inc,
           altref_log2Skew, refalt_log2Skew, altalt_log2Skew,
           altalt_log2Skew_meta, altref_log2Skew_meta, refalt_log2Skew_meta, refref_log2Skew_meta,
           int_emVar_meta_new,
           fdr, beta, n_oligo,
           fdr2, beta2,
           cts_fdr, cts_beta, cts_n_oligo,
           cts_fdr2, cts_beta2,
           any_active_update, int_emVar_meta_new, lowest_source) %>% 
  dplyr::filter(int_emVar_meta_new == T) %>% 
  dplyr::mutate(diff = log2Skew_v1incv2dec + log2Skew_v1decv2inc - log2Skew_v1incv2inc) %>% 
  dplyr::select(-starts_with('log2Skew_v')) %>% 
  dplyr::mutate(window.meta.1 = (fdr > -log10(0.05) & (beta > 0.25 | beta < -0.25) & n_oligo > 1),
                window.meta.2 = (fdr2 > -log10(0.05) & (beta2 > 0.25 | beta2 < -0.25) & n_oligo > 1),
                ct.meta.1 = (cts_fdr > -log10(0.05) & (cts_beta > 0.25 | cts_beta < -0.25) & cts_n_oligo > 1),
                ct.meta.2 = (cts_fdr2 > -log10(0.05) & (cts_beta2 > 0.25 | cts_beta2 < -0.25)) & cts_n_oligo > 1) %>%
  dplyr::filter(window.meta.1 == F, window.meta.2 == F) %>% 
  distinct(v1, v2, int_emVar_meta_new, ct.meta.1, ct.meta.2)

# and get the cell-type-window-meta-analysis additive and double skews
meta_ct_haplos <- vroom('data/preprocess/haplo/meta_windows_celltypes.txt.gz') %>%
  dplyr::mutate(v2 = gsub("chr:", "chr", v2)) %>% 
  dplyr::mutate(refref_log2Skew_meta = 0) %>%
  # get the lowest log2Skew
  dplyr::mutate(lowest = pmin(altalt_log2Skew_meta, altref_log2Skew_meta, refalt_log2Skew_meta, refref_log2Skew_meta)) %>%
  # denote which measurement was lowest
  dplyr::mutate(lowest_source = case_when(altalt_log2Skew_meta == lowest ~ 'aa',
                                          altref_log2Skew_meta == lowest ~ 'ar',
                                          refalt_log2Skew_meta == lowest ~ 'ra',
                                          refref_log2Skew_meta == lowest ~ 'rr')) %>%
  # shift everything by the lowest
  dplyr::mutate(altalt_log2Skew_recode = altalt_log2Skew_meta - lowest,
                refalt_log2Skew_recode = refalt_log2Skew_meta - lowest,
                altref_log2Skew_recode = altref_log2Skew_meta - lowest,
                refref_log2Skew_recode = refref_log2Skew_meta - lowest) %>%
  # reassign decdec, incdec, decinc, incinc
  dplyr::mutate(log2Skew_v1decv2dec = case_when(lowest_source == 'rr' ~ refref_log2Skew_recode,
                                                lowest_source == 'aa' ~ altalt_log2Skew_recode,
                                                lowest_source == 'ra' ~ refalt_log2Skew_recode,
                                                lowest_source == 'ar' ~ altref_log2Skew_recode),
                log2Skew_v1incv2dec = case_when(lowest_source == 'rr' ~ altref_log2Skew_recode,
                                                lowest_source == 'aa' ~ refalt_log2Skew_recode,
                                                lowest_source == 'ra' ~ altalt_log2Skew_recode,
                                                lowest_source == 'ar' ~ refref_log2Skew_recode),
                log2Skew_v1decv2inc = case_when(lowest_source == 'rr' ~ refalt_log2Skew_recode,
                                                lowest_source == 'aa' ~ altref_log2Skew_recode,
                                                lowest_source == 'ra' ~ refref_log2Skew_recode,
                                                lowest_source == 'ar' ~ altalt_log2Skew_recode),
                log2Skew_v1incv2inc = case_when(lowest_source == 'rr' ~ altalt_log2Skew_recode,
                                                lowest_source == 'aa' ~ refref_log2Skew_recode,
                                                lowest_source == 'ra' ~ altref_log2Skew_recode,
                                                lowest_source == 'ar' ~ refalt_log2Skew_recode)) %>%
  # get rid of intermediate columns
  dplyr::select(-altalt_log2Skew_recode,-refref_log2Skew_recode,-altref_log2Skew_recode,-refalt_log2Skew_recode, -lowest) %>%
  distinct(v1, v2 ,log2Skew_v1decv2inc,log2Skew_v1incv2dec, log2Skew_v1incv2inc, any_emVar_meta) %>%
  inner_join(ct_meta_ints, by = c('v1','v2')) 

# Plot fig 3b
plt <- haplos_new_df %>% 
  dplyr::filter(v1 %ni% ct_meta_ints$v1 & v2 %ni% ct_meta_ints$v2) %>% 
  distinct(v1,v2, cell_type ,log2Skew_v1decv2inc,log2Skew_v1incv2dec, log2Skew_v1incv2inc, any_emVar_meta) %>% 
  dplyr::filter(any_emVar_meta == T) %>%
  #dplyr::summarize(cor = cor(altref_log2Skew_meta + refalt_log2Skew_meta, altalt_log2Skew_meta))
  ggplot(aes(x = log2Skew_v1incv2dec + log2Skew_v1decv2inc,
             y = log2Skew_v1incv2inc)) +
  rasterize(geom_pointdensity(), dpi=1000) +
  scale_color_gradientn(colors = BuenColors::jdb_palette("solar_blues")) +
  xlab("additive single alternative alleles") + 
  ylab("double alternative alleles ") + 
  rasterize(geom_point(data = haplos_new_df %>% dplyr::filter(v1 %ni% ct_meta_ints$v1 & v2 %ni% ct_meta_ints$v2) %>% 
                         distinct(v1,v2, cell_type, log2Skew_v1decv2inc,log2Skew_v1incv2dec, log2Skew_v1incv2inc, int_emVar_meta_new) %>% 
                         dplyr::filter(int_emVar_meta_new == T), 
                       color =  BuenColors::jdb_palette("brewer_celsius")[9]), dpi=1000) +
  rasterize(geom_point(data = meta_ct_haplos, 
                       color =  BuenColors::jdb_palette("brewer_celsius")[8]), dpi=1000) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  pretty_plot() +
  ylim(0, 4) + xlim(0, 4) +
  my_theme +
  theme(legend.position = "none")






plt

cowplot::save_plot(
  'figures/fig3/fig3b_int_emvar_meta_plot_meta_cts_highlight.pdf',
  plt,
  base_height = 5,
  base_width = 5,
  device = cairo_pdf)

                  
