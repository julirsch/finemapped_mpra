# look at correlations between windows
# fig 12 supplements h,i

# Load libraries
source("code/utils.R")
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
library("RColorBrewer")
options(stringsAsFactors = FALSE)

# Read in CRE-only haplotypes
haplos_df <- vroom('data/preprocess/haplos/haplos_master_table.txt.gz')

# Define and filter to best observation for each haplotype
haplos_df <- haplos_df %>%
  dplyr::group_by(v1,v2,window,center_var, cell_type)%>%
  dplyr::mutate(best_library = ifelse((sqrt(mean_Plasmid_altalt) + sqrt(mean_Plasmid_refref)+sqrt(mean_Plasmid_altref) + sqrt(mean_Plasmid_refalt)) ==
                                        max(sqrt(mean_Plasmid_altalt) + sqrt(mean_Plasmid_refref)+sqrt(mean_Plasmid_altref) + sqrt(mean_Plasmid_refalt)),1,0)) %>%
  dplyr::filter(best_library == 1) %>%
  dplyr::ungroup()

# Wide to long dataframe
haplos_long_df <- haplos_df %>%
  distinct(v1, v2, center_var, window, cell_type, 
           altalt_log2Skew, altref_log2Skew, refalt_log2Skew,
           altalt_Log2FC, altref_Log2FC,refalt_Log2FC) %>%
  pivot_longer(cols = c('altalt_log2Skew','altref_log2Skew','refalt_log2Skew',
                        'altalt_Log2FC', 'altref_Log2FC','refalt_Log2FC'),
               names_to = 'name',
               values_to = 'value') %>%
  tidyr::separate(col = 'name',
                  into = c('haplo','measurement'),
                  sep = '_') 

# Get correlations across cell-types
haplots_ctcorr_df <- haplos_long_df %>%
  pivot_wider(names_from = "cell_type",
              values_from = "value") %>%
  dplyr::group_by(measurement) %>%
  rstatix::cor_test("K562", "HCT116", "HEPG2", "SKNSH", "A549")

# Get correlations across windows
haplots_wincorr_df <- haplos_long_df %>%
  pivot_wider(names_from = "window",
              values_from = "value") %>%
  dplyr::group_by(measurement) %>%
  rstatix::cor_test("left", "right", "middle")


p1 <- haplots_ctcorr_df %>% 
  dplyr::mutate(
    var1 = factor(var1, levels = c('A549','HCT116','HEPG2','K562','SKNSH')),
    var2 = factor(var2, levels = c('SKNSH','K562','HEPG2','HCT116','A549'))
  ) %>% 
  ggplot(., aes(x=var1, var2, fill= cor)) +
  geom_tile() +
  geom_text(aes(label = cor))+
  #scale_fill_gradient2(low = COL1("YlOrBr", n = 200)[200], mid = COL1("YlOrBr", n = 200)[100], high = COL1("YlOrBr", n = 200)[1], midpoint = 0.75, limits = c(0, 1), oob = scales::squish)+
  scale_fill_distiller(type = 'seq', palette = 17, limits = c(0.4, 0.95), direction = 1, oob = scales::squish)+
  theme_minimal() +
  theme(panel.grid = element_blank())+
  ggtitle('Correlations across cell types')+
  ylab('')+
  xlab('')+
  facet_grid(~measurement)+
  coord_equal()
p1 

p2 <- haplots_wincorr_df %>% 
  dplyr::mutate(
    var1 = factor(var1, levels = c('left','middle','right')),
    var2 = factor(var2, levels = c('right','middle','left'))
  ) %>%
  ggplot(., aes(x=var1, var2, fill= cor)) +
  geom_tile() +
  geom_text(aes(label = cor))+
  #scale_fill_gradient2(low = COL1("YlOrBr", n = 200)[200], mid = COL1("YlOrBr", n = 200)[100], high = COL1("YlOrBr", n = 200)[1], midpoint = 0.75, limits = c(0, 1), oob = scales::squish)+
  scale_fill_distiller(type = 'seq', palette = 17, limits = c(0.4, 0.95), direction = 1, oob = scales::squish)+
  theme_minimal() +
  theme(panel.grid = element_blank())+
  ggtitle('Correlations across windows')+
  ylab('')+
  xlab('')+
  facet_grid(~measurement)+
  coord_equal()

p2

cowplot::save_plot(
  paste0('figures/fig3/s12h.pdf'),
  p1,
  base_height = 4,
  base_width = 8,
  device = cairo_pdf)

cowplot::save_plot(
  paste0('figures/fig3/s12i.pdf'),
  p2,
  base_height = 4,
  base_width = 8,
  device = cairo_pdf)

