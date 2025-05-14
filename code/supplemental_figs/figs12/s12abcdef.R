# look at correlations between windows
# fig 3 supplements a-f

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

# define colors
clrs = c(
  pnw_palette('Bay',8,type='continuous')[5],
  pnw_palette('Bay',8,type='continuous')[3],
  pnw_palette('Starfish')[5],
  pnw_palette('Bay',8,type='continuous')[2],
  pnw_palette('Bay',8,type='continuous')[8]
)


# plot skew

# Get correlations across windows
haplots_wincorr_df <- haplos_long_df %>%
  pivot_wider(names_from = "window",
              values_from = "value") %>%
  dplyr::group_by(measurement) %>%
  rstatix::cor_test("left", "right", "middle") %>% 
  dplyr::filter(measurement == 'log2Skew', var1 != var2)
a =haplos_long_df %>%
  pivot_wider(names_from = "window",
              values_from = "value") %>%
  dplyr::filter(measurement=='log2Skew') %>%
  ggplot(. ,aes(x=left,y=middle, group = cell_type, color = cell_type))+
  rasterize(geom_point(alpha = 0.5, size = 0.5), dpi = 1000)+
  pretty_plot()+
  scale_color_manual(values = clrs)+
  geom_smooth(method = "lm", fill = NA)+
  annotate("text", 
           label = paste0("R = ", as.character(round(haplots_wincorr_df %>% dplyr::filter(var1=='left',var2=='middle') %>% .$cor,3))), 
           x = 0.3, y = 2.5, size = 4, colour = "black")+
  theme(legend.position = 'None')+
  coord_fixed(ratio = 1)+
  xlim(-3,3)+
  ylim(-3,3)+
  ggtitle('Skew correlations')
a


b = haplos_long_df %>%
  pivot_wider(names_from = "window",
              values_from = "value") %>%
  dplyr::filter(measurement == 'log2Skew') %>%
  ggplot(. ,aes(x=right,y=middle, group = cell_type, color = cell_type))+
  rasterize(geom_point(alpha = 0.5, size = 0.5), dpi = 1000)+
  pretty_plot()+
  scale_color_manual(values = clrs)+
  geom_smooth(method = "lm", fill = NA)+
  annotate("text", 
           label = paste0("R = ", as.character(round(haplots_wincorr_df %>% dplyr::filter(var1=='right',var2=='middle') %>% .$cor,3))), 
           x = 0.3, y = 2.5, size = 4, colour = "black")+
  theme(legend.position = 'None')+
  coord_fixed(ratio = 1)+
  xlim(-3,3)+
  ylim(-3,3)+
  ggtitle('Skew correlations')
b

c =  haplos_long_df %>%
  pivot_wider(names_from = "window",
              values_from = "value") %>%
  dplyr::filter(measurement == 'log2Skew') %>%
  ggplot(. ,aes(x=left,y=right, group = cell_type, color = cell_type))+
  rasterize(geom_point(alpha = 0.5, size = 0.5), dpi = 1000)+
  pretty_plot()+
  scale_color_manual(values = clrs)+
  geom_smooth(method = "lm", fill = NA)+
  annotate("text", 
           label = paste0("R = ", as.character(round(haplots_wincorr_df %>% dplyr::filter(var1=='left',var2=='right') %>% .$cor,3))), 
           x = 0.3, y = 2.5, size = 4, colour = "black")+
  theme(legend.position = 'None')+
  coord_fixed(ratio = 1)+
  xlim(-3,3)+
  ylim(-3,3)+
  ggtitle('Skew correlations')
c

# save skew plots
cowplot::save_plot(
  paste0('figures/fig3/s12d.pdf'),
  a,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf)

cowplot::save_plot(
  paste0('figures/fig3/s12e.pdf'),
  b,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf)

cowplot::save_plot(
  paste0('figures/fig3/s12f.pdf'),
  c,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf)


# Get correlations across windows
# plot activity
haplots_wincorr_df <- haplos_long_df %>%
  pivot_wider(names_from = "window",
              values_from = "value") %>%
  dplyr::group_by(measurement) %>%
  rstatix::cor_test("left", "right", "middle") %>% 
  dplyr::filter(measurement == 'Log2FC', var1 != var2)

a =haplos_long_df %>%
  pivot_wider(names_from = "window",
              values_from = "value") %>%
  dplyr::filter(measurement=='Log2FC') %>%
  ggplot(. ,aes(x=left,y=middle, group = cell_type, color = cell_type))+
  rasterize(geom_point(alpha = 0.5, size = 0.5), dpi = 1000)+
  pretty_plot()+
  scale_color_manual(values = clrs)+
  geom_smooth(method = "lm", fill = NA)+
  annotate("text", 
           label = paste0("R = ", as.character(round(haplots_wincorr_df %>% dplyr::filter(var1=='left',var2=='middle') %>% .$cor,3))), 
           x = 0.3, y = 6, size = 4, colour = "black")+
  theme(legend.position = 'None')+
  coord_fixed(ratio = 1)+
  xlim(-2,7)+
  ylim(-2,7)+
  ggtitle('Activity correlations')
a

b = haplos_long_df %>%
  pivot_wider(names_from = "window",
              values_from = "value") %>%
  dplyr::filter(measurement == 'Log2FC') %>%
  ggplot(. ,aes(x=right,y=middle, group = cell_type, color = cell_type))+
  rasterize(geom_point(alpha = 0.5, size = 0.5), dpi = 1000)+
  pretty_plot()+
  scale_color_manual(values = clrs)+
  geom_smooth(method = "lm", fill = NA)+
  annotate("text", 
           label = paste0("R = ", as.character(round(haplots_wincorr_df %>% dplyr::filter(var1=='right',var2=='middle') %>% .$cor,3))), 
           x = 0.3, y = 6, size = 4, colour = "black")+
  theme(legend.position = 'None')+
  coord_fixed(ratio = 1)+
  xlim(-2,7)+
  ylim(-2,7)+
  ggtitle('Activity correlations')
b

c =  haplos_long_df %>%
  pivot_wider(names_from = "window",
              values_from = "value") %>%
  dplyr::filter(measurement == 'Log2FC') %>%
  ggplot(. ,aes(x=left,y=right, group = cell_type, color = cell_type))+
  rasterize(geom_point(alpha = 0.5, size = 0.5), dpi = 1000)+
  pretty_plot()+
  scale_color_manual(values = clrs)+
  geom_smooth(method = "lm", fill = NA)+
  annotate("text", 
           label = paste0("R = ", as.character(round(haplots_wincorr_df %>% dplyr::filter(var1=='left',var2=='right') %>% .$cor,3))), 
           x = 0.3, y = 6, size = 4, colour = "black")+
  theme(legend.position = 'None')+
  coord_fixed(ratio = 1)+
  xlim(-2,7)+
  ylim(-2,7)+
  ggtitle('Activity correlations')
c

cowplot::save_plot(
  paste0('figures/fig3/s12.pdf'),
  a,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf)

cowplot::save_plot(
  paste0('figures/fig3/s12b.pdf'),
  b,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf)

cowplot::save_plot(
  paste0('figures/fig3/s12c.pdf'),
  c,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf)

