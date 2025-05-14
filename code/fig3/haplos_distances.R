
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

# read in mpra
haplos_new_df <- vroom('data/preprocess/haplos/haplos_master_table.txt.gz')
ints <- vroom('data/preprocess/haplos/haplos_ints.txt.gz')

# fig 3c
pc<- haplos_new_df %>%
  group_by(v1,v2) %>% 
  # require testing in all 3 windows to not skew for longer pairs designed
  dplyr::filter(n_distinct(window) == 3) %>%
  ungroup() %>% 
  # among things that are emVars, are interacting emVars closer together?
  dplyr::filter(any_emVar_meta == T) %>%
  distinct(v1,v2,distance, int_emVar_meta_new) %>%
  group_by(v1,v2) %>%
  dplyr::mutate(int_emVar_meta_any = any(int_emVar_meta_new==T)) %>%
  ungroup() %>%
  distinct(v1,v2,distance,int_emVar_meta_any) %>%
  ggplot(aes(x = distance, group = int_emVar_meta_any, fill = int_emVar_meta_any, color = int_emVar_meta_any, alpha = 0.1)) +
  geom_boxplot()+
  scale_fill_manual(values= c("grey","maroon"))+
  scale_color_manual(values = c("grey","maroon"))+
  pretty_plot()+
  theme(legend.position = "none")

pc

test <-haplos_new_df %>%
  group_by(v1,v2) %>% 
  dplyr::filter(n_distinct(window) == 3) %>%
  ungroup() %>% 
  dplyr::filter(any_emVar_meta == T) %>%
  distinct(v1,v2,distance, int_emVar_meta_new) %>%
  group_by(v1,v2) %>%
  dplyr::mutate(int_emVar_meta_any = any(int_emVar_meta_new==T)) %>%
  ungroup() %>%
  distinct(v1,v2,distance,int_emVar_meta_any) 

test
t.test(test %>% dplyr::filter(int_emVar_meta_any ==T) %>% .$distance, 
       test %>% dplyr::filter(int_emVar_meta_any ==F) %>% .$distance,
       alternative = "two.sided", var.equal = FALSE)

cowplot::save_plot(
  'figures/fig3/fig3c.pdf',
  pc,
  base_height = 1,
  base_width = 5,
  device = cairo_pdf)


test <- haplos_new_df %>% 
  dplyr::inner_join(ints %>% distinct(v1,v2,int_type) ) %>% 
  distinct(v1,v2,distance,window, int_type) %>% 
  group_by(v1,v2) %>% 
  dplyr::filter(n_distinct(window) == 3) %>%
  ungroup() %>% 
  distinct(v1,v2,distance, int_type)

test
t.test(test %>% dplyr::filter(int_type =='amplifying') %>% .$distance, 
       test %>% dplyr::filter(int_type =='dampening') %>% .$distance,
       alternative = "two.sided", var.equal = FALSE)
