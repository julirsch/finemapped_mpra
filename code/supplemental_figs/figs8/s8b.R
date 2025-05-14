
# Load libraries
setwd('/mnt/sdb/gwas_eqtl_mpra/reviews/code/')
source("/mnt/sdb/gwas_eqtl_mpra/code/utils.R")
library(plyranges)
library(tidyverse)
library(DESeq2)
library(GenomicRanges)
library(BuenColors)
library(rtracklayer)
library(reshape2)
library(binom)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE)
library(cowplot)
library(patchwork)
library("PNWColors")


mpra_df <- vroom('/mnt/sdb/gwas_eqtl_mpra/data/preprocess/core_mpra.txt.gz') %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0))

test <-mpra_df %>% 
  ungroup() %>% 
  dplyr::filter(cohort %in% c("UKBB","BBJ")) %>% 
  distinct(variant, trait, pip, emVar, cell_type) %>%
  dplyr::filter(pip > 0.1)  %>%
  distinct(variant, trait, emVar, cell_type) %>%
  group_by(trait, cell_type) %>%
  dplyr::count(trait, cell_type, emVar) %>%
  dplyr::mutate(ratio = ifelse(n_distinct(emVar == 2), n[emVar==T] /(n[emVar==F] + n[emVar == T]), 0))

zeros <- test %>% dplyr::filter(n_distinct(emVar)==1) %>%
  dplyr::mutate(emVar=T) %>%
  dplyr::mutate(n = 0) %>%
  dplyr::mutate(ratio = 0)

props_nums <- bind_rows(test %>% dplyr::filter(emVar==T) %>% ungroup(), zeros)


# read in domain data
domain <- c('Metabolic','Lipids','Cardiovascular','Immunological','Hematopoietic','Hepatic','Renal','Skeletal','Neurological','Behavioral','Other')
colors <- c('#4C72AE','#E4812F','#F3DE67','#6BBCCC','#5AA245','#BA2E2C','#8C61B6','#80584E','#CE72BB','#BACD3C','#7F7F7F')
colors.df <- cbind(domain, colors) %>%
  as_tibble()
traits.df <- read_delim("/mnt/sdb/gtex_mpra/release/UKBB_96traits_release1.traits", delim = "\t", col_names = T) %>%
  dplyr::mutate(prevalence = n_cases / n) %>%
  merge(colors.df, by = "domain") %>%
  as_tibble()
traits <- traits.df$trait
domain.df <- traits.df %>% 
  dplyr::select(trait, domain, colors)

props_nums <- props_nums %>% left_join(domain.df %>% distinct(trait,domain), by = 'trait')

props_nums$trait <- factor(props_nums$trait, levels=unique(props_nums$trait[order(props_nums$domain,props_nums$trait)]), ordered=TRUE)

p <-  ggplot(props_nums, aes(cell_type, forcats::fct_rev(trait), fill = ratio, size =n)) +
  geom_point(shape = 21, stroke = 0) +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(1, 15), breaks = c(25, 50, 100, 200 )) +
  scale_fill_gradient2(high = BuenColors::jdb_palette('brewer_spectra')[1], low = "#F9F5FF",breaks = c(0.00, 0.075, 0.15)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25),
                             label.position = "bottom",
                             title.position = "right",
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  labs(size = "Area = Prop emVars", fill = "Num emVars:", x = NULL, y = NULL)

p
plt_combined <- p

cowplot::save_plot(
  'figures/s8b.pdf',
  plt_combined,
  base_height = 5.8,
  base_width = 6,
  device = cairo_pdf)

props_nums



