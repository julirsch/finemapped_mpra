
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
library(VennDiagram)
options(stringsAsFactors = FALSE)

# Read in processed MPRA data
mpra_df <- vroom('data/preprocess/core_mpra.txt.gz') %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0))

# PIP bin breaks
pip_bin_breaks <- c(0, 0.01, 0.1, 0.5, 0.9, 1.0)

# Supp Fig 1a pip breakdown
bar_plot_pips <-  mpra_df %>% 
  dplyr::mutate(cohort = case_when(cohort == 'GTEx' ~ 'eQTLs',
                                   cohort %in% c("UKBB", "BBJ") ~ 'Complex Traits',
                                   cohort == 'control' ~ 'Controls')) %>% 
  distinct(variant, pip, cohort, type)

# separate out controls
bar_plot_pips_counts <- bar_plot_pips %>% dplyr::filter(cohort == 'Controls') %>%
  dplyr::mutate(pip_bin = case_when(type %in% c('loc_CS',
                                                'loc_PIP10') ~ 'Complex Traits: Location',
                                    type %in% c('annot_CS',
                                                'annot_PIP10') ~ 'Complex Traits: Annotation',
                                    type %in% c('null_PIP10',
                                                'null_CS') ~ 'Complex Traits: Null',
                                    type %in% c('3tissue_locctrl',
                                                '49tissue_locctrl') ~ 'eQTLs: Location',
                                    type %in% c('3tissue_annotctrl',
                                                '49tissue_annotctrl') ~ 'eQTLs: Annotation')) %>%
  dplyr::mutate(pip = NA) %>%
  distinct(variant, cohort, pip, pip_bin)

bar_plot_pips<- bar_plot_pips %>%
  dplyr::filter(cohort != 'Controls') %>%
  distinct(variant,pip,cohort) %>% 
  group_by(variant, cohort) %>%
  dplyr::filter(pip == max(pip)) %>%
  ungroup() %>%
  dplyr::mutate(pip_bin = cut(pip, pip_bin_breaks, include.lowest = T, ordered_result = T))

bar_plot_pips <- rbind(bar_plot_pips, bar_plot_pips_counts)

bar_plot_pips$cohort <- factor(bar_plot_pips$cohort, levels = c("Controls","Complex Traits", "eQTLs"))
bar_plot_pips<- bar_plot_pips %>% dplyr::mutate(dummy = 'dummy')
bar_plot_pips <- bar_plot_pips %>% dplyr::mutate(cn = paste(cohort, pip_bin, sep = ' '))

bar_plot_pips_counts <- bar_plot_pips %>% group_by(cn) %>% count() %>% dplyr::mutate(dummy = 'dummy')

bar_plot_pips_counts$cn <- factor(bar_plot_pips_counts$cn, levels = c("Controls eQTLs: Location",
                                                                      "Controls eQTLs: Annotation",
                                                                      "Controls Complex Traits: Null",
                                                                      "Controls Complex Traits: Location",
                                                                      "Controls Complex Traits: Annotation",
                                                                      "Complex Traits [0,0.01]",
                                                                      "Complex Traits (0.01,0.1]",
                                                                      "Complex Traits (0.1,0.5]",
                                                                      "Complex Traits (0.5,0.9]",
                                                                      "Complex Traits (0.9,1]",
                                                                      "eQTLs [0,0.01]",
                                                                      "eQTLs (0.01,0.1]",
                                                                      "eQTLs (0.1,0.5]",
                                                                      "eQTLs (0.5,0.9]",
                                                                      "eQTLs (0.9,1]"))

colours <- c('#8B85A7','#B9B5C3','#5F7982','#7FA1AD','#AAC5CF', # controls
             "#DEDBEC","#BEB8D9","#9E95C7","#7E72B4","#5E4FA2", # complex traits 
             "#BEDDE9","#9DCCDD","#7CBBD2","#5CA9C7","#3F96B7") # eQTLs


p <- ggplot(bar_plot_pips_counts, aes(fill = cn, x = dummy, y = n)) +
  geom_bar(position='stack',stat='identity') +
  scale_color_manual(values = colours,aesthetics = c("color","fill")) +
  coord_flip() +
  BuenColors::pretty_plot(fontsize = 20)+
  xlab(NULL)+
  ylab(NULL)+
  labs(fill = "cn")+
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())
p

plt_combined <- p + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot('figures/s1a.pip.structure.pdf',
                   #paste0("../tables/", pdfname, "_tracks.pdf"),
                   plt_combined,
                   base_height = 5,
                   base_width = 10,
                   device = cairo_pdf
)

