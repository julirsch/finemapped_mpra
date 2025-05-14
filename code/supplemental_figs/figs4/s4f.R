
source("code/utils.R")
library(vroom)
library(tidyverse)
library(BuenColors)
library(ggplot2)
library(dplyr)
library(binom)
library(tidyr)
options(stringsAsFactors = FALSE)

# read in sei annotations
# only keep silencer-related annotations
sei <- vroom('/mnt/sdb/gwas_eqtl_mpra/data/preprocess/sei_mpra.txt.gz') 

sei <- sei %>%
  dplyr::select('variant', contains('Poly'), contains("Het"))

# binarize
sei[,2:6] <- ifelse(sei[,2:6] > 0, 1, 0)


# read in processed MPRA data
# remove coding variants
# join in sei data
mpra_df <- vroom('data/preprocess/core_mpra.txt.gz') %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0)) %>%
  group_by(variant) %>%
  dplyr::filter(consequence %ni% c("LoF", "missense", "synonymous"), promoter == 0) %>%
  ungroup() %>% 
  distinct(variant, cohort, trait, region, gene, tissue, pip, cell_type, active, active_any, type, emVar, emVar_any) %>% 
  left_join(sei, by = 'variant')

# remove test variants that overlap cohorts
# define emVar and active
mpra_df <- mpra_df %>%
  dplyr::filter(cohort != 'controls') %>%
  dplyr::mutate(cohort = ifelse(cohort %in% c('UKBB','BBJ'), 'Complex Traits', cohort)) %>%
  group_by(variant) %>%
  dplyr::filter(all(cohort == 'GTEx') | all(cohort == 'Complex Traits')) %>% 
  dplyr::mutate(active_any = ifelse(any(active ==T), T, F)) %>%
  dplyr::mutate(emVar_any = ifelse(any(emVar ==T), T, F)) %>%
  ungroup() %>%
  dplyr::select(variant, cohort, active_any, emVar_any, starts_with('SEI')) %>%
  distinct()

# get counts in each cohort
counts <- mpra_df %>%
  group_by(cohort, active_any) %>%
  dplyr::mutate(tot = n()) %>%
  distinct(cohort, active_any, tot)

# get counts per SEI category
to_plot <- mpra_df %>% 
  group_by(cohort, active_any) %>%
  dplyr::select(-variant,-emVar_any) %>%
  summarise_all(sum) %>%
  ungroup() %>%
  left_join(counts, by = c('cohort','active_any')) %>%
  pivot_longer(cols = starts_with('SEI'),
               names_to = 'SEI',
               names_prefix = 'SEI.',
               values_to = 'num') %>%
  dplyr::mutate(prop = num/tot,
                lower = binom.confint(num,tot, method="wilson")$lower,
                upper = binom.confint(num,tot, method="wilson")$upper) %>%
  dplyr::mutate(cohort = ifelse(cohort == 'GTEx', 'eQTLs', cohort))  

# plot eQTLs
pg <- to_plot %>% 
  #dplyr::filter(SEI %ni% c('Low.signal','Weak.Polycomb')) %>%
  dplyr::filter(SEI %in% c("Polycomb","Heterochromatin")) %>%
  dplyr::filter(cohort == 'eQTLs') %>% 
  dplyr::rename('active' = 'active_any') %>%
  ggplot(.,aes(x=SEI, y= prop, group = active, fill = active)) +
  geom_col(position = 'dodge')+
  scale_fill_manual(values = c('orange','dark green'))+
  geom_errorbar(aes(ymin = lower, ymax =upper), position = position_dodge(0.9), width = 0)+
  pretty_plot()+
  ggtitle('eQTL: Silencing annotations')+
  ylab('Proportion of variants with annotation')+
  xlab('Annotation')+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(0,0.06)

pg

# plot complex traits
pu <- to_plot %>% 
  #dplyr::filter(SEI %ni% c('Low.signal','Weak.Polycomb')) %>%
  dplyr::filter(SEI %in% c("Polycomb","Heterochromatin")) %>%
  dplyr::filter(cohort == 'Complex Traits') %>% 
  dplyr::rename('active' = 'active_any') %>%
  ggplot(.,aes(x=SEI, y= prop, group = active, fill = active)) +
  geom_col(position = 'dodge')+
  scale_fill_manual(values = c('orange','dark green'))+
  geom_errorbar(aes(ymin = lower, ymax =upper), position = position_dodge(0.9), width = 0)+
  pretty_plot()+
  ggtitle('Complex Traits: Silencing annotations')+
  ylab('Proportion of variants with annotation')+
  xlab('Annotation')+scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(0,0.06)

pu
# save 
plt_combined <- ggarrange(pg, pu, 
                          ncol = 1, nrow = 2)

plt_combined

# save off files and figures

cowplot::save_plot(paste0('figures/s4f.pdf'),
                   plt_combined,
                   base_height = 3,
                   base_width = 6,
                   device = cairo_pdf
)

