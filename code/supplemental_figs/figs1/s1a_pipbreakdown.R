
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
library(VennDiagram)
options(stringsAsFactors = FALSE)

mpra_df <- vroom('data/preprocess/core_mpra.txt.gz') %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0))

# define pip threshold
pip_thresh = 0.5

# fig 1b: trait domains and UKBB/BBJ variants 
# PIP threshold of 0.5
# get domain color info
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

# Complex traits pip breakdown
ukbb <-  mpra_df %>% 
  dplyr::mutate(cohort = case_when(cohort == 'GTEx' ~ 'eQTLs',
                                   cohort %in% c("UKBB", "BBJ") ~ 'Complex Traits',
                                   cohort == 'control' ~ 'Controls',
                                   is.na(cohort) ~ 'Controls')) %>%
  dplyr::filter(cohort == 'Complex Traits', pip > pip_thresh) %>%
  distinct(variant, trait, cohort) %>% 
  left_join(., domain.df %>% dplyr::select(-colors), by='trait') %>%
  dplyr::mutate(domain = ifelse(trait == 'Glaucoma', 'Other',domain))

# count number of domains per variant
ukbb_counts <- ukbb %>% 
  group_by(variant) %>%
  dplyr::mutate(num_domain = n_distinct(domain) ) %>% 
  ungroup()
ukbb_counts <- ukbb_counts %>% 
  dplyr::mutate(domain = ifelse(num_domain > 1, 'More than 1', domain)) %>%
  distinct(variant, domain) %>%
  group_by(domain) %>% 
  count() %>% 
  dplyr::mutate(dummy = 'dummy') %>%
  left_join(., domain.df %>% distinct(domain,colors), by = 'domain') %>%
  dplyr::mutate(colors = ifelse(domain == 'More than 1', 'black', colors))

# prepare for plotting
domaincolors <- ukbb_counts$colors
names(domaincolors) <- ukbb_counts$domain
ukbb_counts$domain <- factor(ukbb_counts$domain, levels = rev(c('More than 1', 
                                                            'Hematopoietic', 'Other','Metabolic',
                                                            'Hepatic',"Renal",'Lipids','Skeletal',
                                                            'Cardiovascular','Immunological',
                                                            'Behavioral','Neurological')))

# plot as bar plot
p1 <- ggplot(ukbb_counts, aes(fill = domain, x = dummy, y = n)) +
  geom_bar(position='stack',stat='identity') +
  scale_color_manual(values = domaincolors,aesthetics = c("color","fill")) +
  coord_flip() +
  BuenColors::pretty_plot(fontsize = 20)+
  xlab(NULL)+
  ylab(NULL)+
  labs(fill = "cn")+
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

p1

# save 
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot(paste0('figures/ukbb.pip.',pip_thresh,'.pdf'),
                   plt_combined,
                   base_height = 3,
                   base_width = 3,
                   device = cairo_pdf
)

# define colors
cohortcolors <-c('#5e4fa2','#a298c4','#C5BFDC')
names(cohortcolors) <-c('Both','UKBB','BBJ')

# separating out UKBB, BBJ, or Both
breakdown <- mpra_df %>% 
  dplyr::filter(cohort %in% c('UKBB','BBJ'), pip > pip_thresh) %>% 
  distinct(variant, cohort, trait) %>% 
  left_join(., domain.df %>% dplyr::select(-colors), by='trait') %>%
  dplyr::mutate(domain = ifelse(trait == 'Glaucoma', 'Other',domain)) %>%
  distinct(variant, cohort, domain) %>% 
  dplyr::group_by(variant) %>% 
  dplyr::mutate(num_domain = n_distinct(domain) ) %>% 
  dplyr::mutate(cohort = ifelse(n_distinct(cohort)>1, 'Both',cohort)) %>%
  ungroup() %>% 
  dplyr::mutate(domain = ifelse(num_domain > 1, 'More than 1', domain)) %>%
  distinct(variant, cohort, domain) %>% 
  dplyr::count(cohort, domain) %>%
  dplyr::mutate(dummy = 'dummy')

# prepare for plotting
breakdown <- breakdown %>% 
  dplyr::mutate(combo = paste(domain, cohort, sep = ' '))
breakdown$combo<- factor(breakdown$combo, levels = c('More than 1 Both', 
                                                     'More than 1 BBJ', 
                                                     'More than 1 UKBB', 
                                                     'Hematopoietic Both', 
                                                     'Hematopoietic BBJ', 
                                                     'Hematopoietic UKBB', 
                                                     'Other Both',
                                                     'Other BBJ',
                                                     'Other UKBB',
                                                     'Metabolic Both',
                                                     'Metabolic BBJ',
                                                     'Metabolic UKBB',
                                                     'Hepatic Both',
                                                     'Hepatic BBJ',
                                                     'Hepatic UKBB',
                                                     "Renal Both",
                                                     "Renal BBJ",
                                                     "Renal UKBB",
                                                     'Lipids Both',
                                                     'Lipids BBJ',
                                                     'Lipids UKBB',
                                                     'Skeletal Both',
                                                     'Skeletal BBJ',
                                                     'Skeletal UKBB',
                                                     'Cardiovascular Both',
                                                     'Cardiovascular BBJ',
                                                     'Cardiovascular UKBB',
                                                     'Immunological Both',
                                                     'Immunological BBJ',
                                                     'Immunological UKBB',
                                                     'Behavioral Both',
                                                     'Behavioral BBJ',
                                                     'Behavioral UKBB',
                                                     'Neurological Both',
                                                     'Neurological BBJ',
                                                     'Neurological UKBB'))
breakdown$domain <- factor(breakdown$domain, levels = rev(c('More than 1', 
                                                        'Hematopoietic', 'Other','Metabolic',
                                                        'Hepatic',"Renal",'Lipids','Skeletal',
                                                        'Cardiovascular','Immunological',
                                                        'Behavioral','Neurological')))
# plot as bar plot
p2<- ggplot(breakdown, aes(fill = cohort, x = dummy, y = n, group = combo)) +
  geom_bar(position='stack',stat='identity') +
  scale_color_manual(values = cohortcolors,aesthetics = c("color","fill")) +
  coord_flip() +
  BuenColors::pretty_plot(fontsize = 20)+
  labs(fill = "Cohort")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
p2

# save 
plt_combined <- p2 + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot(paste0('figures/ukbb.bbj.',pip_thresh,'.pdf'),
                   #paste0("../tables/", pdfname, "_tracks.pdf"),
                   plt_combined,
                   base_height = 3,
                   base_width = 6,
                   device = cairo_pdf
)

# Separate out GTEx tisue systems
gtex <-  mpra_df %>% ungroup() %>% dplyr::mutate(cohort = case_when(cohort == 'GTEx' ~ 'eQTLs',
                                                                 cohort %in% c("UKBB", "BBJ") ~ 'Complex Traits',
                                                                 cohort == 'control' ~ 'Controls')) %>%
  dplyr::filter(cohort == 'eQTLs') %>%
  distinct(variant, pip, cohort, system) %>%
  dplyr::filter(pip > pip_thresh) 

# get counts
gtex_counts <- gtex %>% group_by(variant) %>%
  dplyr::mutate(num_sys = n_distinct(system) ) %>% 
  ungroup() %>%
  dplyr::mutate(system = ifelse(num_sys > 1, 'More than 1', system)) %>%
  distinct(variant, system) %>%
  group_by(system) %>% 
  count() %>% 
  dplyr::mutate(dummy = 'dummy')

# define colors and factors
systemcolors <- c("#86CB66", "#4C72AE", "#F3DE67", "#CE72BB", "#b565b1", "#6BBCCC", "#78a996", "#80584E", "#FD6467", "#80584E", "#FDD380", "black")
names(systemcolors) <- c("Integumentary", "Endocrine", "Circulatory", "Nervous", "Exocrine", "Immune", 
                         "Digestive", "Renal", "Respiratory", "Muscular", "Reproductive", "More than 1")
gtex_counts$system <- factor(gtex_counts$system,  levels = rev(c("More than 1",
                                                             "Nervous",
                                                             "Digestive",
                                                             "Integumentary",
                                                             "Circulatory",
                                                             "Reproductive",
                                                             "Endocrine",
                                                             "Muscular",
                                                             "Respiratory",
                                                             "Immune",
                                                             "Exocrine",
                                                             "Renal"))
                             
)


p3 <- ggplot(gtex_counts, aes(fill = system, x = dummy, y = n)) +
  geom_bar(position='stack',stat='identity') +
  scale_color_manual(values = systemcolors,aesthetics = c("color","fill")) +
  coord_flip() +
  BuenColors::pretty_plot(fontsize = 20)+
  xlab(NULL)+
  ylab(NULL)+
  labs(fill = "System")+
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())
p3

plt_combined <- p3 + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot(paste0('figures/gtex.pip.',pip_thresh,'.pdf'),
                   plt_combined,
                   base_height = 3,
                   base_width = 10,
                   device = cairo_pdf
)

