
# Load in libraries
source("code/utils.R")
library(vroom)
library(tidyverse)
library(BuenColors)
library(ggplot2)
library(dplyr)
library(binom)
library(PNWColors)
options(stringsAsFactors = FALSE)

# read in mpra data
mpra <- vroom('data/preprocess/core_mpra.txt.gz') %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0))

# select pip > 0.5 and active in each cell type
# count domains
tissue.props <- mpra %>%
  dplyr::filter(pip > 0.5) %>% 
  dplyr::mutate(cohort = case_when(cohort %in% c('UKBB','BBJ') ~ 'Complex traits',
                                   cohort == 'GTEx' ~ 'eQTLs',
                                   T ~ cohort)) %>%
  dplyr::filter(cohort == "eQTLs") %>%
  # dplyr::mutate(tissue = case_when(tissue %in% c('Whole_Blood','Liver','Colon_Transverse') ~ tissue,
  #                                  system == 'Nervous' ~ 'Nervous',
  #                                  T ~ NA)) %>%
  dplyr::filter(!is.na(system)) %>%
  distinct(variant, cell_type, emVar, system) %>%
  count(system, cell_type, emVar) %>% 
  group_by(system, cell_type) %>%
  dplyr::mutate(tot = sum(n),
                prop = n/tot,
                lower = binom.confint(n,tot, method="wilson")$lower,
                upper= binom.confint(n, tot, method="wilson")$upper) %>%
  ungroup() %>%
  dplyr::filter(emVar== T)



tissue.props.any <- mpra %>%
  dplyr::filter(pip > 0.5) %>% 
  dplyr::mutate(cohort = case_when(cohort %in% c('UKBB','BBJ') ~ 'Complex traits',
                                   cohort == 'GTEx' ~ 'eQTLs',
                                   T ~ cohort)) %>%
  distinct(cohort, variant, emVar, system) %>%
  group_by(variant, cohort) %>%
  dplyr::mutate(emVar_any = ifelse(any(emVar),T,F)) %>%
  ungroup() %>%
  dplyr::filter(cohort == "eQTLs") %>%
  # dplyr::mutate(tissue = case_when(tissue %in% c('Whole_Blood','Liver','Colon_Transverse') ~ tissue,
  #                                  system == 'Nervous' ~ 'Nervous',
  #                                  T ~ NA)) %>%
  dplyr::filter(!is.na(system)) %>%
  distinct(variant, emVar_any, system) %>%
  count(system, emVar_any) %>% 
  group_by(system) %>%
  dplyr::mutate(tot = sum(n),
                prop = n/tot,
                lower = binom.confint(n,tot, method="wilson")$lower,
                upper= binom.confint(n, tot, method="wilson")$upper) %>%
  ungroup() %>%
  dplyr::filter(emVar_any == T) %>%
  dplyr::mutate(cell_type = 'any')

clrs = c(
  pnw_palette('Bay',8,type='continuous')[5],
  pnw_palette('Bay',8,type='continuous')[3],
  pnw_palette('Starfish')[5],
  pnw_palette('Bay',8,type='continuous')[2],
  pnw_palette('Bay',8,type='continuous')[8],
  "grey",
  'black'
)


names(clrs) = c("A549","HCT116",'HEPG2','K562','SKNSH','any','all')


tissue.props <- rbind(tissue.props,
                      tissue.props.any %>% dplyr::rename(emVar = emVar_any))

tissue.props$cell_type <- factor(tissue.props$cell_type,
                                 levels = c('HCT116','HEPG2','K562','SKNSH','any'))

p3 <- ggplot(tissue.props,aes(x=system, y=prop, group = cell_type, fill = cell_type)) +
  geom_col(position = 'dodge')+
  scale_fill_manual(values=clrs)+
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9), width = 0)+
  pretty_plot() +
  # geom_label(aes(label = paste(n, tot, sep = ' / ')))+ 
  ggtitle('emVars per tissue/system, pip > 0.5')+
  xlab('Cell type')+
  ylab('Proportion of variants')+
  theme(axis.text.x = element_text(angle = 45))

p3

# save
cowplot::save_plot('figures/s8c.pdf',
                   #paste0("../tables/", pdfname, "_tracks.pdf"),
                   p3,
                   base_height = 4,
                   base_width = 6,
                   device = cairo_pdf
)



