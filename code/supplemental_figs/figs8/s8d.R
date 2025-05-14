
# Load in libraries
source("code/utils.R")
library(vroom)
library(tidyverse)
library(BuenColors)
library(ggplot2)
library(dplyr)
library(binom)
library(tidyr)
library(scales)
library(ggpubr)
library(grid)
library(PNWColors)
options(stringsAsFactors = FALSE)

# read in mpra data
mpra <- vroom('/mnt/sdb/gwas_eqtl_mpra/data/preprocess/core_mpra.txt.gz') %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0)) %>%
  dplyr::filter(CRE >0)


# Get variants that are emVars in only a single cell type
mpra_ct_emvar <- mpra %>% 
  dplyr::filter(cohort == 'GTEx') %>% 
  distinct(variant, cell_type, emVar) %>%
  dplyr::group_by(variant,cell_type) %>%
  dplyr::mutate(emVar = ifelse(any(emVar ==T),T,F)) %>%
  ungroup() %>%
  distinct(variant, cell_type, emVar) %>%
  pivot_wider(names_from = cell_type,
              values_from = emVar) %>%
  dplyr::filter((K562+SKNSH+HCT116+HEPG2) == 1) 


# for each cell type x each tissue system, get OR in matched cell type vs others
mpra.tissues <- mpra %>% 
  dplyr::filter(cohort == 'GTEx') %>% 
  distinct(variant, pip,tissue) %>%
  dplyr::filter(tissue %in% c('Colon_Transverse',
                              'Whole_Blood',
                              'Liver')) %>%
  dplyr::filter(pip > 0.1) %>%
  distinct(variant, tissue) %>%
  dplyr::mutate(val = 1) %>%
  pivot_wider(names_from = tissue,
              values_from = val) %>%
  replace(is.na(.),0)

mpra.sys <- mpra %>% 
  dplyr::filter(cohort == 'GTEx') %>% 
  distinct(variant, pip,system) %>%
  dplyr::filter(system == 'Nervous') %>%
  dplyr::filter(pip > 0.1) %>%
  distinct(variant, system)

mpra.annot <- left_join(mpra_ct_emvar, mpra.tissues) %>%
  replace(is.na(.),0) %>%
  dplyr::mutate(Nervous = ifelse(variant %in% unique(mpra.sys$variant), 1,0)) 

prop.table <- mpra.annot %>%
  pivot_longer(cols = c('K562','HCT116','HEPG2','SKNSH'),
               names_to = 'cell_type',
               values_to = 'emVar')  %>%
  pivot_longer(cols = c('Colon_Transverse','Whole_Blood','Liver','Nervous'),
               names_to = 'Tissue',
               values_to = 'in.tissue') %>% 
  count(cell_type, emVar, Tissue, in.tissue) %>%
  arrange(Tissue, cell_type, emVar, in.tissue) %>%
  group_by(cell_type, Tissue) %>%
  dplyr::summarise(p.val = prop.test(x=matrix(c(n[in.tissue==1 & emVar == T],
                                                n[in.tissue==1 & emVar == F],
                                                n[in.tissue==0 & emVar == T],
                                                n[in.tissue==0 & emVar == F]),
                                              nrow=2),
                                     alternative = 'g')$p.value,
                   p_x = prop.test(x=matrix(c(n[in.tissue==1 & emVar == T],
                                              n[in.tissue==1 & emVar == F],
                                              n[in.tissue==0 & emVar == T],
                                              n[in.tissue==0 & emVar == F]),
                                            nrow=2),
                                   alternative = 'g')$estimate[[1]],
                   p_y = prop.test(x=matrix(c(n[in.tissue==1 & emVar == T],
                                              n[in.tissue==1 & emVar == F],
                                              n[in.tissue==0 & emVar == T],
                                              n[in.tissue==0 & emVar == F]),
                                            nrow=2),
                                   alternative = 'g')$estimate[[2]],
                   n_x = sum(n[emVar==T]),
                   n_y = sum(n[emVar==F]),
                   se = sqrt(p_x * (1 - p_x) / n_x + p_y * (1 - p_y) / n_y),
                   z_crit =  qnorm(p = 0.05 / 2, mean = 0, sd = 1, lower.tail = FALSE),
                   d_hat = p_x - p_y,
                   lcl = d_hat - z_crit * se,
                   ucl = d_hat + z_crit * se)





# define colors
clrs = c(
  #pnw_palette('Bay',8,type='continuous')[5],
  pnw_palette('Bay',8,type='continuous')[3],
  pnw_palette('Starfish')[5],
  pnw_palette('Bay',8,type='continuous')[2],
  pnw_palette('Bay',8,type='continuous')[8]
)


names(clrs) = c("Colon_Transverse", "Liver", "Whole_Blood","Nervous")

prop.table <- prop.table %>%
  dplyr::mutate(Tissue = case_when(Tissue == 'Whole_Blood' ~ 'Whole Blood',
                                   Tissue == 'Colon_Transverse' ~ 'Transverse Colon',
                                   T ~ Tissue)) 

prop.table$Tissue <- factor(prop.table$Tissue, levels = rev(c('Transverse Colon',
                                                              'Liver',
                                                              'Whole Blood',
                                                              'Nervous')))

p1 <- ggplot(prop.table,aes(x=cell_type, y=Tissue))+
  #geom_tile(col="black", fill="white") + 
  geom_tile(col="black",fill = NA,  na.value = NA, linewidth = 0.5) + 
  geom_point(aes(size = -log10(p.val), color=d_hat), shape=15) +
  #geom_text(aes(label = round(d_hat,2)))+
  theme_classic() +
  scale_colour_gradient2(low = muted('blue'),
                         mid = 'white',
                         high = muted('red'),
                         midpoint = 0,
                         breaks=c(-0.06,-0.03,0,0.03,0.06),
                         name = 'Difference in proportion')+
  #scale_color_gradientn(colors = jdb_palette("brewer_celsius")) +
  #scale_color_continuous(type = 'viridis', limits = c(0, 0.02), oob = scales::squish)+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(text=element_text(family="Helvetica"),
        axis.title.x=element_blank(),
        panel.background = element_blank(),
        # axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title.y=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  scale_size(
    limits = c(0,10),
    name = "-log10(p-value)",
    range = c(1, 10),
    breaks= c(2.5,5,7.5),
    guide = "legend",
  )+
  coord_equal()

p1



cowplot::save_plot('figures/s8d.pdf',
                   #paste0("../tables/", pdfname, "_tracks.pdf"),
                   plt_combined,
                   base_height = 4,
                   base_width = 4,
                   
                   device = cairo_pdf
)

