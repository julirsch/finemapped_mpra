

# Load libraries
source("code/utils.R")
library(tidyverse)
library(BuenColors)
library(binom)
library(ggplot2)
library(cowplot)
library(patchwork)
library("PNWColors")
library(vroom)
options(stringsAsFactors = FALSE)

mpra_df <- vroom("/mnt/sdb/gwas_eqtl_mpra/data/preprocess/core_mpra.txt.gz") 

# throw out anything in multiple cohorts for this figure
mpra_df<- mpra_df %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0)) %>% 
  group_by(variant) %>%
  dplyr::filter(n_distinct(cohort) == 1 | all(cohort %in% c('UKBB','BBJ'))) %>% 
  ungroup()


pip_bin_breaks <- seq(0,1,0.01)

# filter to test, non-coding variants in CREs, 
# sort into pip bins
mpra_pre <- mpra_df %>% 
  dplyr::filter(!grepl("ctrl", type), cohort != 'control') %>%
  dplyr::filter(consequence %ni% c("LoF", "missense", "synonymous")) %>%
  dplyr::filter(CRE>0) %>%
  group_by(variant, cohort) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  ungroup()  %>%
  distinct(variant, cohort, pip, cell_type, emVar, emVar_any,log2Skew, Skew_logPadj, consequence) %>%
  dplyr::mutate(pip_bin = cut(pip, pip_bin_breaks, include.lowest = T, ordered_result = T))

# make same DF but for controls
mpra_pre_ctrls <- mpra_df %>% ungroup()  %>%
  dplyr::filter(cohort == 'control') %>%
  dplyr::filter(CRE > 0) %>%
  #### throw out coding variants
  dplyr::filter(consequence %ni% c("LoF", "missense", "synonymous")) %>%
  distinct(variant, cohort, cell_type, emVar, emVar_any,log2Skew, Skew_logPadj, type) %>% 
  dplyr::mutate(cohort = case_when( type %in% c('loc_CS', 'loc_PIP10', 'annot_PIP10','annot_CS', 'null_PIP10','null_CS') ~ 'UKBB',
                                    type %in% c('3tissue_locctrl','49tissue_locctrl','3tissue_annotctrl', '49tissue_annotctrl') ~ 'GTEx')) %>% 
  dplyr::mutate(pip_bin = type)


mpra_pre_ukbb <- mpra_pre %>% dplyr::filter(cohort %in% c("UKBB","BBJ"))%>% dplyr::mutate(cohort = 'UKBB')
mpra_pre_gtex <- mpra_pre %>% dplyr::filter(cohort == 'GTEx') 
mpra_pre_ukbb_ctrl <- mpra_pre_ctrls %>% dplyr::filter(type %in% c('loc_CS', 'loc_PIP10', 'annot_PIP10','annot_CS', 'null_PIP10','null_CS')) 
mpra_pre_gtex_ctrl <- mpra_pre_ctrls %>% dplyr::filter(type %in% c('3tissue_locctrl','49tissue_locctrl','3tissue_annotctrl', '49tissue_annotctrl')) 

# filter to variants that sruvived testing in all 4 cell types
filt_func <- function(mpra_pre){
  mpra_pre <- mpra_pre %>% 
    group_by(variant) %>% 
    dplyr::filter(n_distinct(cell_type) == 4)  %>%
    # take emVars first
    arrange(-emVar,-Skew_logPadj) %>%  
    filter(row_number() == 1) %>%
    ungroup() 
  return(mpra_pre)
}

filt_func_ctrls <- function(mpra_pre){
  mpra_pre <- mpra_pre %>% 
    group_by(variant, type) %>% 
    dplyr::filter(n_distinct(cell_type) == 4)  %>%
    # take emVars first
    arrange(-emVar,-Skew_logPadj) %>%  
    filter(row_number() == 1) %>%
    ungroup() 
  return(mpra_pre)
}


mpra_pre_ukbb <- filt_func(mpra_pre_ukbb)
mpra_pre_gtex <- filt_func(mpra_pre_gtex)

mpra_pre_ukbb_ctrl <- filt_func_ctrls(mpra_pre_ukbb_ctrl)
mpra_pre_gtex_ctrl <- filt_func_ctrls(mpra_pre_gtex_ctrl)

# combine
mpra_pre_ukbb <- rbind(mpra_pre_ukbb %>% dplyr::select(-c(pip,consequence)), mpra_pre_ukbb_ctrl %>% dplyr::select(-type))
mpra_pre_gtex <- rbind(mpra_pre_gtex %>% dplyr::select(-c(pip,consequence)), mpra_pre_gtex_ctrl %>% dplyr::select(-type))


make_mpra_spec <- function(mpra_pre, typ, attr,name){
  mpra_spec <- mpra_pre %>%
    dplyr::count(pip_bin, !! sym(attr)) %>%
    dplyr::group_by(pip_bin) %>%
    dplyr::mutate(tot = sum(n)) %>%
    dplyr::mutate(prop = n / tot,
                  lower = binom.confint(n, tot, method='wilson')$lower,
                  upper = binom.confint(n, tot, method='wilson')$upper) %>%
    ungroup() %>% 
    dplyr::filter(!! sym(attr) == T) %>%
    dplyr::select(-!! sym(attr)) %>% 
    arrange(-prop) 
  print(mpra_spec)
  return(mpra_spec)
}


# get emVar proportion across pip bins
mpra_spec_ukbb <- make_mpra_spec(mpra_pre_ukbb, "UKBB",'emVar','emVar')
mpra_spec_ukbb <- mpra_spec_ukbb %>% 
  dplyr::mutate(pip_bin = as.character(pip_bin)) %>% 
  dplyr::mutate(pip_bin = ifelse(pip_bin == 'annot_PIP10','annot-matched',pip_bin)) %>%
  dplyr::mutate(pip_bin = ifelse(pip_bin == 'loc_PIP10','loc-matched',pip_bin)) %>%
  dplyr::mutate(pip_bin = ifelse(pip_bin == 'null_PIP10','Null',pip_bin)) %>%
  dplyr::filter(pip_bin %ni% c('null_CS', 'annot_CS','loc_CS'))


mpra_spec_gtex <- make_mpra_spec(mpra_pre_gtex, "GTEx",'emVar','emVar')
mpra_spec_gtex <- mpra_spec_gtex %>% 
  dplyr::mutate(pip_bin = as.character(pip_bin)) %>% 
  dplyr::mutate(pip_bin = ifelse(pip_bin == '49tissue_annotctrl','annot-matched',pip_bin)) %>%
  dplyr::mutate(pip_bin = ifelse(pip_bin == '49tissue_locctrl','loc-matched',pip_bin)) %>%
  dplyr::filter(pip_bin %ni% c('3tissue_annotctrl', '3tissue_locctrl'))


# keep only coding ones in the end
mpra_spec <- bind_rows(mpra_spec_ukbb %>% dplyr::mutate(cohort = 'Complex Traits'),
                       mpra_spec_gtex %>% dplyr::mutate(cohort = 'eQTLs')) %>%
  dplyr::filter(pip_bin %ni% c("Null", 'annot-matched','loc-matched'))

mpra_spec

# plot
p <- mpra_spec  %>% 
  group_by(pip_bin) %>%
  mutate(Group = cur_group_id()/100+0.01) %>% 
  ggplot(., aes(x=Group, y=prop, group=cohort, color = cohort))+
  geom_smooth(method='lm', se=F, linetype = 'dashed')+
  geom_smooth(method='loess', se=F, span = 1)+
  pretty_plot()+
  scale_x_continuous(trans='log2', breaks = c(0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1))+
  scale_color_manual(values = jdb_palette('brewer_spectra'))+
  xlab('log2(PIP)')+
  ylab('Proportion of emVars')



p


p

cowplot::save_plot('figures/s8e.pdf',
                   #paste0("../tables/", pdfname, "_tracks.pdf"),
                   p,
                   base_height = 3,
                   base_width = 5,
                   device = cairo_pdf
)

