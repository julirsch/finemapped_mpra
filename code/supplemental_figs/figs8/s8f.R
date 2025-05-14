unixtools:: set.tempdir("/mnt/sdb/rstudio/sirajl/")
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

mpra_df <- vroom("data/preprocess/core_mpra.txt.gz")  %>% 
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0))

# throw out anything in multiple cohorts for this figure
mpra_df<- mpra_df %>%
  group_by(variant) %>%
  dplyr::filter(n_distinct(cohort) == 1 | all(cohort %in% c('UKBB','BBJ'))) %>% 
  ungroup()

# define pip bin breaks
pip_bin_breaks <- c(0, 0.01, 0.1, 0.5, 0.9, 1.0)

# filter to non-coding test variants
mpra_pre <- mpra_df %>% 
  dplyr::filter(!grepl("ctrl", type), cohort != 'control') %>%
  dplyr::filter(consequence %ni% c("LoF", "missense", "synonymous")) %>%
  group_by(variant, cohort) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  ungroup()  %>%
  distinct(variant, cohort, pip, cell_type, emVar, emVar_any,log2Skew, Skew_logPadj, consequence, consequence2, CRE) %>%
  dplyr::mutate(pip_bin = cut(pip, pip_bin_breaks, include.lowest = T, ordered_result = T))

# repeat for non-coding controls
mpra_pre_ctrls <- mpra_df %>% ungroup()  %>%
  dplyr::filter(cohort == 'control') %>%
  dplyr::filter(consequence %ni% c("LoF", "missense", "synonymous")) %>%
  distinct(variant, cohort, CRE, cell_type, emVar, emVar_any,log2Skew, consequence2, Skew_logPadj, type) %>% 
  dplyr::mutate(cohort = case_when( type %in% c('loc_CS', 'loc_PIP10', 'annot_PIP10','annot_CS', 'null_PIP10','null_CS') ~ 'UKBB',
                                    type %in% c('3tissue_locctrl','49tissue_locctrl','3tissue_annotctrl', '49tissue_annotctrl') ~ 'GTEx')) %>% 
  dplyr::mutate(pip_bin = type)

# separate by cohort
mpra_pre_ukbb <- mpra_pre %>% dplyr::filter(cohort %in% c("UKBB","BBJ"))%>% dplyr::mutate(cohort = 'UKBB')
mpra_pre_gtex <- mpra_pre %>% dplyr::filter(cohort == 'GTEx') 
mpra_pre_ukbb_ctrl <- mpra_pre_ctrls %>% dplyr::filter(type %in% c('loc_CS', 'loc_PIP10', 'annot_PIP10','annot_CS', 'null_PIP10','null_CS')) 
mpra_pre_gtex_ctrl <- mpra_pre_ctrls %>% dplyr::filter(type %in% c('3tissue_locctrl','49tissue_locctrl','3tissue_annotctrl', '49tissue_annotctrl')) 

# filter to variants tested in all 4 cell types
# reduce down to one line per variant, emVars first
filt_func <- function(mpra_pre){
  mpra_pre <- mpra_pre %>% 
    group_by(variant) %>% 
    dplyr::filter(n_distinct(cell_type) == 4)  %>%
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

# filter
mpra_pre_ukbb <- filt_func(mpra_pre_ukbb)
mpra_pre_gtex <- filt_func(mpra_pre_gtex)

mpra_pre_ukbb_ctrl <- filt_func_ctrls(mpra_pre_ukbb_ctrl)
mpra_pre_gtex_ctrl <- filt_func_ctrls(mpra_pre_gtex_ctrl)

mpra_pre_ukbb <- rbind(mpra_pre_ukbb %>% dplyr::select(-c(pip,consequence)), mpra_pre_ukbb_ctrl %>% dplyr::select(-type))
mpra_pre_gtex <- rbind(mpra_pre_gtex %>% dplyr::select(-c(pip,consequence)), mpra_pre_gtex_ctrl %>% dplyr::select(-type))


# get annotations of interest
# and only select highest and lowest pip bins
mpra_pre_ukbb1 <- mpra_pre_ukbb %>%   
  dplyr::mutate(cat = case_when(CRE > 0 ~ 'CRE',
                                consequence2 == 'promoter' ~ 'Promoter',
                                CRE == 0 & consequence2 != 'promoter' ~'None')) %>%
  dplyr::filter(pip_bin %in% c('(0.9,1]', '[0,0.01]') )
mpra_pre_ukbb2 <- mpra_pre_ukbb %>%   
  dplyr::mutate(cat = 'all') %>%
  dplyr::filter(pip_bin %in% c('(0.9,1]', '[0,0.01]') )
mpra_pre_gtex1 <- mpra_pre_gtex %>% 
  dplyr::mutate(cat = case_when(CRE > 0 ~ 'CRE',
                                consequence2 == 'promoter' ~ 'Promoter',
                                CRE == 0 & consequence2 != 'promoter' ~'None')) %>%
  dplyr::filter(pip_bin %in% c('(0.9,1]', '[0,0.01]') )
mpra_pre_gtex2 <- mpra_pre_gtex %>%   
  dplyr::mutate(cat = 'all') %>%
  dplyr::filter(pip_bin %in% c('(0.9,1]', '[0,0.01]') )

# get proportions across pip bins
make_mpra_spec <- function(mpra_pre, typ, attr,name){
  mpra_spec <- mpra_pre %>%
    dplyr::count(pip_bin, cat, !! sym(attr)) %>%
    dplyr::group_by(pip_bin, cat) %>%
    dplyr::mutate(tot = sum(n)) %>%
    dplyr::filter(!! sym(attr) == T) %>%
    dplyr::select(-!! sym(attr)) %>% 
    dplyr::mutate(prop = n / tot,
                  lower = binom.confint(n, tot, method='wilson')$lower,
                  upper = binom.confint(n, tot, method='wilson')$upper) %>%
    arrange(-prop) %>% 
    ungroup()
  print(mpra_spec)
  return(mpra_spec)
}


mpra_spec_ukbb1 <- make_mpra_spec(mpra_pre_ukbb1, "UKBB",'emVar','emVar')
mpra_spec_ukbb2 <- make_mpra_spec(mpra_pre_ukbb2, "UKBB",'emVar','emVar')

mpra_spec_gtex1 <- make_mpra_spec(mpra_pre_gtex1, "GTEx",'emVar','emVar')
mpra_spec_gtex2 <- make_mpra_spec(mpra_pre_gtex2, "GTEx",'emVar','emVar')


# combine 
mpra_spec_ukbb <- rbind(mpra_spec_ukbb1 %>% dplyr::filter(cat != 'None'), mpra_spec_ukbb2)
mpra_spec_gtex <- rbind(mpra_spec_gtex1 %>% dplyr::filter(cat != 'None'), mpra_spec_gtex2)

mpra_spec <- rbind(mpra_spec_ukbb %>% dplyr::mutate(cohort = 'Complex Traits'),
                       mpra_spec_gtex %>% dplyr::mutate(cohort = 'eQTLs')) 


# get enrichment of highest PIP bin over lowest PIP bin
mpra_enrich <- mpra_spec %>% 
  dplyr::mutate(se= sqrt((prop)*(1-prop)/n)) %>%
  dplyr::mutate(se = (se/prop)^2) %>%
  dplyr::select(-lower,-upper) %>%
  group_by(cat, cohort) %>%
  summarise(prop = prop[pip_bin == "(0.9,1]"] / prop[pip_bin == "[0,0.01]"],
            se = sqrt(se[pip_bin == "(0.9,1]"] + se [pip_bin == "[0,0.01]"])) %>%
  ungroup() %>%
  dplyr::mutate(se = se * prop) %>%
  dplyr::mutate(upper = prop + se*1.96,
                lower = prop - se*1.96) %>%
  dplyr::mutate(lower = pmax(lower,0))


# plot
p <- mpra_enrich %>%
  dplyr::filter(cat != 'Promoter') %>%
  ggplot(., aes(x = cat, y = prop, group=cohort, fill=cohort)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                size = 0.5, width = 0, col = "black",
                position = position_dodge(0.9, preserve = "single")) +
  pretty_plot() +
  theme(#legend.title = element_blank(),
    axis.text.x=element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"))+
  # legend.text = element_text(size = 6),
  # legend.position = "right",
  # legend.justification = c(1, 1),
  #legend.title = element_text(margin = margin(0, 0, 0, 0)),
  #legend.background = element_blank()) +
  scale_fill_manual(values = jdb_palette('brewer_spectra')) +
  #guides(colour = guide_legend(show = FALSE))+
  ylab("Proportion of variants") +
  xlab("")+
  theme(legend.position = "none")
p

# save
plt_combined <- p
plt_combined
cowplot::save_plot('figures/s8f.pdf',
                   #paste0("../tables/", pdfname, "_tracks.pdf"),
                   plt_combined,
                   base_height = 3,
                   base_width = 3,
                   device = cairo_pdf
)





