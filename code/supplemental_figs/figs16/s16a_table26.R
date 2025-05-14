# Load in libraries

source("code/utils.R")
library(tidyverse)
library(BuenColors)
library(reticulate)
library(binom)
library(ggplot2)
library(cowplot)
library(patchwork)
library(vroom)
library(GenomicRanges)
library(ComplexHeatmap)
library(ggpubr)
library(dplyr)
options(stringsAsFactors = FALSE)

# filtering function to only keep variants tested in 4 cell types
filt_func <- function(mpra_pre){
  mpra_pre <- mpra_pre %>%
    group_by(variant) %>%
    dplyr::filter(n_distinct(cell_type) == 4)  %>%
    # take emVars first
    arrange(-active,-logPadj_BF) %>%
    filter(row_number() == 1) %>%
    ungroup()
  return(mpra_pre)
}
# filtering function for controls
filt_func_ctrls <- function(mpra_pre){
  mpra_pre <- mpra_pre %>%
    group_by(variant, type) %>%
    dplyr::filter(n_distinct(cell_type) == 4)  %>%
    # take emVars first
    arrange(-active,-logPadj_BF) %>%
    filter(row_number() == 1) %>%
    ungroup()
  return(mpra_pre)
}

# read in processed MPRA data
mpra <- vroom('data/preprocess/core_mpra.txt.gz') %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0)) 


# get annotations
#DHS_Meuleman.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/meuleman/meuleman_en_gr.rds")
atac <- readRDS("/mnt/sdb/ukbb-finemapping/data/annotations/meuleman/meuleman_cpm_gr.rds")


# Annotate GTEx variants
variant.gr <- mpra %>%
  distinct(variant, chromosome, position) %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "position", 
                           end.field = "position",
                           keep.extra.columns = T)



# metadata from https://www.meuleman.org/DHS_Index_and_Vocabulary_metadata.html
metadata <- vroom::vroom("/mnt/sdb/gwas_eqtl_mpra/reviews/data/meuleman_dhs_metadata.txt") %>%
  dplyr::rename('biosample' = 'Biosample name',
                'altius' = "Altius Aggregation ID",
                'sample.type'= "Biosample type",
                'state'="Biological state")

# annotate each peak with euclidean norm, Meuleman component
test <- as.data.frame(atac) %>%
  dplyr::select(-c(seqnames,start,end,width,strand)) 
euclid_norm <- function(x){sqrt(sum(x**2))}
norms <- as.data.frame(apply(test[,1:dim(test)[2]],1,euclid_norm)) 
norms <- norms %>% dplyr::rename('enorm' = 'apply(test[, 1:dim(test)[2]], 1, euclid_norm)')
norms <- tibble::rownames_to_column(norms, "peak")
# import full DHS index
index <- vroom::vroom("/mnt/sdb/gwas_eqtl_mpra/reviews/data/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz")
norms <- norms %>% left_join(index %>% dplyr::select(identifier, component, numsamples),
                             by = c('peak' = 'identifier'))
norms$numsamples <- as.numeric(norms$numsamples)
norms <-  norms %>% dplyr::mutate(color.codes = case_when(component == 'Placental / trophoblast' ~ '#fadc13',
                                                          component == 'Lymphoid' ~ '#ffa85a',
                                                          component == 'Myeloid / erythroid' ~'#de2e1e',
                                                          component == 'Cardiac' ~ '#33b51b',
                                                          component == 'Musculoskeletal' ~ '#77993d',
                                                          component == 'Vascular / endothelial' ~ '#5a8727',
                                                          component == 'Primitive / embryonic' ~ '#6bcde1',
                                                          component == 'Neural' ~ '#4c6db8',
                                                          component == 'Digestive' ~ '#28988d',
                                                          component == 'Stromal A' ~ '#ae49bb',
                                                          component == 'Stromal B' ~ '#7e4bb7',
                                                          component == 'Organ devel. / renal' ~ '#743217',
                                                          component == 'Cancer / epithelial' ~ '#222957',
                                                          component == 'Pulmonary devel.' ~ '#be4b18',
                                                          component == 'Renal / cancer' ~ '#516a7a',
                                                          component == 'Tissue invariant' ~ '#ccc8c8'))


# get the MPRA variants selected
# throw out anything in multiple cohorts for this analysis
mpra<- mpra %>%
  group_by(variant) %>%
  dplyr::filter(n_distinct(cohort) == 1 | all(cohort %in% c('UKBB','BBJ'))) %>% 
  ungroup()

pip_bin_breaks <- c(0, 0.01, 0.1, 0.5, 0.9, 1.0)

# sort non-coding test variants in CREs into PIP bins
mpra_pre <- mpra %>% 
  dplyr::filter(CRE > 0) %>%
  dplyr::filter(cohort != 'control') %>%
  dplyr::filter(consequence %ni% c("LoF", "missense", "synonymous")) %>%
  group_by(variant, cohort) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  ungroup()  %>%
  distinct(variant, chromosome, position, cohort, cell_type, pip, active, active_any, logPadj_BF, trait, region, tissue, gene, system) %>%
  dplyr::mutate(pip_bin = cut(pip, pip_bin_breaks, include.lowest = T, ordered_result = T))

# same for controls
mpra_pre_ctrls <- mpra %>% ungroup()  %>%
  dplyr::filter(CRE >0) %>%
  dplyr::filter(cohort == 'control') %>%
  dplyr::filter(consequence %ni% c("LoF", "missense", "synonymous")) %>%
  distinct(variant, chromosome, position, cohort, cell_type, active, active_any, logPadj_BF, trait, region, tissue, gene, system, type) %>% 
  dplyr::mutate(cohort = case_when( type %in% c('loc_CS', 'loc_PIP10', 'annot_PIP10','annot_CS', 'null_PIP10','null_CS') ~ 'Complex Traits',
                                    type %in% c('3tissue_locctrl','49tissue_locctrl','3tissue_annotctrl', '49tissue_annotctrl') ~ 'eQTLs')) %>% 
  #### make pip bins based on the type
  dplyr::mutate(pip_bin = type)

mpra_pre_ukbb <- mpra_pre %>% dplyr::filter(cohort %in% c("UKBB","BBJ"))%>% dplyr::mutate(cohort = 'Complex Traits') 
mpra_pre_gtex <- mpra_pre %>% dplyr::filter(cohort == 'GTEx') %>% dplyr::mutate(cohort = 'eQTLs')
mpra_pre_ukbb_ctrl <- mpra_pre_ctrls %>% dplyr::filter(type %in% c('loc_CS', 'loc_PIP10', 'annot_PIP10','annot_CS', 'null_PIP10','null_CS')) 
mpra_pre_gtex_ctrl <- mpra_pre_ctrls %>% dplyr::filter(type %in% c('3tissue_locctrl','49tissue_locctrl','3tissue_annotctrl', '49tissue_annotctrl')) 

# filter down to just variants tested in all 4 cell types
mpra_pre_ukbb <- filt_func(mpra_pre_ukbb)
mpra_pre_gtex <- filt_func(mpra_pre_gtex)
mpra_pre_ukbb_ctrl <- filt_func_ctrls(mpra_pre_ukbb_ctrl)
mpra_pre_gtex_ctrl <- filt_func_ctrls(mpra_pre_gtex_ctrl)


mpra_df <- rbind(mpra_pre_ukbb %>% dplyr::select(-pip),
                 mpra_pre_gtex %>% dplyr::select(-pip), 
                 mpra_pre_gtex_ctrl %>% dplyr::select(-type),
                 mpra_pre_ukbb_ctrl %>% dplyr::select(-type)) 



get_dhs <- function(b,c, a, in_df){
  variant.gr <- in_df %>% 
    dplyr::filter(pip_bin == b, cohort == c, active_any == a) %>%
    distinct(variant, chromosome, position) %>%
    dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                  position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
    makeGRangesFromDataFrame(., 
                             seqnames.field = "chromosome", 
                             start.field =  "position", 
                             end.field = "position",
                             keep.extra.columns = T)
  
  # make metadata
  df <- as.data.frame(countOverlaps(atac,variant.gr)) %>%
    dplyr::rename("overlap" = "countOverlaps(atac, variant.gr)")
  df <- tibble::rownames_to_column(df, "identifier")
  df <- df %>% left_join(norms, by=c('identifier'='peak'))
  df <- df %>% 
    group_by(component) %>% 
    dplyr::summarise('num_in_group' = sum(overlap)) %>%
    ungroup() %>%
    dplyr::mutate(tot = sum(num_in_group)) %>%
    dplyr::mutate(prop = num_in_group/tot) %>%
    dplyr::mutate(bin = b,
                  cohort = c,
                  active_any = a)
  return(df)
}


mpra_df <- mpra_df %>% 
  dplyr::mutate(pip_bin = as.character(pip_bin)) %>% 
  dplyr::mutate(pip_bin = ifelse(pip_bin == 'annot_PIP10','annot-matched',pip_bin)) %>%
  dplyr::mutate(pip_bin = ifelse(pip_bin == 'loc_PIP10','loc-matched',pip_bin)) %>%
  dplyr::mutate(pip_bin = ifelse(pip_bin == 'null_PIP10','Null',pip_bin)) %>%
  dplyr::filter(pip_bin %ni% c('null_CS', 'annot_CS','loc_CS')) %>%
  dplyr::mutate(pip_bin = ifelse(pip_bin == '49tissue_annotctrl','annot-matched',pip_bin)) %>%
  dplyr::mutate(pip_bin = ifelse(pip_bin == '49tissue_locctrl','loc-matched',pip_bin)) %>%
  dplyr::filter(pip_bin %ni% c('3tissue_annotctrl', '3tissue_locctrl'))

mpra_out <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(mpra_out) <- c("component",'num_in_group','tot','prop','bin','cohort','active_any')
for (row in 1:nrow(mpra_df %>% distinct(pip_bin, cohort, active_any))){
  print(row)
  tmp <- mpra_df %>% ungroup() %>%droplevels() %>% distinct(pip_bin, cohort, active_any)
  b <- as.character(tmp[row,'pip_bin'][[1]])
  c <- tmp[row,'cohort'][[1]]
  a <- tmp[row, 'active_any'][[1]]
  mpra_out <- rbind(mpra_out, get_dhs(b,c,a,mpra_df))
}

mpra_out <- mpra_out %>%
  dplyr::rename(tot_in_dhs = tot) %>%
  left_join(.,
            mpra_df %>% 
              dplyr::distinct(pip_bin, cohort, active_any, variant,chromosome,position) %>% 
              dplyr::group_by(pip_bin,cohort,active_any) %>%
              dplyr::summarise(tot = n()) %>%
              ungroup(),
            by=c('bin' = 'pip_bin', 'cohort','active_any'))

#mpra_out %>% vroom::vroom_write(paste0('/mnt/sdb/gwas_eqtl_mpra/reviews/tables/',Sys.Date(),'.dhs.meuleman.components.table.bin.cohort.emvar.txt'))

#prepare for plotting
tmp <- mpra_out %>%
  group_by(cohort, bin, active_any) %>% 
  dplyr::mutate(lower = binom.confint(num_in_group,tot_in_dhs, method="wilson")$lower,
                upper = binom.confint(num_in_group,tot_in_dhs, method="wilson")$upper)%>%
  distinct(cohort, bin, active_any, num_in_group,tot_in_dhs, prop, lower, upper, component) %>%
  dplyr::mutate(name = paste0(cohort, ': active? ', active_any))%>% 
  dplyr::filter(bin == '(0.9,1]') %>%
  dplyr::mutate(name = case_when((active_any == T & cohort == 'Complex Traits') ~ 'Complex Traits active',
                                 (active_any == F & cohort == 'Complex Traits') ~ 'Complex Traits inactive',
                                 (active_any == T & cohort == 'eQTLs') ~ 'eQTLs active',
                                 (active_any == F & cohort == 'eQTLs') ~ 'eQTLs inactive'))

tmp$name <- factor(tmp$name, levels = c('eQTLs inactive', 'eQTLs active', 'Complex Traits inactive', 'Complex Traits active'))

# Plot
p <- tmp %>%
  ungroup() %>% 
  dplyr::group_by(component, active_any) %>%
  #dplyr::mutate(num_in_group = sum(num_in_group),
  #                 tot_in_dhs = sum(tot_in_dhs)) %>%
  #ungroup() %>% 
  pivot_wider(id_cols = c(cohort, component), names_from = c(active_any), values_from = c(num_in_group, tot_in_dhs)) %>%
  dplyr::group_by(cohort, component) %>% 
  dplyr::summarize(diff = PropCIs::wald2ci(num_in_group_TRUE, tot_in_dhs_TRUE, num_in_group_FALSE, tot_in_dhs_FALSE, conf.level = 0.95, adjust = "AC")$estimate,
                   diff_lower = PropCIs::wald2ci(num_in_group_TRUE, tot_in_dhs_TRUE, num_in_group_FALSE, tot_in_dhs_FALSE, conf.level = 0.95, adjust = "AC")$conf.int[1],
                   diff_upper = PropCIs::wald2ci(num_in_group_TRUE, tot_in_dhs_TRUE, num_in_group_FALSE, tot_in_dhs_FALSE, conf.level = 0.95, adjust = "AC")$conf.int[2]) %>% 
  ggplot(., aes(y = reorder(component, diff), x = diff, group = cohort, color = component, shape = cohort)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(position = position_dodge(width = .75)) +
  geom_errorbarh(aes(xmin = diff_lower, xmax = diff_upper, height = 0), position = position_dodge(width = .75)) +
  scale_color_manual(breaks = c("Musculoskeletal",
                                "Tissue invariant",
                                "Stromal A",
                                "Stromal B",
                                "Primitive / embryonic",
                                "Renal / cancer",
                                "Lymphoid",
                                "Placental / trophoblast",
                                "Organ devel. / renal",
                                "Vascular / endothelial",
                                "Pulmonary devel.",
                                "Cardiac",
                                "Myeloid / erythroid",
                                "Digestive",
                                "Neural",
                                "Cancer / epithelial"),
                     values = c( "#77993d",
                                 "#ccc8c8",
                                 "#ae49bb",
                                 "#7e4bb7",
                                 "#6bcde1",
                                 "#516a7a",
                                 "#ffa85a",
                                 "#fadc13",
                                 "#743217",
                                 "#5a8727",
                                 "#be4b18",
                                 "#33b51b",
                                 "#de2e1e",
                                 "#28988d",
                                 "#4c6db8",
                                 "#222957")
  )+
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Difference in proportion of active vs inactive elements")
p


cowplot::save_plot(paste0('figures/s16a.pdf'),
                   p,
                   base_height = 6,
                   base_width = 5,
                   device = cairo_pdf
)

# save off supplementary table 26
tmp %>%
  ungroup() %>% 
  dplyr::group_by(component, active_any) %>%
  #dplyr::mutate(num_in_group = sum(num_in_group),
  #                 tot_in_dhs = sum(tot_in_dhs)) %>%
  #ungroup() %>% 
  pivot_wider(id_cols = c(cohort, component), names_from = c(active_any), values_from = c(num_in_group, tot_in_dhs)) %>%
  dplyr::group_by(cohort, component) %>% 
  dplyr::summarize(diff = PropCIs::wald2ci(num_in_group_TRUE, tot_in_dhs_TRUE, num_in_group_FALSE, tot_in_dhs_FALSE, conf.level = 0.95, adjust = "AC")$estimate,
                   diff_lower = PropCIs::wald2ci(num_in_group_TRUE, tot_in_dhs_TRUE, num_in_group_FALSE, tot_in_dhs_FALSE, conf.level = 0.95, adjust = "AC")$conf.int[1],
                   diff_upper = PropCIs::wald2ci(num_in_group_TRUE, tot_in_dhs_TRUE, num_in_group_FALSE, tot_in_dhs_FALSE, conf.level = 0.95, adjust = "AC")$conf.int[2]) %>%
  vroom::vroom_write('tables/stable26.txt', delim = '\t')

