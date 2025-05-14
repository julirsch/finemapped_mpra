## produces /mnt/sdb/gwas_eqtl_mpra/data/haplos/haplos_CS_annots.txt.gz
## supp figure 12 j,k,l
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
library(ggpointdensity)
library(ggpubr)
library(ggrastr)
library(circlize)
library(ComplexHeatmap)
library(metafor)
library(GenomicRanges)
options(stringsAsFactors = FALSE)

# read in mpra
mpra <- vroom('data/preprocess/core_mpra.txt.gz')
haplos_new_df <- vroom('data/haplos/haplos_master_table.txt.gz')

# read in credible set metadata
files <- list.files("/mnt/sdb/gtex_mpra/data/merged_CSs/", pattern = "*.txt.gz", full.names = T)
csm_full <- vroom::vroom(files, col_names = T)

# get credible set info 
mpra_csm <- mpra %>% 
  dplyr::filter(pip > 0.1) %>% 
  ungroup() %>%
  dplyr::mutate(uid = case_when( cohort %in% c('UKBB','BBJ') ~ paste0(cohort,';',trait, ";", region, ";", cs_id),
                                 cohort == 'GTEx' ~ paste0(cohort,';',tissue, ";", gene, ";", cs_id),
                                 T ~ NA)) %>% 
  dplyr::select(uid, trait, region, tissue, gene, cs_id, variant, pip, cohort, CRE) %>%
  distinct() %>% 
  left_join( csm_full, by = c("uid"='id'))

# merge in credible set info for haplotype variants
out_meta_info <-haplos_new_df %>%
  distinct(v1,v2) %>%
  ungroup() %>%
  left_join(mpra_csm,
             by = c("v1" = "variant")) %>%
  left_join(mpra_csm,
             by = c("v2" = "variant")) %>%
  # throw out anything with no credible set ID - just for credible set corroboration
  dplyr::filter(!is.na(csm_id.x)) %>%
  dplyr::filter(!is.na(csm_id.y)) %>% 
  dplyr::mutate(sameCS = ifelse(csm_id.x == csm_id.y, T, F)) %>% 
  ## is the pair of variants in the same CS?
  dplyr::group_by(v1,v2) %>%
  dplyr::mutate(any_same_CS = ifelse(any(sameCS == T, na.rm = T), T, F)) %>%
  # if every pair contains at least one NA for csm, mark it as NA
  dplyr::mutate(any_same_CS = ifelse(all(is.na(sameCS)), NA, any_same_CS)) %>% 
  ungroup()

out_meta_info <- out_meta_info %>% 
  dplyr::mutate(cohort.v1 = case_when(cohort.x %in% c('UKBB','BBJ') ~ 'Complex Trait',
                                                      cohort.x == 'GTEx' ~ 'eQTL',
                                                      T ~ NA)) %>% 
  dplyr::mutate(cohort.v2 = case_when(cohort.y %in% c('UKBB','BBJ') ~ 'Complex Trait',
                                      cohort.y == 'GTEx' ~ 'eQTL',
                                      T ~ NA)) %>%
  dplyr::mutate(same_CS_cohort = ifelse(cohort.v1 == cohort.v2 & sameCS==T, T, F)) %>%
  dplyr::group_by(v1,v2) %>%
  dplyr::mutate(any_same_CS_cohort = ifelse(any(same_CS_cohort == T , na.rm = T), T, F)) %>%
  # if every pair contains at least one NA for csm, mark it as NA
  dplyr::mutate(any_same_CS_cohort = ifelse(all(is.na(same_CS_cohort)), NA, any_same_CS_cohort)) %>% 
  ungroup() %>%
  dplyr::mutate(cohort = ifelse(cohort.v1 == cohort.v2, cohort.v1, NA)) 

# pairs that are only cross UKBB/GTEx:
out_meta_info %>% 
  distinct(v1,v2,cohort.v1, cohort.v2, cohort, any_same_CS, sameCS) %>% 
  dplyr::rename(overall_cohort = cohort) %>% 
  group_by(v1,v2) %>% dplyr::filter(all(is.na(overall_cohort))) %>% 
  ungroup() %>%
  dplyr::filter(!is.na(any_same_CS))
  # left_join(mpra %>% distinct(variant, cohort, pip), by = c("v1" = "variant") ) %>% 
  # left_join(mpra %>% distinct(variant, cohort, pip), by = c("v2" = "variant") ) %>%
  # print(n=500)

# save out_meta_info
out_meta_info %>% vroom::vroom_write('data/haplos/haplos_CS_annots_NEW.txt.gz',delim = "\t")

# variants in credible sets vs not  (number in text)
out_meta_info %>% 
  distinct(v1,v2,any_same_CS) %>%
  count(any_same_CS)

# supplemental figure of proportion of variant pairs in same or diff CS's that have a variant as an emVar
# 59 pairs in both cohorts - need to make sure we have any_same_CS specific to the cohort
df1 <- out_meta_info %>% distinct(v1,v2,cohort, any_same_CS_cohort) %>% 
  dplyr::filter(!is.na(cohort)) %>% 
  left_join(haplos_new_df %>% distinct(v1,v2,any_emVar_meta), by = c('v1','v2')) %>%  
  na.omit() %>% 
  group_by(v1,v2,cohort, any_same_CS_cohort) %>%
  dplyr::mutate(any_emVar_meta = any(any_emVar_meta)) %>%
  ungroup() %>% 
  distinct(v1,v2,cohort, any_same_CS_cohort, any_emVar_meta) %>% 
  count(cohort, any_same_CS_cohort, any_emVar_meta) %>%
  group_by(cohort, any_same_CS_cohort) %>%
  dplyr::mutate(tot = sum(n)) %>%
  dplyr::mutate(prop = n / tot,
                lower = binom.confint(n, tot, method='wilson')$lower,
                upper = binom.confint(n, tot, method='wilson')$upper) %>%
  dplyr::filter(any_emVar_meta == T)

df1
p1 <- df1 %>% ggplot(., aes(x = any_same_CS_cohort, y = prop, fill = cohort, group = cohort)) +
  geom_col(position = 'dodge')+
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                linewidth = 0.5, width = 0, col = "black",
                position = position_dodge(0.9, preserve = "single")) +
  pretty_plot()+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c(jdb_palette('brewer_spectra')[c(1,2)]),
                    aesthetics = c("colour", "fill")) +
  ylab("Proportion of emVars") +
  xlab("Same Credible Set")+
  theme(legend.position = "none")+
  ylim(0,0.8)
p1

cowplot::save_plot(
  'figures/s12j.pdf',
  p1,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf)

# proportion of variant pairs in same or diff CS that are interacting
# of variant pairs that have 1 emVar, what proportion in same/diff CS's are interacting?
df2 <- out_meta_info %>% distinct(v1,v2,any_same_CS_cohort, cohort) %>% 
  dplyr::filter(!is.na(cohort)) %>% 
  left_join(haplos_new_df %>% 
              group_by(v1,v2,cell_type, any_emVar_meta, int_emVar_meta_new) %>% 
              # int_emvar_meta_any = is this pair an interacting emVar in any window/center_var/library
              dplyr::mutate(int_emVar_meta_any = any(int_emVar_meta_new)) %>% 
              ungroup() %>%
              distinct(v1,v2,cell_type, int_emVar_meta_any, any_emVar_meta),
            by = c('v1','v2')) %>% 
  group_by(v1,v2) %>%
  dplyr::filter(any(any_emVar_meta == T)) %>%
  ungroup() %>% 
  group_by(v1,v2, cohort, any_same_CS_cohort) %>% 
  # int_emvar_meta_any_ct = is this pair an interacting emVar in any cell type
  dplyr::mutate(int_emVar_meta_any_ct = any(int_emVar_meta_any)) %>% 
  ungroup() %>% 
  distinct(v1,v2,int_emVar_meta_any_ct, cohort, any_same_CS_cohort) %>% 
  count(int_emVar_meta_any_ct, cohort, any_same_CS_cohort) %>%
  group_by(cohort, any_same_CS_cohort) %>%
  dplyr::mutate(tot = sum(n)) %>%
  dplyr::mutate(prop = n / tot,
                lower = binom.confint(n, tot, method='wilson')$lower,
                upper = binom.confint(n, tot, method='wilson')$upper) %>%
  dplyr::filter(int_emVar_meta_any_ct == T)

df2

p2 <- df2 %>% ggplot(., aes(x = any_same_CS_cohort, y = prop, fill = cohort, group = cohort)) +
  geom_col(position = 'dodge')+
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                size = 0.5, width = 0, col = "black",
                position = position_dodge(0.9, preserve = "single")) +
  pretty_plot()+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c(jdb_palette('brewer_spectra')[c(1,2)]),
                    aesthetics = c("colour", "fill")) +
  ylab("Proportion of emVars") +
  xlab("Same Credible Set")+
  theme(legend.position = "none")+
  ylim(0,0.8)
p2

cowplot::save_plot(
  'figures/s12k.pdf',
  p2,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf)


haplos_new_df %>% distinct(v1,v2,center_var, window, library, altref_emVar, refalt_emVar) %>% 
  group_by(v1,v2) %>%
  dplyr::mutate(altref_emVar_any = any(altref_emVar)) %>%
  dplyr::mutate(refalt_emVar_any = any(refalt_emVar)) %>%
  ungroup() %>% 
  distinct(v1,v2,refalt_emVar_any, altref_emVar_any) %>% 
  dplyr::filter(refalt_emVar_any == T, altref_emVar_any == T) %>% 
  left_join(out_meta_info, by = c("v1","v2")) %>%
  dplyr::filter(any_same_CS_cohort==F) %>%
  distinct(v1,v2)


# interacting emVars
ints <- vroom('/mnt/sdb/gwas_eqtl_mpra/data/haplos/haplos_ints.txt.gz')
p <- left_join(ints, out_meta_info %>% distinct(v1, v2, cohort, any_same_CS_cohort), by = c('v1','v2')) %>%
  count(cohort, any_same_CS_cohort, int_type) %>%
  na.omit() %>% 
  dplyr::mutate(group = paste(cohort, any_same_CS_cohort, sep = ' ')) %>% 
  ggplot(., aes(x=group, y=n, fill = int_type)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c('red','blue'))+
  xlab('same CS') + 
  pretty_plot()+
  theme(legend.position = "none") 

p

cowplot::save_plot(
  'figures/s12l.pdf',
  p,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf)
