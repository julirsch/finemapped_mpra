# Set temp storage
unixtools:: set.tempdir("/mnt/sdb/rstudio/sirajl/")

# Set working directory
setwd("/mnt/sdb/gwas_eqtl_mpra/code/github_ready")

# Load libraries
source("code/utils.R")
library(vroom)
library(tidyverse)

# Read in processed MPRA results
mpra_df <- vroom('data/preprocess/core_mpra.txt.gz') 

## Lines 156-160
# How many variants in complex traits, eqtls, or both
mpra_df %>% 
  dplyr::filter(cohort != 'control') %>%
  dplyr::filter((pip > 0.1 | cs_id > 0)) %>%
  distinct(variant, cohort) %>% 
  dplyr::mutate(count = 1) %>%
  pivot_wider(names_from = cohort, values_from = count) %>%
  dplyr::mutate(UKBB = ifelse(is.na(UKBB), 0, UKBB),
                GTEx = ifelse(is.na(GTEx), 0, GTEx),
                BBJ = ifelse(is.na(BBJ), 0, BBJ)) %>% 
  dplyr::mutate(complex = ifelse(UKBB > 0 |BBJ > 0, 1, 0)) %>%
  dplyr::select(-UKBB,BBJ) %>%
  dplyr::count(GTEx, complex)
# How many complex traits?
mpra_df %>% 
  dplyr::filter(pip > 0.1 | cs_id > 0 | cohort == 'control') %>% 
  distinct(trait) %>%
  dplyr::filter(!is.na(trait))
# How many genes and tissues?
mpra_df %>% 
  dplyr::filter(pip > 0.1 | cs_id > 0 | cohort == 'control') %>% 
  distinct(gene) %>%
  dplyr::filter(!is.na(gene))
mpra_df %>% 
  dplyr::filter(pip > 0.1 | cs_id > 0 | cohort == 'control') %>% 
  distinct(tissue) %>%
  dplyr::filter(!is.na(tissue))

## Line 161
# Number of credible sets
mpra_df %>% 
  dplyr::filter(pip > 0.1 | cs_id > 0 | cohort == 'control') %>% 
  dplyr::filter(cohort != 'control') %>% 
  distinct(variant, cs_uid, pip) %>%
  dplyr::filter(!is.na(cs_uid)) %>%
  group_by(cs_uid) %>%
  dplyr::mutate(total_pip = sum(pip)) %>%
  distinct(total_pip, cs_uid) %>%
  ungroup() %>%
  count(total_pip > 0.8)
# Number of unique high-PIP variants (PIP > 0.5)
mpra_df %>% 
  dplyr::filter(pip > 0.1 | cs_id > 0 | cohort == 'control') %>% 
  dplyr::filter(cohort != 'control') %>% 
  dplyr::filter(pip > 0.5) %>%
  distinct(variant)

## Lines 162-165
# Controls breakdoown
mpra_df %>%
  dplyr::filter(cohort == 'control') %>%
  distinct(variant)
# Null trait controls
mpra_df %>% 
  dplyr::filter(cohort == 'control') %>%
  dplyr::filter(grepl('null',type)) %>%
  distinct(variant)
# Location controls
mpra_df %>% 
  dplyr::filter(cohort == 'control') %>%
  dplyr::filter(grepl('loc',type)) %>%
  distinct(variant)
# Annotation controls
mpra_df %>% 
  dplyr::filter(cohort == 'control') %>%
  dplyr::filter(grepl('annot',type)) %>%
  distinct(variant)

## Line 167
# How many overall distinct variants
mpra_df %>% 
  dplyr::filter(pip > 0.1 | cs_id > 0 | cohort == 'control') %>% 
  distinct(variant) 

## Lines 174-177
# How many active
mpra_df %>%
  dplyr::filter(pip > 0.1 | cs_id > 0 | cohort == 'control') %>% 
  dplyr::filter(active_any==T) %>%
  distinct(variant)
# How many active repressive
mpra_df %>% 
  dplyr::filter(pip > 0.1 | cs_id > 0 | cohort == 'control') %>% 
  dplyr::filter(active==T, log2FC < -1) %>%
  distinct(variant) %>% count()
# Test variants: emVar vs active
mpra_df %>%
  dplyr::filter(pip > 0.1 | cs_id > 0) %>% 
  dplyr::filter(cohort != 'control') %>%
  distinct(emVar_any, active_any, variant) %>%
  dplyr::count(emVar_any,active_any)
# Controls - emVar vs active
mpra_df %>%
  dplyr::filter(cohort == 'control') %>%
  distinct(emVar_any, active_any, variant) %>%
  dplyr::count(emVar_any,active_any)
# Total - emVar vs active
mpra_df %>%
  dplyr::filter(cohort == 'control' | pip > 0.1 | cs_id > 0) %>%
  distinct(emVar_any, active_any, variant) %>%
  dplyr::count(emVar_any,active_any)

## Line 213
# Number of emVars among active variants
mpra_df %>% 
  dplyr::filter(cohort == 'control' | pip > 0.1 | cs_id > 0) %>%
  dplyr::filter(active_any == T) %>%
  distinct(variant, emVar_any) %>%
  count(emVar_any==T)


# Effect sizes of emVars overall
mpra_df %>% 
  dplyr::filter(pip > 0.1 | cs_id > 0 | cohort == 'control') %>% 
  # get specific lines where they are active and emVars
  dplyr::filter(active==T, emVar == T) %>%
  distinct(variant, cell_type, log2Skew) %>%
  dplyr::mutate(log2Skew_mag = abs(log2Skew)) %>%
  group_by(cell_type) %>%
  dplyr::mutate(med_skew = median(abs(log2Skew))) %>%
  ungroup() %>%
  distinct(cell_type, med_skew) %>%
  dplyr::summarize(mean = mean(med_skew))
# Effect sizes of controls / test
mpra_df %>% 
  dplyr::filter(cohort == 'control') %>% 
  # get specific lines where they are active and emVars
  dplyr::filter(active==T, emVar == T) %>%
  distinct(variant, cell_type, log2Skew) %>%
  dplyr::mutate(log2Skew_mag = abs(log2Skew)) %>%
  group_by(cell_type) %>%
  dplyr::mutate(med_skew = median(abs(log2Skew))) %>%
  ungroup() %>%
  distinct(cell_type, med_skew) %>%
  dplyr::summarize(mean = median(med_skew))
mpra_df %>% 
  dplyr::filter(pip > 0.1 | cs_id > 0) %>% 
  dplyr::filter(cohort != 'control') %>% 
  # get specific lines where they are active and emVars
  dplyr::filter(active==T, emVar == T) %>%
  distinct(variant, cell_type, log2Skew) %>%
  dplyr::mutate(log2Skew_mag = abs(log2Skew)) %>%
  group_by(cell_type) %>%
  dplyr::mutate(med_skew = median(abs(log2Skew))) %>%
  ungroup() %>%
  distinct(cell_type, med_skew) %>%
  dplyr::summarize(mean = median(med_skew))



# table 2
t2 <- vroom('/mnt/sdb/gwas_eqtl_mpra/tables/stable2.txt.gz')
# unique oligos:
t2 %>%
  dplyr::select(-duplicate) %>% 
  pivot_longer(cols = c('allele1_oligo','allele2_oligo')) %>%
  distinct(value)

# table 3
#t3 <- vroom('/mnt/sdb/gwas_eqtl_mpra/tables/stable3_v2.txt.gz')
t3 <- vroom('/mnt/sdb/gwas_eqtl_mpra/tables/stable3_v2_unpivoted.txt.gz')

t3 %>% 
  distinct(variant)
