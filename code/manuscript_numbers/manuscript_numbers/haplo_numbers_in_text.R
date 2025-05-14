
source("code/utils.R")
# Load libraries
library(tidyverse)
library(binom)
library(vroom)
options(stringsAsFactors = FALSE)

haplos_new_df <- vroom('/mnt/sdb/gwas_eqtl_mpra/data/haplos/haplos_master_table.txt.gz')
haplos_meta <- vroom('/mnt/sdb/gwas_eqtl_mpra/data/haplos/haplos_CS_annots_NEW.txt.gz')
ints <- vroom('/mnt/sdb/gwas_eqtl_mpra/data/haplos/haplos_ints.txt.gz')

# Count number of pairs (number in text)
haplos_new_df %>%
  distinct(v1, v2)
  


# how many haplos have each variant in the pair with allele-specific regulatory activity
# in different CSs
#(number in text)
haplos_new_df %>% distinct(v1,v2,center_var, window, library, altref_emVar, refalt_emVar) %>% 
  group_by(v1,v2) %>%
  dplyr::mutate(altref_emVar_any = any(altref_emVar)) %>%
  dplyr::mutate(refalt_emVar_any = any(refalt_emVar)) %>%
  ungroup() %>% 
  distinct(v1,v2,refalt_emVar_any, altref_emVar_any) %>% 
  dplyr::filter(refalt_emVar_any == T, altref_emVar_any == T) %>% 
  left_join(haplos_meta, by = c("v1","v2")) %>%
  dplyr::filter(any_same_CS_cohort==F) %>%
  distinct(v1,v2)

# how many haplos in a CRE with at least 1 emVar?
haplos_new_df %>% dplyr::filter(any_emVar_meta == T) %>% distinct(v1,v2)

## checking interacting emVar number (number in text)
haplos_new_df %>% 
  dplyr::filter(int_emVar_meta_new == T) %>% 
  distinct(v1, v2)

# amplifying vs dampening
ints %>% count(int_type)
binom.test(139, 180, 0.5, alternative="two.sided") 
