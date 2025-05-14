# Load in libraries
source("code/utils.R")
library(vroom)
library(tidyverse)
library(BuenColors)
library(ggplot2)
library(dplyr)
library(binom)
library(ggpubr)
options(stringsAsFactors = FALSE)

# Read in processed MPRA data
mpra_df <- vroom('data/preprocess/core_mpra.txt.gz') %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0))


# Read in spliceAI data
vcf2 <- data.table::fread('/mnt/sdb/gwas_eqtl_mpra/reviews/data/spliceAI2/0001.vcf') %>%
  dplyr::mutate(`#CHROM` = paste0('chr',`#CHROM`)) %>%
  dplyr::rename(CHROM = `#CHROM`) %>%
  dplyr::mutate(variant = paste(CHROM, POS, REF, ALT, sep = ':')) %>%
  separate(.,col = 'INFO', into = c('ALLELE','SYMBOL','AG_score','AL_score','DG_score','DL_score','AG_pos','AL_pos','DG_pos','DL_pos'),
           sep = '\\|',
           remove = F)
vcf2$AG_score <- as.numeric(vcf2$AG_score)
vcf2$AL_score <- as.numeric(vcf2$AL_score)
vcf2$DG_score <- as.numeric(vcf2$DG_score)
vcf2$DL_score <- as.numeric(vcf2$DL_score)


# Look at high-PIP variant
df <- mpra_df %>% 
  dplyr::filter(pip > 0.9, cohort != 'control') %>%
  dplyr::mutate(cohort = case_when(cohort =='GTEx' ~'eQTLs',
                                   cohort %in% c('BBJ','UKBB') ~ 'Complex Traits',
                                   T ~ NA)) %>%
  group_by(variant, cohort) %>%
  dplyr::mutate(status = case_when(all(active == F) & all(emVar == F) ~ 'inactive',
                                   any(active == T) & all(emVar == F) ~ 'active non-emVar',
                                   any(active == T) & any(emVar == T) ~ 'emVar',
                                   T ~ NA)) %>%
  ungroup() %>% 
  left_join(vcf2 %>%
              dplyr::select(-ID,-QUAL,-FILTER,-INFO, -ends_with('pos')) %>%
              pivot_longer(cols = c('AG_score','AL_score','DG_score','DL_score'), 
                           names_to = 'spliceAI_type',
                           values_to ='score') %>%
              group_by(variant) %>%
              dplyr::mutate(maxscore = max(score)) %>%
              pivot_wider(names_from = spliceAI_type, values_from = score) %>%
              dplyr::filter(maxscore > 0.1),
            by = 'variant')

# Assign consequence and count
df2 <- df %>% 
  distinct(variant, status, consequence, consequence2, maxscore) %>%
  group_by(variant) %>% 
  dplyr::mutate(consequence3 = case_when(any(!is.na(maxscore)) ~ 'SpliceAI',
                                         any(consequence == 'LoF') | any(consequence == 'missense') ~ 'LoF / missense',
                                         any(consequence == 'UTR3') ~ 'UTR3',
                                         any(consequence == 'UTR5') ~ 'UTR5',
                                         any(consequence2 =='promoter') ~ 'promoter',
                                         any(consequence2 == 'CRE') ~ 'Distal CRE',
                                         T ~ NA)
  ) %>%
  ungroup() %>%
  distinct(variant, status, consequence3) %>%
  dplyr::count(status, consequence3) %>%
  group_by(status) %>%
  dplyr::mutate(tot=sum(n),
                prop = n/tot,
                upper= binom.confint(n, tot, method="wilson")$upper,
                lower= binom.confint(n, tot, method="wilson")$lower) %>%
  ungroup() %>%
  dplyr::filter(!is.na(consequence3)) 

# get TARV proportion
df3 <-  df %>% 
  distinct(variant, status, consequence, consequence2, maxscore) %>%
  group_by(variant) %>% 
  dplyr::mutate(consequence3 = case_when(any(!is.na(maxscore)) ~ 'SpliceAI',
                                         any(consequence == 'LoF') | any(consequence == 'missense') ~ 'LoF / missense',
                                         any(consequence == 'UTR3') ~ 'UTR3',
                                         any(consequence == 'UTR5') ~ 'UTR5',
                                         any(consequence2 =='promoter') ~ 'promoter',
                                         any(consequence2 == 'CRE') ~ 'Distal CRE',
                                         T ~ NA)
  ) %>%
  ungroup() %>%
  dplyr::mutate(consequence3 = ifelse(consequence3 %in% c('promoter','Distal CRE'), 'Any CRE', NA)) %>%
  distinct(variant, status, consequence3) %>%
  dplyr::count(status, consequence3) %>%
  group_by(status) %>%
  dplyr::mutate(tot=sum(n),
                prop = n/tot,
                upper= binom.confint(n, tot, method="wilson")$upper,
                lower= binom.confint(n, tot, method="wilson")$lower) %>%
  ungroup() %>%
  dplyr::filter(!is.na(consequence3)) 

df_overall <- rbind(df2,df3)

# Prepare for plotting 
df_overall$status <- factor(df_overall$status, levels= c('inactive','active non-emVar','emVar'))
df_overall$consequence3 <- factor(df_overall$consequence3, levels= c('LoF / missense','SpliceAI','UTR3','UTR5','promoter','Distal CRE','Any CRE'))

# Plot
p1 <- ggplot(df_overall,aes(x=consequence3, y= prop, group = status, fill = status)) +
  geom_col(position = 'dodge')+
  scale_fill_manual(values = c('#8376BC',
                               '#FCCCB1',
                               BuenColors::jdb_palette('brewer_spectra')[7]))+
  geom_errorbar(aes(ymin = lower, ymax = upper),position = position_dodge(0.9), width = 0)+
  pretty_plot()+
  ylab('proportion of variants PIP > 0.9')+
  xlab('Annotation')+
  ggtitle('High PIP (PIP > 0.9) variant breakdown')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
p1

# Save
cowplot::save_plot(paste0('figures/s4c.pdf'),
                   p1,
                   base_height = 6,
                   base_width = 6,
                   device = cairo_pdf
)

# Now, examining credible sets
df_cs <- mpra_df %>% 
  dplyr::filter(cs_id != -1, cohort != 'control') %>%
  dplyr::mutate(uid = paste(cohort,tissue, trait, gene, region, cs_id, sep = ';')) %>%
  group_by(uid, cohort) %>%
  dplyr::mutate(status = case_when(all(active == F) & all(emVar == F) ~ 'inactive',
                                   any(active == T) & all(emVar == F) ~ 'active non-emVar',
                                   any(active == T) & any(emVar == T) ~ 'emVar',
                                   T ~ NA)) %>%
  ungroup() %>% 
  left_join(vcf2 %>%
              dplyr::select(-ID,-QUAL,-FILTER,-INFO, -ends_with('pos')) %>%
              pivot_longer(cols = c('AG_score','AL_score','DG_score','DL_score'), 
                           names_to = 'spliceAI_type',
                           values_to ='score') %>%
              group_by(variant) %>%
              dplyr::mutate(maxscore = max(score)) %>%
              pivot_wider(names_from = spliceAI_type, values_from = score) %>%
              dplyr::filter(maxscore > 0.1),
            by = 'variant') %>%
  dplyr::mutate(TARV = ifelse(emVar == T & CRE> 0 ,1,0)) %>%
  distinct(variant, uid, status, consequence, consequence2, maxscore, TARV) 

# assign consequence per UID 
df_cs2 <- df_cs %>%
  dplyr::mutate(consequence3 = case_when(!is.na(maxscore) ~ 'SpliceAI',
                                         consequence %in% c('LoF','missense')~ 'LoF / missense',
                                         consequence == 'UTR3' ~ 'UTR3',
                                         consequence == 'UTR5' ~ 'UTR5',
                                         consequence2 =='promoter' ~ 'promoter',
                                         consequence2 == 'CRE' ~ 'Distal CRE',
                                         T ~ NA)) %>%
  group_by(uid) %>%
  dplyr::mutate(`LoF / missense` = ifelse(any(consequence3 == 'LoF / missense'),1,0),
                promoter = ifelse(any(consequence3 == 'promoter'),1,0),
                `Distal CRE` = ifelse(any(consequence3 == 'Distal CRE'),1,0),
                UTR5=ifelse(any(consequence3 == 'UTR5'),1,0),
                UTR3=ifelse(any(consequence3 == 'UTR3'),1,0),
                spliceAI =ifelse(any(consequence3 == 'SpliceAI'),1,0),
                TARV_cs = ifelse(any(TARV > 0), 1,0)) %>%
  dplyr::mutate(`Any CRE` = ifelse((`Distal CRE` >0 | promoter>0), 1,0)) %>%
  ungroup() %>%
  distinct(uid, status, spliceAI, `LoF / missense`, UTR3, UTR5, promoter, `Distal CRE`, TARV_cs, `Any CRE`) %>%
  replace(is.na(.),0) %>%
  dplyr::mutate(val = 1) %>%
  pivot_wider(id_cols = c(uid, spliceAI, `LoF / missense`, UTR3, UTR5, promoter, `Distal CRE`, TARV_cs, `Any CRE`),
              names_from = 'status',
              values_from = val) %>%
  replace(is.na(.),0) 

# count consequence
df_cs3 <- df_cs2 %>%
  pivot_longer(cols = c(spliceAI, `LoF / missense`, UTR3, UTR5, promoter, `Distal CRE`, `Any CRE`, TARV_cs, inactive, `active non-emVar`,emVar),
               names_to = 'annot',
               values_to = 'n') %>%
  dplyr::mutate(annot = ifelse(annot == 'TARV_cs','TARV',annot)) %>%
  dplyr::filter(n==1) %>%
  count(annot) %>%
  dplyr::mutate(prop = n/dim(df_cs2)[1],
                upper= binom.confint(n, dim(df_cs2)[1], method="wilson")$upper,
                lower= binom.confint(n, dim(df_cs2)[1], method="wilson")$lower)

# Prepare for plotting
df_cs3 <- df_cs3 %>% dplyr::mutate(c = ifelse(annot %in% c('inactive','active non-emVar','emVar','TARV'),annot,'CS'))

df_cs3$annot <- factor(df_cs3$annot, levels= c('LoF / missense','spliceAI','UTR3','UTR5','promoter','Distal CRE', 'Any CRE','inactive','active non-emVar','emVar','TARV'))

#Plot
#EB5C0A for tarvs
##584A96 for the rest of the CS
p2 <-ggplot(df_cs3,aes(x=annot, y= prop, fill = c)) +
  geom_col(position = 'identity')+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
  pretty_plot()+
  scale_fill_manual(values = c('#FCCCB1',
                               '#584A96',
                               BuenColors::jdb_palette('brewer_spectra')[7],
                               '#8376BC',
                               '#EB5C0A'))+
  ylab('proportion of credible sets')+
  xlab('Annotation')+
  ggtitle('Credible set breakdown')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
p2

# Save
cowplot::save_plot(paste0('figures/s4d.pdf'),
                   p2,
                   base_height = 6,
                   base_width = 6,
                   device = cairo_pdf
)


plt_combined <- ggarrange(p1, p2,
                          labels = c("c","d"),
                          ncol = 2, nrow = 1)

plt_combined

cowplot::save_plot(paste0('figures/s4cd.pdf'),
                   plt_combined,
                   base_height = 6,
                   base_width = 15,
                   device = cairo_pdf
)

# p val for UTR3
prop.test(matrix(c(c(571,167+127),c(9562-571, 3817+3796-127-167)), nrow = 2, ncol = 2), alternative = 't')

