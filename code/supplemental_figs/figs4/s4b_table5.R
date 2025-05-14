
# Load in libraries
source("code/utils.R")
library(vroom)
library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(plyranges)
library(BuenColors)
library(rtracklayer)
library(reshape2)
library(ComplexHeatmap)
library(binom)
library(ggExtra)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE)

# read in data

chip <- vroom('/mnt/sdb/gwas_eqtl_mpra/data/annotations/motif/variant_chip_lite.txt.gz') %>%
  dplyr::mutate(instance = 1) %>%
  pivot_wider(names_from = antigen, values_from = instance) %>%
  replace(is.na(.), 0)

mpra_df <- vroom("data/preprocess/core_mpra.txt.gz") 

# keep the test variants
mpra_df <-
  mpra_df %>% 
  dplyr::filter(pip > 0.1 | cs_id > 0 ) %>%
  dplyr::filter(cohort != 'control')

#check for completeness
mpra_df %>% distinct(variant, active_any)
chip %>% distinct(variant)
mpra_df %>% 
  distinct(variant,CRE) %>%
  dplyr::filter(variant %ni% (chip %>% distinct(variant) %>% .$variant)) %>% dplyr::count(CRE > 0)


mpra_lite <- mpra_df %>%
  distinct(variant, active_any) %>%
  inner_join(.,chip, by = "variant")

#mpra_lite %>% dplyr::filter(is.na(ESR1))

df <- mpra_lite %>% 
  dplyr::select(-variant) %>% 
  gather(key, value, -active_any) %>% 
  xtabs(value ~ key + active_any, data = .)
dimnames(df)$active_any <- c("yesTF_noActive","yesTF_yesActive")
df <- as.data.frame.matrix(df)

# set up the complete columns for odds ratio
noactive=mpra_lite %>% dplyr::filter(active_any==FALSE) %>% nrow()
yesactive=mpra_lite %>% dplyr::filter(active_any==TRUE) %>% nrow()
df <- df %>% dplyr::mutate(noTF_noActive = noactive - yesTF_noActive,
                           noTF_yesActive = yesactive - yesTF_yesActive)

#method 1: using formulas
df2 <- df  %>% 
  dplyr::select(yesTF_yesActive, noTF_yesActive, yesTF_noActive, noTF_noActive) %>% 
  dplyr::mutate(OR = (noTF_noActive * yesTF_yesActive) / (noTF_yesActive * yesTF_noActive)) %>%
  dplyr::mutate(lower.CI = exp(log(OR) - 1.96 * sqrt((1/noTF_noActive) + (1/noTF_yesActive) + (1/yesTF_noActive) + (1/yesTF_yesActive))),
                upper.CI = exp(log(OR) + 1.96 * sqrt((1/noTF_noActive) + (1/noTF_yesActive) + (1/yesTF_noActive) + (1/yesTF_yesActive)))) %>%
  dplyr::mutate(sig = ifelse((lower.CI < 1 & upper.CI < 1) | (lower.CI > 1 & upper.CI > 1), T, F))

#method2: using fisher exact test
fisher_estimate <- function(a){
  a =  fisher.test(matrix(c(a[1],a[2],a[3],a[4]), nrow =2, ncol = 2, byrow = T))
  return(a$estimate)
}
fisher_lower_CI <- function(a){
  a =  fisher.test(matrix(c(a[1],a[2],a[3],a[4]), nrow =2, ncol = 2, byrow = T))
  return(a$conf.int[1])
}
fisher_upper_CI <- function(a){
  a =  fisher.test(matrix(c(a[1],a[2],a[3],a[4]), nrow =2, ncol = 2, byrow = T))
  return(a$conf.int[2])
}
fisher_p <- function(a){
  a =  fisher.test(matrix(c(a[1],a[2],a[3],a[4]),  nrow =2, ncol = 2, byrow = T))
  return(a$p)
}


df3 <- df %>% 
  dplyr::mutate(OR = apply(df %>% dplyr::select(yesTF_yesActive, noTF_yesActive, yesTF_noActive, noTF_noActive),
                           MARGIN = 1,
                           FUN = fisher_estimate),
                lower.CI = apply(df %>% dplyr::select(yesTF_yesActive, noTF_yesActive, yesTF_noActive, noTF_noActive),
                                 MARGIN = 1,
                                 FUN = fisher_lower_CI),
                upper.CI = apply(df %>% dplyr::select(yesTF_yesActive, noTF_yesActive, yesTF_noActive, noTF_noActive),
                                 MARGIN = 1,
                                 FUN = fisher_upper_CI),
                p = apply(df %>% dplyr::select(yesTF_yesActive, noTF_yesActive, yesTF_noActive, noTF_noActive),
                          MARGIN = 1,
                          FUN = fisher_p))
df3 <- tibble::rownames_to_column(df3, "TF")

df3 <- df3 %>% 
  dplyr::mutate(p.corrected = p.adjust(df3$p,method = 'bonferroni')) %>%
  dplyr::filter( (yesTF_yesActive + yesTF_noActive) > 20) %>%
  dplyr::mutate(color = ifelse(p.corrected < 0.05, "sig","nonsig")) %>% dplyr::mutate(color = replace_na(color,"nonsig"))

# Plot
p <-ggplot(df3, aes(x=reorder(TF, OR), y=OR), color = color)+
  #geom_point(aes(color=color))+
  #geom_point(data= df3%>> arrang %>% dplyr::filter(color == "nonsig"), aes(fill = "black", color = color))+
  #geom_point(data=df3 %>% dplyr::filter(color == "sig"), aes(fill = "black", color = color))+
  # geom_errorbar(aes(ymin=lower.CI, ymax=upper.CI, color = color), width=1)+
  geom_pointrange(aes(ymin = lower.CI, ymax = upper.CI, color = color), fatten = 0.01, alpha = 0.5) + 
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)
  scale_color_manual(values = c("grey","maroon"))+
  pretty_plot()+
  ylim(0,30)+
  theme(legend.position="none")
plt_combined <- p + geom_hline(yintercept=1, linetype="dashed")

# Save plot
plt_combined
cowplot::save_plot("figures/s4b.pdf",
                   plt_combined,
                   base_height = 6,
                   base_width = 6,
                   device = cairo_pdf
)

# save off Supplementary table
df3 %>% vroom::vroom_write('tables/stable5.txt.gz', delim = '\t',col_names = T)


