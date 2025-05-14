# Load in libraries
source("code/utils.R")
library(stats)
library(vroom)
library(tidyverse)
library(BuenColors)
library(ggplot2)
library(dplyr)
library(binom)
library(tidyr)
library(TFBSTools)
library(rtracklayer)
library(GenomicRanges)
library(genomation)
options(stringsAsFactors = FALSE)

# Read in processed MPRA data
mpra_df <- vroom('data/preprocess/core_mpra.txt.gz') %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0)) 

variant_gr <- mpra_df %>%
  distinct(variant) %>% 
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) 

variant_gr <- as.data.frame(variant_gr)
rownames(variant_gr) <- variant_gr$variant

variant_gr <- variant_gr %>% 
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "position", 
                           end.field = "position",
                           keep.extra.columns = F)


# read in ChIP files
files <- list.files("/mnt/sdb/gwas_eqtl_mpra/reviews/data/genomic_dark_matter/chip", full.names = T)
darktfs <- vroom::vroom('/mnt/sdb/gwas_eqtl_mpra/reviews/data/genomic_dark_matter/tables/darktfs.csv')
tfs <- do.call(rbind,lapply(files, function(x){strsplit(strsplit(x,'/')[[1]][9],'_')[[1]][3]}))
exp <- do.call(rbind,lapply(files, function(x){strsplit(strsplit(x,'/')[[1]][9],'_')[[1]][4]}))
files.df <- data.frame(files,tfs, exp)
files.df <- files.df %>% 
  left_join(darktfs %>% distinct(TF, `TF Type`),
            by=c("tfs" = 'TF'))
files2 <- list.files("/mnt/sdb/gwas_eqtl_mpra/reviews/data/genomic_dark_matter/chip2", full.names = T)
tfs2 <- do.call(rbind,lapply(files2, function(x){strsplit(strsplit(x,'/')[[1]][9],'_')[[1]][2]}))
exp2 <- do.call(rbind,lapply(files2, function(x){strsplit(strsplit(x,'/')[[1]][9],'_')[[1]][3]}))
files2.df <- data.frame(files2,tfs2, exp2) %>%
  dplyr::rename(files = files2,
                tfs = tfs2,
                exp = exp2)
files2.df <- files2.df %>% 
  left_join(darktfs %>% distinct(TF, `TF Type`),
            by=c("tfs" = 'TF'))

files.tot <- rbind(files.df, files2.df) %>% 
  dplyr::filter(`TF Type` == 'Dark TF') %>%
  distinct(files,tfs, exp) %>% 
  dplyr::mutate(name = paste(tfs, exp, sep = '.'))

# prepare ChIp-seq peaks, collapsing by antigen
chip.gr <- GRanges()
for (tf in (files.tot %>% distinct(tfs) %>% .$tfs)){
  print(tf)
  test <- files.tot %>% dplyr::filter(tfs == tf)
  test.gr <- lapply(test$files, import)
  test.grangeslist <- GRangesList(test.gr)
  test.gr <- reduce(unlist(test.grangeslist))
  test.gr$antigen <- tf
  chip.gr <- c(chip.gr, test.gr)
} 
chip.gr.list <- split(chip.gr, as.factor(chip.gr$antigen))
chip.df <- data.frame(lapply(1:length(chip.gr.list), 
                             function(x) {countOverlaps(variant_gr, chip.gr.list[[x]])}))
names(chip.df) <- paste("TF" ,names(chip.gr.list),sep=".")
chip.df[,1:54] <- ifelse(chip.df[,1:54] > 0, 1, 0)

# select test variants
mpra_lite <- mpra_df%>%
  distinct(variant,active_any, active, cohort) %>% 
  dplyr::group_by(variant) %>% 
  dplyr::mutate(active_any = ifelse(any(active ==T), T, F)) %>%
  ungroup() %>%
  dplyr::filter(cohort != 'control') %>%
  distinct(variant, active_any) %>% 
  left_join(tibble::rownames_to_column(chip.df, 'variant'), by='variant') 

# get odds ratios
counts_active <- mpra_lite %>%
  group_by(active_any) %>%
  dplyr::mutate(tot = n()) %>%
  distinct(active_any, tot)

active <- mpra_lite %>% 
  group_by(active_any) %>%
  dplyr::select(-variant) %>%
  summarise_all(sum) %>%
  ungroup() %>%
  left_join(counts_active, by = c('active_any')) %>%
  pivot_longer(cols = colnames(chip.df),
               names_to = 'dark.tfs',
               names_prefix = 'TF.',
               values_to = 'num') %>%
  pivot_wider(names_from = active_any, values_from = c(num, tot))

estimate <- function(num_FALSE,num_TRUE,tot_FALSE,tot_TRUE){
  e <- fisher.test(matrix(c(num_FALSE, 
                            num_TRUE,
                            tot_FALSE-num_FALSE, 
                            tot_TRUE-num_TRUE),
                          dimnames = list(c("inactive","active"),c("in ChIP", "not in ChIP")),
                          nrow=2, ncol=2))$estimate
  return(e)
}
pval <- function(num_FALSE,num_TRUE,tot_FALSE,tot_TRUE){
  e <- fisher.test(matrix(c(num_FALSE, 
                            num_TRUE,
                            tot_FALSE-num_FALSE, 
                            tot_TRUE-num_TRUE),
                          dimnames = list(c("inactive","active"),c("in ChIP", "not in ChIP")),
                          nrow=2, ncol=2))$p
  return(e)
}
upper <- function(num_FALSE,num_TRUE,tot_FALSE,tot_TRUE){
  e <- fisher.test(matrix(c(num_FALSE, 
                            num_TRUE,
                            tot_FALSE-num_FALSE, 
                            tot_TRUE-num_TRUE),
                          dimnames = list(c("inactive","active"),c("in ChIP", "not in ChIP")),
                          nrow=2, ncol=2))$conf.int[2]
  return(e)
}
lower <- function(num_FALSE,num_TRUE,tot_FALSE,tot_TRUE){
  e <- fisher.test(matrix(c(num_FALSE, 
                            num_TRUE,
                            tot_FALSE-num_FALSE, 
                            tot_TRUE-num_TRUE),
                          dimnames = list(c("inactive","active"),c("in ChIP", "not in ChIP")),
                          nrow=2, ncol=2))$conf.int[1]
  return(e)
}

active$lower <- apply(active[,c("num_FALSE","num_TRUE","tot_FALSE","tot_TRUE")], MARGIN=1, function(x) lower(x[1],x[2],x[3],x[4]))
active$upper <- apply(active[,c("num_FALSE","num_TRUE","tot_FALSE","tot_TRUE")], MARGIN=1, function(x) upper(x[1],x[2],x[3],x[4]))
active$OR <- apply(active[,c("num_FALSE","num_TRUE","tot_FALSE","tot_TRUE")], MARGIN=1, function(x) estimate(x[1],x[2],x[3],x[4]))
active$p <- apply(active[,c("num_FALSE","num_TRUE","tot_FALSE","tot_TRUE")], MARGIN=1, function(x) pval(x[1],x[2],x[3],x[4]))

active <- active %>% dplyr::mutate(se = sqrt(1/num_TRUE[1]+1/(tot_TRUE - num_TRUE)+1/num_FALSE+1/(tot_FALSE-num_FALSE)))

active %>% dplyr::filter(num_FALSE + num_TRUE > 20,
                         p < 0.05/54)


p<- active %>%
  dplyr::mutate(cl = ifelse(p < 0.05/54, 'sig','nonsig'))%>%
  dplyr::mutate(upper = ifelse(upper > 2, 2, upper)) %>%
  ggplot(.,aes(x=reorder(dark.tfs,OR), y= OR, color = cl)) +
  geom_point()+
  geom_errorbar(aes(ymin = lower, ymax =upper), position = position_dodge(0.9), width = 0)+
  pretty_plot() +
  scale_color_manual(values = c("gray","maroon"))+
  xlab('Dark TFs') +
  ylab("Odds ratio, inactive vs active")+
  ggtitle('"Dark" Transcription Factors')+
  theme(legend.position = 'None')+
  geom_hline(yintercept = 1)+
  coord_flip()
p

cowplot::save_plot(paste0('figures/s4e.pdf'),
                   #paste0("../tables/", pdfname, "_tracks.pdf"),
                   p,
                   base_height = 6,
                   base_width = 6,
                   device = cairo_pdf
)

