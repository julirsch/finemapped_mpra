

# Load libraries
source("code/utils.R")
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

# Read in CRE-only haplotypes
haplos_df <- vroom('data/preprocess/haplos/haplos_master_table.txt.gz')

# Define and filter to best observation for each haplotype
haplos_df <- haplos_df %>%
  dplyr::group_by(v1,v2,window,center_var, cell_type)%>%
  dplyr::mutate(best_library = ifelse((sqrt(mean_Plasmid_altalt) + sqrt(mean_Plasmid_refref)+sqrt(mean_Plasmid_altref) + sqrt(mean_Plasmid_refalt)) ==
                                        max(sqrt(mean_Plasmid_altalt) + sqrt(mean_Plasmid_refref)+sqrt(mean_Plasmid_altref) + sqrt(mean_Plasmid_refalt)),1,0)) %>%
  dplyr::filter(best_library == 1) %>%
  dplyr::ungroup()


# Wide to long dataframe
haplos_long_df <- haplos_df %>%
  distinct(v1, v2, center_var, window, cell_type, 
           altalt_log2Skew, altref_log2Skew, refalt_log2Skew,
           altalt_Log2FC, altref_Log2FC,refalt_Log2FC) %>%
  pivot_longer(cols = c('altalt_log2Skew','altref_log2Skew','refalt_log2Skew',
                        'altalt_Log2FC', 'altref_Log2FC','refalt_Log2FC'),
               names_to = 'name',
               values_to = 'value') %>%
  tidyr::separate(col = 'name',
                  into = c('haplo','measurement'),
                  sep = '_') 


# Count number of pairs (number in text)
haplos_df %>%
  distinct(v1, v2)

# get percent overlaps of each window
# Get ID to merge
#/mnt/sdb/gtex_mpra/mpra_design/haplotypes_info.txt.gz
haplo_info_df <- vroom::vroom("/mnt/sdb/gwas_eqtl_mpra/data/design/haplotypes_info.txt.gz") %>%
  separate(name, c("tmp1", "center_variant", "variant1_allele", "variant2_allele", "window"), "_", remove = F) %>%
  dplyr::mutate(haplotype = paste0(variant1_allele, "_", variant2_allele)) %>%
  dplyr::mutate(v1v2_construct = paste0(v1, ";", v2, ";", center_variant, ";", window)) %>%
  distinct()

# Read in bed files of haplotypes
#/mnt/sdb/gtex_mpra/results/encode_beds/
files <- list.files("/mnt/sdb/gwas_eqtl_mpra/data/annotations/encode_beds/", pattern = "*.bed.gz", full.names = T)
haplo_beds <- vroom::vroom(files, col_names = F) %>%
  dplyr::filter(grepl("Alt", X4))

# Munge and convert to Granges
haplo_bed_df <- haplos_df %>% distinct(v1,v2,center_var, window, cell_type, v1v2_construct) %>%  
  left_join(haplo_info_df %>%
              distinct(ID, v1v2_construct, distance)) %>% 
  left_join(haplo_beds %>%
              dplyr::select("ID" = X4, "chr" = X1, "start" = X2, "end" = X3)) %>%
  distinct(v1, v2, v1v2_construct, distance, chr, start, end) %>%
  na.omit()
haplo_bed_gr <- haplo_bed_df %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chr",
                           start.field =  "start",
                           end.field = "end",
                           keep.extra.columns = T)

# Get % overlap, match, and create bins
hits <- findOverlaps(haplo_bed_gr, haplo_bed_gr, minoverlap = 1)
overlap <- pintersect(haplo_bed_gr[queryHits(hits)],
                      haplo_bed_gr[subjectHits(hits)]) %>% as_tibble() %>%
  bind_cols("subject" = haplo_bed_gr[subjectHits(hits)]$v1v2_construct,
            "query" = haplo_bed_gr[queryHits(hits)]$v1v2_construct,
            "subject_width" = haplo_bed_gr[subjectHits(hits)]@ranges@width - 1,
            "query_width" = haplo_bed_gr[queryHits(hits)]@ranges@width - 1) %>%
  dplyr::mutate(perc_overlap = (width - 1) / pmin(subject_width, query_width)) %>%
  distinct(subject, query, perc_overlap, distance)

dist_bins <- c(0, 5, 10, 15, 20, 25, 50, 75, 100, 200)
perco_bins <- c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.95, 1)

haplos_long_df <- haplos_long_df %>%
  left_join(haplos_df %>% 
              distinct(v1,v2, center_var, window,v1v2_construct)) 
# for skew
overlap_tmp <- overlap %>%
  left_join(haplos_long_df %>%
              dplyr::select(-center_var, -window) %>% 
              dplyr::rename("query" = v1v2_construct, 
                            "query_cell_type" = cell_type,
                            "query_haplo" = haplo,
                            "query_measurement" = measurement,
                            "query_value" = value),
            by = c('query')) %>%
  left_join(haplos_long_df %>%
              dplyr::select(-center_var, -window) %>% 
              dplyr::rename("subject" = v1v2_construct, 
                            "subject_cell_type" = cell_type,
                            "subject_haplo" = haplo,
                            "subject_measurement" = measurement,
                            "subject_value" = value),
            by = c('subject', 'v1','v2')) %>%
  dplyr::filter(query_cell_type == subject_cell_type,
                query_measurement == subject_measurement,
                query_haplo == subject_haplo,
                query != subject) %>% 
  distinct() %>%
  dplyr::mutate(perc_overlap_bin = cut(perc_overlap, perco_bins, include.lowest = F))


skew <- lapply(30:100, function(x){
  overlap_tmp %>%
    dplyr::filter(query != subject,
                  query_measurement == 'log2Skew' & subject_measurement == 'log2Skew',
                  distance > 0,
                  perc_overlap >= x/100-2/100,
                  perc_overlap <= x/100+2/100) %>%
    dplyr::summarize(cor = cor(query_value, subject_value, use = "pairwise.complete.obs"),
                     perc_overlap = round(mean(perc_overlap), 2),
                     n = n_distinct(v1,v2)
                     ) %>%
    na.omit() %>%
    dplyr::mutate(perc_bin = x)
}) %>%
  bind_rows() 


activity <- lapply(30:100, function(x){
  overlap_tmp %>%
    dplyr::filter(query != subject,
                  query_measurement == 'Log2FC' & subject_measurement == 'Log2FC',
                  distance > 0,
                  perc_overlap >= x/100-2/100,
                  perc_overlap <= x/100+2/100) %>%
    dplyr::summarize(cor = cor(query_value, subject_value, use = "pairwise.complete.obs"),
                     perc_overlap = round(mean(perc_overlap), 2),
                     n = n_distinct(v1,v2)
                     ) %>%
    na.omit() %>%
    dplyr::mutate(perc_bin = x)
}) %>%
  bind_rows()



combo <- full_join(activity, skew, by = c('perc_overlap','perc_bin','n'), suffix = c('.activity','.skew')) %>%
pivot_longer(cols = c('cor.activity','cor.skew'),
             names_to = 'type',
             values_to = 'cor',
             names_prefix = 'cor.')
  
  
p <- ggplot(combo, #combo[rep(row.names(combo), times = combo$n), ],
            aes(x=perc_overlap*100, y = cor, group = type, color = type, fill = type))+
  geom_point(aes(size=sqrt(n)),alpha = 0.5)+
  geom_smooth(method = "gam",
              aes(weight = n/mean(n)),
              se = TRUE,
              formula = y ~ s(x, bs = "cs", fx = TRUE, k = 5))+
  scale_colour_manual(values = c('#cc5c76','#f9ad2a'))+
  scale_fill_manual(values = c('#cc5c76','#f9ad2a'))+
  pretty_plot() +
  xlab('Percent overlap of windows')+
  ylab('Smoothed Empirical Correlations')+
  ylim(0,1)+
  scale_size_continuous(range = c(0,3))

p


cowplot::save_plot(
  paste0('figures/fig3/s12g'),
  p,
  base_height = 4,
  base_width = 5,
  device = cairo_pdf)  


#############################

# check to make sure that the only vairaitns we're dropping are those represented by only 1 window
tmp <-  haplos_df %>% dplyr::filter(v1v2_construct %ni% overlap_tmp$subject | v1v2_construct %ni% overlap_tmp$query) %>% distinct(v1,v2, v1v2_construct, center_var, window, cell_type) 
tmp %>% group_by(v1,v2) %>% dplyr::filter(n_distinct(v1v2_construct)>1)

