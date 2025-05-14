## produces haplos_master_table.txt.gz

# Load libraries
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

# Read in haplotype results
#/mnt/sdb/gtex_mpra/final_mpra_data/processed/haplotypes/meta_windows.txt.gz
meta_haplos <- read_delim('data/preprocess/haplos/meta_windows.txt.gz', delim = '\t') %>%
  dplyr::mutate(v2 = gsub("chr:", "chr", v2)) %>%
  dplyr::filter(cell_type != "GM12878")
#/mnt/sdb/gtex_mpra/final_mpra_data/processed/haplotypes/all_windows.txt.gz
haplos <- read_delim('data/preprocess/haplos/all_windows.txt.gz', delim = '\t')

# filter out GM12878s
haplos <- haplos %>%
  separate(v1v2_construct, into = c("v1","v2","center_var","window"), sep = ';', remove=FALSE) %>%
  dplyr::mutate(v2 = gsub("chr:", "chr", v2),
                v1 = gsub("chr:", "chr", v1)) %>%
  dplyr::filter(cell_type != "GM12878")
haplos <- haplos %>%
  left_join(meta_haplos %>%
              dplyr::select(v1, v2, cell_type, refref_log2FC_meta, refalt_log2FC_meta, altref_log2FC_meta, altalt_log2FC_meta, refref_log2FC_meta_SE, refalt_log2FC_meta_SE, altref_log2FC_meta_SE, altalt_log2FC_meta_SE, altref_log2Skew_meta,refalt_log2Skew_meta, altalt_log2Skew_meta, int_log2Skew_meta, any_emVar_meta, int_emVar_meta),
            by = c("v1", "v2", "cell_type"))

# Filter to trait-associated, fine-mapped variants in CREs
haplos_df <- haplos %>% 
  inner_join(mpra %>% 
               dplyr::filter(cohort %in% c("UKBB", "BBJ", "GTEx"), pip > 0.1) %>%
               distinct(variant, CRE) %>%
               dplyr::rename("CRE_v1" = "CRE"),
             by = c("v1" = "variant")) %>%
  inner_join(mpra %>% 
               dplyr::filter(cohort %in% c("UKBB", "BBJ", "GTEx"), pip > 0.1) %>%
               distinct(variant, CRE) %>%
               dplyr::rename("CRE_v2" = "CRE"),
             by = c("v2" = "variant")) %>% 
  dplyr::filter(CRE_v1 > 0 | CRE_v2 > 0) %>%
  distinct()

# Count number of pairs (number in text)
haplos_df %>%
  distinct(v1, v2)

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
haplo_bed_df <- haplos_df %>% 
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

overlap_tmp <- overlap %>%
  left_join(haplos_df %>%
              distinct(v1v2_construct, cell_type, int_log2Skew) %>%
              dplyr::rename("query" = v1v2_construct, 
                            "query_cell_type" = cell_type,
                            "query_int_log2Skew" = int_log2Skew)) %>%
  left_join(haplos_df %>%
              distinct(v1v2_construct, cell_type, int_log2Skew) %>%
              dplyr::rename("subject" = v1v2_construct, 
                            "subject_cell_type" = cell_type,
                            "subject_int_log2Skew" = int_log2Skew)) %>%
  dplyr::filter(query_cell_type == subject_cell_type) %>%
  distinct() %>%
  dplyr::mutate(perc_overlap_bin = cut(perc_overlap, perco_bins, include.lowest = F),
                distance_bin = cut(distance, dist_bins, include.lowest = F))

# Compute empirical correlation  between interaction effects based upon oligo % overlap
overlap_tmp_mat <- overlap_tmp %>%
  dplyr::filter(query != subject,
                distance > 0) %>%
  group_by(perc_overlap_bin, distance_bin) %>%
  dplyr::summarize(cor = cor(query_int_log2Skew, subject_int_log2Skew, use = "pairwise.complete.obs"),
                   distance = round(mean(distance), 1),
                   perc_overlap = round(mean(perc_overlap), 2)) %>%
  na.omit()

# View empirical correlation
overlap_tmp_mat %>%
  tidytext::cast_sparse(distance_bin, perc_overlap_bin, cor)

# Smooth empirical correlations and check result
mod <- mgcv::gam(cor ~ te(distance, perc_overlap), data = overlap_tmp_mat)
newdata <- expand.grid(distance = dist_bins, perc_overlap = perco_bins)
bind_cols(newdata, "pred" = predict(mod, newdata) %>% as.vector()) %>%
  as_tibble() %>% 
  tidytext::cast_sparse(distance, perc_overlap, pred)


# Get smoothed empirical correlation for each pair
overlap_tmp$pred <- predict(mod, overlap_tmp) %>%
  as.vector()

overlap_cor_mat <- overlap_tmp %>%
  tidytext::cast_sparse(query, subject, pred)

heatmap(overlap_cor_mat, Colv = NA, Rowv = NA)

# Loop through and meta-analyze (within cell-types, across windows) using expected empirical covariance matrix
haplos_list <- haplos_df %>%
  arrange(cell_type, v1v2_construct) %>%
  dplyr::select(cell_type, v1, v2, v1v2_construct, int_log2Skew, int_log2SkewSE) %>%
  collapse::rsplit(., ~ cell_type + v1 + v2, flatten = T, keep.by = T)

haplo_meta_out <- lapply(haplos_list, function(x) {
  
  cell_type <- x$cell_type[1]
  v1 <- x$v1[1]
  v2 <- x$v2[1]
  beta <- x$int_log2Skew
  se <- x$int_log2SkewSE
  
  if(dim(x)[1] == 1) {
    pval <- -1 * pnorm(abs(beta / se), lower.tail = F, log.p = T)
    return(c("cell_type" = cell_type, "v1" = v1, "v2" = v2, "beta" = beta, "se" = se, "pval" = pval, "beta2" = beta, "se2" = se, "pval2" = pval, "n_oligo" = 1))
  }
  
  V <- overlap_cor_mat[x$v1v2_construct, x$v1v2_construct]
  diag(V) <- 1
  rownames(V) <- NULL
  colnames(V) <- NULL
  V <- diag(se) %*% V %*% diag(se)
  out <- metafor::rma.mv(beta, V)
  V <- se^2
  out2 <- metafor::rma.mv(beta, V)
  c("cell_type" = cell_type, "v1" = v1, "v2" = v2, 
    "beta" = out$beta, "se" = out$se, "pval" = -1 * log10(out$pval), 
    "beta2" = out2$beta, "se2" = out2$se, "pval2" = -1 * log10(out2$pval),
    "n_oligo" = length(beta))
}) %>%
  bind_rows() %>%
  dplyr::mutate(beta = as.numeric(beta),
                se = as.numeric(se),
                pval = as.numeric(pval),
                beta2 = as.numeric(beta2),
                se2 = as.numeric(se2),
                pval2 = as.numeric(pval2))

# FDR correction
haplo_meta_out$fdr <- -1 * log10(qvalue::qvalue(10^-haplo_meta_out$pval)$qvalue)
haplo_meta_out$fdr2 <- -1 * log10(qvalue::qvalue(10^-haplo_meta_out$pval2)$qvalue)

# Loop through and meta-analyze (across cell-types, across windows) using expected empirical covariance matrix
haplos_list <- haplos_df %>%
  arrange(v1v2_construct) %>%
  dplyr::select(cell_type, v1, v2, v1v2_construct, int_log2Skew, int_log2SkewSE) %>%
  collapse::rsplit(., ~ v1 + v2, flatten = T, keep.by = T)


haplo_meta_ct_out <- lapply(haplos_list, function(x) {
  
  v1 <- x$v1[1]
  v2 <- x$v2[1]
  beta <- x$int_log2Skew
  se <- x$int_log2SkewSE
  
  if(dim(x)[1] == 1) {
    pval <- -1 * pnorm(abs(beta / se), lower.tail = F, log.p = T)
    return(c("v1" = v1, "v2" = v2, "beta" = beta, "se" = se, "pval" = pval, "beta2" = beta, "se2" = se, "pval2" = pval, "n_oligo" = 1))
  }
  
  V <- overlap_cor_mat[x$v1v2_construct, x$v1v2_construct]
  diag(V) <- 1
  rownames(V) <- NULL
  colnames(V) <- NULL
  V <- diag(se) %*% V %*% diag(se)
  out <- metafor::rma.mv(beta, V)
  V <- se^2
  out2 <- metafor::rma.mv(beta, V)
  c("v1" = v1, "v2" = v2, 
    "beta" = out$beta, "se" = out$se, "pval" = -1 * log10(out$pval), 
    "beta2" = out2$beta, "se2" = out2$se, "pval2" = -1 * log10(out2$pval),
    "n_oligo" = length(beta))
}) %>%
  bind_rows() %>%
  dplyr::mutate(beta = as.numeric(beta),
                se = as.numeric(se),
                pval = as.numeric(pval),
                beta2 = as.numeric(beta2),
                se2 = as.numeric(se2),
                pval2 = as.numeric(pval2))

# FDR correction
haplo_meta_ct_out$fdr <- -1 * log10(qvalue::qvalue(10^-haplo_meta_ct_out$pval)$qvalue)
haplo_meta_ct_out$fdr2 <- -1 * log10(qvalue::qvalue(10^-haplo_meta_ct_out$pval2)$qvalue)
haplo_meta_ct_out <- haplo_meta_ct_out %>%
  `colnames<-`(c(colnames(haplo_meta_ct_out)[1:2], paste0("cts_", colnames(haplo_meta_ct_out)[3:11])))


# Combine results 
haplos_new_df <- haplos_df %>% 
  left_join(meta_haplos %>%
              distinct(v1, v2, cell_type, int_log2Skew_meta_nlogp, int_log2Skew_meta_padj)) %>%
  left_join(haplo_meta_out) %>%
  left_join(haplo_meta_ct_out) %>%
  dplyr::mutate(any_active_update = ifelse(refref_active + refalt_active + altref_active + altalt_active > 0, T, F)) %>%
  dplyr::mutate(int_emVar_meta_new = (fdr > -log10(0.05) & (beta > 0.25 | beta < -0.25) & n_oligo > 1) | 
                  (fdr2 > -log10(0.05) & (beta2 > 0.25 | beta2 < -0.25) & n_oligo > 1) | 
                  (cts_fdr > -log10(0.05) & (cts_beta > 0.25 | cts_beta < -0.25) & cts_n_oligo > 1) | 
                  (cts_fdr2 > -log10(0.05) & (cts_beta2 > 0.25 | cts_beta2 < -0.25)) & cts_n_oligo > 1) %>%
  dplyr::mutate(int_emVar_meta_new = ifelse(int_emVar_meta_new + any_active_update == 2, T, F))

## checking interacting emVar number (number in text)
haplos_new_df %>% 
  dplyr::filter(int_emVar_meta_new == T) %>% 
  distinct(v1, v2)

# recode alleles as increasing/decreasing relative to baseline
haplos_new_df <- haplos_new_df %>% 
  # set refref log2skew as 0
  dplyr::mutate(refref_log2Skew_meta = 0) %>%
  # get the lowest log2Skew
  dplyr::mutate(lowest = pmin(altalt_log2Skew_meta, altref_log2Skew_meta, refalt_log2Skew_meta, refref_log2Skew_meta)) %>%
  # denote which measurement was lowest
  dplyr::mutate(lowest_source = case_when(altalt_log2Skew_meta == lowest ~ 'aa',
                                          altref_log2Skew_meta == lowest ~ 'ar',
                                          refalt_log2Skew_meta == lowest ~ 'ra',
                                          refref_log2Skew_meta == lowest ~ 'rr')) %>%
  # shift everything by the lowest
  dplyr::mutate(altalt_log2Skew_recode = altalt_log2Skew_meta - lowest,
                refalt_log2Skew_recode = refalt_log2Skew_meta - lowest,
                altref_log2Skew_recode = altref_log2Skew_meta - lowest,
                refref_log2Skew_recode = refref_log2Skew_meta - lowest) %>%
  # reassign decdec, incdec, decinc, incinc
  dplyr::mutate(log2Skew_v1decv2dec = case_when(lowest_source == 'rr' ~ refref_log2Skew_recode,
                                                lowest_source == 'aa' ~ altalt_log2Skew_recode,
                                                lowest_source == 'ra' ~ refalt_log2Skew_recode,
                                                lowest_source == 'ar' ~ altref_log2Skew_recode),
                log2Skew_v1incv2dec = case_when(lowest_source == 'rr' ~ altref_log2Skew_recode,
                                                lowest_source == 'aa' ~ refalt_log2Skew_recode,
                                                lowest_source == 'ra' ~ altalt_log2Skew_recode,
                                                lowest_source == 'ar' ~ refref_log2Skew_recode),
                log2Skew_v1decv2inc = case_when(lowest_source == 'rr' ~ refalt_log2Skew_recode,
                                                lowest_source == 'aa' ~ altref_log2Skew_recode,
                                                lowest_source == 'ra' ~ refref_log2Skew_recode,
                                                lowest_source == 'ar' ~ altalt_log2Skew_recode),
                log2Skew_v1incv2inc = case_when(lowest_source == 'rr' ~ altalt_log2Skew_recode,
                                                lowest_source == 'aa' ~ refref_log2Skew_recode,
                                                lowest_source == 'ra' ~ altref_log2Skew_recode,
                                                lowest_source == 'ar' ~ refalt_log2Skew_recode)) %>%
  # get rid of intermediate columns
  dplyr::select(-altalt_log2Skew_recode,-refref_log2Skew_recode,-altref_log2Skew_recode,-refalt_log2Skew_recode, -lowest)


## save off
haplos_new_df %>% vroom::vroom_write('data/preprocess/haplos/haplos_master_table.txt.gz')
