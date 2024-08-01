# Load libraries
library(tidyverse)
library(vroom)
source("code/utils.R")
options(stringsAsFactors = FALSE)

# Read in processed MPRA results
mpra_df <- vroom("data/preprocess/core_mpra.txt.gz")

# Meta-analyze MPRA
mpra_meta_df <- mpra_df %>%
  ungroup() %>%
  dplyr::select(variant, cohort, cell_type, A_log2FC, A_log2FC_SE, log2FC, log2FC_SE, log2Skew, Skew_SE) %>%
  distinct() %>%
  na.omit() %>%
  #group_by(variant)
  dplyr::group_by(variant, cohort) %>%
  dplyr::mutate(log2FC_A_meta = sum(A_log2FC * (1 / A_log2FC_SE^2)) / sum(1 / A_log2FC_SE^2),
                log2FC_meta = sum(log2FC * (1 / log2FC_SE^2)) / sum(1 / log2FC_SE^2),
                log2FC_meta_SE = sqrt(1 / sum(1 / log2FC_SE^2)),
                log2Skew_meta =  sum(log2Skew * (1 / Skew_SE^2)) / sum(1 / Skew_SE^2),
                log2Skew_meta_SE = sqrt(1 / sum(1 / Skew_SE^2)),
                log2FC_meta_nlogp = -1 * pnorm(abs(log2FC_meta / log2FC_meta_SE), lower.tail = F, log.p = T),
                log2Skew_meta_nlogp = -1 * pnorm(abs(log2Skew_meta / log2Skew_meta_SE), lower.tail = F, log.p = T)) %>%
  ungroup() %>%
  distinct(variant, cohort, log2FC_A_meta, log2FC_meta, log2FC_meta_SE, log2Skew_meta, log2Skew_meta_SE, log2FC_meta_nlogp, log2Skew_meta_nlogp) %>%
  dplyr::mutate(lof2FC_meta_padj = -1 * log10(p.adjust(10^(-1 * log2FC_meta_nlogp), "bonferroni")),
                lof2Skew_meta_padj = -1 * log10(p.adjust(10^(-1 * log2Skew_meta_nlogp), "fdr"))) %>%
  dplyr::mutate(active_meta = ifelse(lof2FC_meta_padj >= -log10(0.01) & abs(log2FC_meta) >= 1, T, F),
                emVar_meta = ifelse(active_meta & lof2Skew_meta_padj >= -log10(0.1) & !is.na(lof2Skew_meta_padj) & abs(log2Skew_meta) >= 0, T, F))

# Write out
mpra_meta_df %>%
  vroom_write("data/preprocess/mpra_meta.txt.gz")


                emVar_meta = ifelse(active_meta & lof2Skew_meta_padj >= -log10(0.1) & !is.na(lof2Skew_meta_padj) & abs(log2Skew_meta) >= 0, T, F))