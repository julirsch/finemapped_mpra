# Load libraries
library(tidyverse)
library(vroom)
source("code/utils.R")
options(stringsAsFactors = FALSE)

# Read in MPRA data
mpra_df <- vroom("data/preprocess/core_mpra.txt.gz")
mpra_meta_df <- vroom("data/preprocess/mpra_meta.txt.gz")
mpra_df <- mpra_df %>%
  left_join(mpra_meta_df,
            by = c("variant", "cohort")) %>%
  ungroup() 

# Read in and process saturation mutagenesis data
satmut_seq_df <- vroom("data/satmut/satmut_sequences.txt.gz", delim = "\t") %>%
  dplyr::select(final_id, final_lsid, mut_pos, mut_base, sequence) %>%
  dplyr::mutate(final_id = gsub("4763691", "4763491", final_id),
                final_id = gsub("ALT", "Alt", final_id),
                final_lsid = gsub("ALT", "Alt", final_lsid),
                final_id = gsub("6:41924998:C:G:A:wC", "6:41924998:C:G:A:wX", final_id),
                final_id = gsub("6:41924998:C:G:R:wC", "6:41924998:C:G:R:wX", final_id),
                final_lsid = gsub("6:41924998:C:G:A:wC", "6:41924998:C:G:A:wX", final_lsid),
                final_lsid = gsub("6:41924998:C:G:R:wC", "6:41924998:C:G:R:wX", final_lsid)) %>%
  dplyr::mutate(parent_id = gsub(":m.*", "", final_id)) %>%
  distinct(final_id, sequence)
satmut_meta_df <- vroom("data/satmut/satmut_sequences.txt.gz", delim = "\t") %>%
  dplyr::select(final_id, final_lsid, mut_pos, mut_base) %>%
  dplyr::mutate(final_id = gsub("4763691", "4763491", final_id),
                final_id = gsub("ALT", "Alt", final_id),
                final_lsid = gsub("ALT", "Alt", final_lsid),
                final_id = gsub("6:41924998:C:G:A:wC", "6:41924998:C:G:A:wX", final_id),
                final_id = gsub("6:41924998:C:G:R:wC", "6:41924998:C:G:R:wX", final_id),
                final_lsid = gsub("6:41924998:C:G:A:wC", "6:41924998:C:G:A:wX", final_lsid),
                final_lsid = gsub("6:41924998:C:G:R:wC", "6:41924998:C:G:R:wX", final_lsid)) %>%
  dplyr::mutate(parent_id = gsub(":m.*", "", final_id))
satmut_meta2_df <- vroom("data/satmut/satmut_meta.txt", delim = "\t")
satmut_meta2_df <- satmut_meta2_df %>%
  bind_rows(satmut_meta2_df %>%
              dplyr::filter(!grepl("Alt_", final_id)) %>%
              dplyr::mutate(final_id = gsub("R:w", "A:w", final_id),
                            final_lsid = gsub("R:w", "A:w", final_lsid))) %>%
  dplyr::mutate(parent_id = gsub(":m0", "", final_id),
                parent_lsid = gsub(":m0", "", final_lsid))
satmut_meta_df <- satmut_meta_df %>%
  inner_join(satmut_meta2_df %>%
               dplyr::select(-final_id, -final_lsid),
             by = "parent_id")
files <- list.files("data/satmut/", pattern = "*_20230216.out.gz", full.names = T)
celltypes <- c("HEPG2", "K562")
satmut_df <- lapply(1:2, function(x) {
  vroom(files[x], delim = "\t") %>%
    dplyr::mutate(cell_type = celltypes[x])
}) %>%
  bind_rows() %>%
  #dplyr::mutate(parent_id = gsub(":m.*", "", ID)) %>%
  dplyr::filter(project == "UKBB-GTEx_Sat") %>%
  dplyr::filter(!is.na(log2FoldChange)) %>%
  inner_join(satmut_meta_df,
             by = c("ID" = "final_id")) %>%
  #dplyr::filter(is.na(indel)) %>%
  dplyr::select(ID, cell_type, sat_ref_parent, sat_ref, sat_mut, mut_pos, mut_base, indel, DNA_mean, log2FoldChange, lfcSE, padj, var1_position, var2_position, centervar) %>%
  dplyr::mutate(is_haplo = ifelse(grepl("Alt_", sat_ref), T, F)) %>%
  separate(sat_ref, c("var1", "var2"), ":Alt_", remove = F) %>%
  separate(var1, c("var1", "window"), ":w", remove = T) %>%
  dplyr::mutate(window = paste0("w", window),
                var1_tmp = var1,
                var2_tmp = var2,
                var1 = ifelse(centervar == "var1", var1_tmp, var2_tmp),
                var2 = ifelse(centervar == "var1", var2_tmp, var1_tmp)) %>%
  separate(var1, c("chr", "pos", "ref", "alt", "allele"), ":", remove = T) %>%
  separate(var2, c("chr2", "pos2", "ref2", "alt2", "allele2"), ":", remove = T) %>%
  dplyr::mutate(var1 = paste0(chr, ":", pos, ":", ref, ":", alt),
                var2 = paste0(chr2, ":", pos2, ":", ref2, ":", alt2),
                var2 = ifelse(var2 == "NA:NA:NA:NA", NA, var2)) %>%
  dplyr::mutate(mut_pos = ifelse(sat_mut == "m0", var1_position, mut_pos),
                mut_base = ifelse(sat_mut == "m0", ref, mut_base)) %>%
  dplyr::mutate(pos = as.numeric(pos),
                pos2 = as.numeric(pos2),
                pos = pos + (mut_pos - var1_position)) %>%
  dplyr::mutate(allele = ifelse(is.na(allele2), allele, paste0(allele, allele2))) %>%
  dplyr::mutate(sat_ref_parent = ifelse(!is_haplo, sat_ref_parent, paste0(var1, ":", window, ":Alt_", var2))) %>%
  dplyr::mutate(sat_ref_parent = ifelse(!is_haplo, sat_ref_parent, paste0(var1, ":", window, ":Alt_", gsub(":A$|:R$", "", var2)))) %>%
  dplyr::select(ID, cell_type, sat_ref_parent, sat_ref, allele, var1, var2, centervar, window, chr, pos, ref, alt, ref2, alt2, mut_pos, mut_base, is_haplo, "oligomut" = sat_mut, indel, DNA_mean, "log2FC" = log2FoldChange, "log2FC_SE" = lfcSE, padj) %>%
  distinct() %>%
  dplyr::filter(!grepl("4763368.*wC|4763368.*wR|4763491.*wC|4763491.*wR", sat_ref)) %>%
  ungroup() 

# Filter to only elements for this project
keepvars_df <- vroom("data/satmut/satmut_final_128elements.txt", delim = "\t")
satmut_df <- satmut_df %>%
  dplyr::filter(sat_ref_parent %in% keepvars$sat_ref_parent)

# Add emVar annotations
emvars_ct_df <- mpra_df %>% 
  dplyr::filter(cell_type %in% c("K562", "HEPG2"), emVar == T) %>% 
  distinct(variant, cell_type) 
emvars_df <- emvars_ct_df %>% 
  distinct(variant) 
satmut_df <- satmut_df %>%
  dplyr::mutate(var1 = paste0("chr", var1),
                var2 = ifelse(!is.na(var2), paste0("chr", var2), var2)) %>%
  separate(var1, c("remove1", "pos1"), ":", remove = F) %>%
  separate(var2, c("remove2", "pos2"), ":", remove = F) %>%
  dplyr::select(-remove1, -remove2) %>%
  dplyr::mutate(var1_emVar = ifelse(var1 %in% emvars_df$variant, T, F),
                var2_emVar = ifelse(var2 %in% emvars_df$variant, T, F),
                is_var1 = ifelse(pos == pos1, T, F),
                is_var2 = ifelse(pos == pos2, T, F),
                var_emVar = ifelse(var1_emVar == T & is_var1 == T, var1, NA),
                var_emVar = ifelse(var2_emVar == T & is_var2 == T, var2, var_emVar))

# Empirical Bayes
satmut_df <- satmut_df %>%
  group_by(cell_type, sat_ref) %>%
  dplyr::mutate(qntl_lower = quantile(log2FC, 0.25),
                qntl_higher = quantile(log2FC, 0.75),
                prior_mu = mean(log2FC[log2FC > qntl_lower & log2FC < qntl_higher]),
                prior_sigmasq = var(log2FC[log2FC > qntl_lower & log2FC < qntl_higher]),
                qntl_lower = quantile(log2FC, 0.05),
                qntl_higher = quantile(log2FC, 0.95),
                prior_mu_full = mean(log2FC[log2FC > qntl_lower & log2FC < qntl_higher]),
                prior_sigmasq_full = var(log2FC[log2FC > qntl_lower & log2FC < qntl_higher])) %>%
  ungroup() %>%
  dplyr::mutate(obs_mu = log2FC,
                obs_sigmasq = log2FC_SE ^ 2 * 5, 
                post_sigmasq = 1 / (1 / prior_sigmasq + 5 / obs_sigmasq),
                post_mu = post_sigmasq * (prior_mu / prior_sigmasq + (5 * obs_mu) / obs_sigmasq),
                post_sigmasq_full = 1 / (1 / prior_sigmasq_full + 5 / obs_sigmasq),
                post_mu_full = post_sigmasq_full * (prior_mu_full / prior_sigmasq_full + (5 * obs_mu) / obs_sigmasq)) %>%
  dplyr::select(-obs_mu, -obs_sigmasq, -qntl_lower, -qntl_higher)

# Compute log2Skew
satmut_df <- satmut_df %>% 
  dplyr::group_by(cell_type, sat_ref) %>%
  dplyr::mutate(n_obs = length(prior_mu),
                log2FC_baseline = prior_mu[oligomut == "m0"],
                log2FC_SE_baseline = sqrt(prior_sigmasq[oligomut == "m0"] / n_obs),
                log2Skew = log2FC - log2FC_baseline,
                log2Skew_SE = sqrt(log2FC_SE^2 + log2FC_SE_baseline^2),
                post_log2Skew = post_mu_full - log2FC_baseline,
                post_log2FC = post_mu_full) %>%
  dplyr::select(-prior_mu, -prior_sigmasq, -prior_mu_full, -prior_sigmasq_full, -post_sigmasq, -post_mu, -post_sigmasq_full, -post_mu_full, -n_obs) %>%
  ungroup()

# P-values for skew 
satmut_df <- satmut_df %>%
  group_by(sat_ref, cell_type) %>%
  dplyr::mutate(log2Skew_pval = 2 * pnorm(q = abs(log2Skew) / sqrt(log2Skew_SE), lower.tail = FALSE),
                log2Skew_fdr = qvalue::qvalue(log2Skew_pval)$qvalue) %>%
  ungroup()

# Write preprocessed SatMut data
satmut_df %>%
  vroom_write("data/satmut/satmut.txt.gz")




