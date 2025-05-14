# Load libraries
library(tidyverse)
library(GenomicRanges)
library(BuenColors)
library(rtracklayer)
library(reshape2)
library(binom)
library(magrittr)
library(MESS)
library(RcppRoll)
library(corrplot)
library(patchwork)
library(MotifDb)
library(PNWColors)
library(vroom)

# Set parameters
outdir <- "figures/"
code_path <- ""

# Load utility functions
source(paste0(code_path, "code/Preprocess/utils.R"))

# Load motifs
load(paste0(code_path, "data/preprocess/motifs_processed_20230225.RData"))
motifs_pfml <- c(hocomoco_pfml, jaspar_pfml, cisbp_pfml, vierstra_pfml)

# Read in MPRA data
mpra <- vroom::vroom(paste0(code_path, "data/preprocess/core_mpra.txt.gz")) %>%
  distinct() %>% 
  dplyr::mutate(cohort = ifelse(is.na(cohort), "control", cohort))

# Define file paths
final_sequences <- paste0(code_path, "data/preprocess/final.sequences.txt")
satmut_metadata <- paste0(code_path, "data/SatMut_MPRA/satmut_metadata_v4.txt")
satmut_final_128 <- paste0(code_path, "data/SatMut_MPRA/satmut_final_128elements.txt")

# Setup theme
my_theme <-
  BuenColors::pretty_plot(fontsize = 8) +
  BuenColors::L_border() +
  theme(
    plot.background = element_blank(),
    plot.margin = margin(0, 0.1, 0, 0.1, unit = "cm"),
    plot.title = element_text(hjust = 4e-3, margin = margin(b = -12)),
    legend.position = c(0.3, 1),
    legend.justification = c(1, 1),
    legend.title = element_text(margin = margin(0, 0, 0, 0)),
    legend.background = element_blank(),
    legend.key.size = unit(0.2, "cm")
  )

# Meta-analyze MPRA
mpra_meta <- mpra %>%
  ungroup() %>%
  dplyr::select(variant, cohort, cell_type, A_log2FC, A_log2FC_SE, log2FC, log2FC_SE, log2Skew, Skew_SE) %>%
  distinct() %>%
  na.omit() %>%
  dplyr::group_by(variant, cohort) %>%
  dplyr::mutate(log2FC_A_meta = sum(A_log2FC * (1 / A_log2FC_SE^2)) / sum(1 / A_log2FC_SE^2),
                log2FC_meta = sum(log2FC * (1 / log2FC_SE^2)) / sum(1 / log2FC_SE^2),
                log2FC_meta_SE = sqrt(1 / sum(1 / log2FC_SE^2)),
                log2Skew_meta = sum(log2Skew * (1 / Skew_SE^2)) / sum(1 / Skew_SE^2),
                log2Skew_meta_SE = sqrt(1 / sum(1 / Skew_SE^2)),
                log2FC_meta_nlogp = -1 * pnorm(abs(log2FC_meta / log2FC_meta_SE), lower.tail = F, log.p = T),
                log2Skew_meta_nlogp = -1 * pnorm(abs(log2Skew_meta / log2Skew_meta_SE), lower.tail = F, log.p = T)) %>%
  ungroup() %>%
  distinct(variant, cohort, log2FC_A_meta, log2FC_meta, log2FC_meta_SE, log2Skew_meta, log2Skew_meta_SE, log2FC_meta_nlogp, log2Skew_meta_nlogp) %>%
  dplyr::mutate(lof2FC_meta_padj = -1 * log10(p.adjust(10^(-1 * log2FC_meta_nlogp), "bonferroni")),
                lof2Skew_meta_padj = -1 * log10(p.adjust(10^(-1 * log2Skew_meta_nlogp), "fdr"))) %>%
  dplyr::mutate(active_meta = ifelse(lof2FC_meta_padj >= -log10(0.01) & abs(log2FC_meta) >= 1, T, F),
                emVar_meta = ifelse(active_meta & lof2Skew_meta_padj >= -log10(0.1) & !is.na(lof2Skew_meta_padj) & abs(log2Skew_meta) >= 0, T, F))

mpra <- mpra %>%
  left_join(mpra_meta,
            by = c("variant", "cohort")) %>%
  ungroup() 

# Read in saturation mutagenesis data
satmut_seq <- read_delim(final_sequences, delim = "\t") %>%
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

satmut_meta <- read_delim(final_sequences, delim = "\t") %>%
  dplyr::select(final_id, final_lsid, mut_pos, mut_base) %>%
  dplyr::mutate(final_id = gsub("4763691", "4763491", final_id),
                final_id = gsub("ALT", "Alt", final_id),
                final_lsid = gsub("ALT", "Alt", final_lsid),
                final_id = gsub("6:41924998:C:G:A:wC", "6:41924998:C:G:A:wX", final_id),
                final_id = gsub("6:41924998:C:G:R:wC", "6:41924998:C:G:R:wX", final_id),
                final_lsid = gsub("6:41924998:C:G:A:wC", "6:41924998:C:G:A:wX", final_lsid),
                final_lsid = gsub("6:41924998:C:G:R:wC", "6:41924998:C:G:R:wX", final_lsid)) %>%
  dplyr::mutate(parent_id = gsub(":m.*", "", final_id))

satmut_meta2 <- read_delim(satmut_metadata, delim = "\t")
satmut_meta2 <- satmut_meta2 %>%
  bind_rows(satmut_meta2 %>%
              dplyr::filter(!grepl("Alt_", final_id)) %>%
              dplyr::mutate(final_id = gsub("R:w", "A:w", final_id),
                            final_lsid = gsub("R:w", "A:w", final_lsid))) %>%
  dplyr::mutate(parent_id = gsub(":m0", "", final_id),
                parent_lsid = gsub(":m0", "", final_lsid))

satmut_meta <- satmut_meta %>%
  inner_join(satmut_meta2 %>%
               dplyr::select(-final_id, -final_lsid),
             by = "parent_id")

# Find SatMut output files
files <- list.files(paste0(code_path,"data/SatMut_MPRA"), pattern = "*_20230216.out", full.names = T)
files <- files[!grepl("emVAR", files)]
celltypes <- c("HEPG2", "K562")

# Create satmut dataframe
satmut <- lapply(1:2, function(x) {
  read_delim(files[x], delim = "\t") %>%
    dplyr::mutate(cell_type = celltypes[x])
}) %>%
  bind_rows() %>%
  dplyr::filter(project == "UKBB-GTEx_Sat") %>%
  dplyr::filter(!is.na(log2FoldChange)) %>%
  inner_join(satmut_meta,
             by = c("ID" = "final_id")) %>%
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
keepvars <- vroom::vroom(satmut_final_128, delim = "\t")
satmut <- satmut %>%
  dplyr::filter(sat_ref_parent %in% keepvars$sat_ref_parent)

# Add emVar annotations
emvars_ct <- mpra %>% 
  dplyr::filter(cell_type %in% c("K562", "HEPG2"), emVar == T) %>% 
  distinct(variant, cell_type) 
emvars <- emvars_ct %>% 
  distinct(variant) 
satmut <- satmut %>%
  dplyr::mutate(var1 = paste0("chr", var1),
                var2 = ifelse(!is.na(var2), paste0("chr", var2), var2)) %>%
  separate(var1, c("remove1", "pos1"), ":", remove = F) %>%
  separate(var2, c("remove2", "pos2"), ":", remove = F) %>%
  dplyr::select(-remove1, -remove2) %>%
  dplyr::mutate(var1_emVar = ifelse(var1 %in% emvars$variant, T, F),
                var2_emVar = ifelse(var2 %in% emvars$variant, T, F),
                is_var1 = ifelse(pos == pos1, T, F),
                is_var2 = ifelse(pos == pos2, T, F),
                var_emVar = ifelse(var1_emVar == T & is_var1 == T, var1, NA),
                var_emVar = ifelse(var2_emVar == T & is_var2 == T, var2, var_emVar))

# Empirical Bayes
satmut <- satmut %>%
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
satmut <- satmut %>% 
  group_by(cell_type, sat_ref) %>%
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
satmut <- satmut %>%
  group_by(sat_ref, cell_type) %>%
  dplyr::mutate(log2Skew_pval = 2 * pnorm(q = abs(log2Skew) / sqrt(log2Skew_SE), lower.tail = FALSE),
                log2Skew_fdr = qvalue::qvalue(log2Skew_pval)$qvalue) %>%
  ungroup()

# Now generate the satmut_nuc_sub_effects.pdf plot

# Pick a single best element by highest m0 log2FC
mpra_plot <- satmut %>%
  group_by(sat_ref_parent) %>%
  dplyr::filter(log2FC_baseline == max(log2FC_baseline)) %>%
  dplyr::filter(oligomut != "m0") %>%
  ungroup()

# Per nucleotide
mpra_plot_col <- mpra_plot %>%
  group_by(sat_ref_parent, pos) %>%
  dplyr::mutate(post_log2Skew = max(abs(post_log2Skew), na.rm = T)) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() 

# Per substitution
mpra_plot_sub <- mpra_plot %>%
  group_by(sat_ref_parent, pos, mut_base) %>%
  dplyr::mutate(post_log2Skew = max(abs(post_log2Skew), na.rm = T)) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup()

# Distribution of effects across elements
mpra_plot_col_in <- mpra_plot_col %>%
  group_by(sat_ref_parent) %>%
  dplyr::summarize(effect = quantile(2^(-1 * abs(post_log2Skew)), 1 - seq(0, 0.99, 0.01)),
                   percentile = round((1 - seq(0, 0.99, 0.01)) * 100)) %>%
  group_by(percentile) %>%
  dplyr::summarize(mean = mean(effect),
                   lower = quantile(effect, 0.1),
                   upper = quantile(effect, 0.9)) %>%
  dplyr::mutate(type = "position")

mpra_plot_sub_in <- mpra_plot_sub %>%
  group_by(sat_ref_parent) %>%
  dplyr::summarize(effect = quantile(2^(-1 * abs(post_log2Skew)), 1 - seq(0, 0.99, 0.01)),
                   percentile = round((1 - seq(0, 0.99, 0.01)) * 100)) %>%
  group_by(percentile) %>%
  dplyr::summarize(mean = mean(effect),
                   lower = quantile(effect, 0.1),
                   upper = quantile(effect, 0.9)) %>%
  dplyr::mutate(type = "substitution")

# Create the plot
p1 <- bind_rows(mpra_plot_col_in, mpra_plot_sub_in) %>%
  dplyr::mutate(percentile = 100 - as.numeric(percentile)) %>%
  ggplot(aes(x = percentile, y = mean * 100, color = type, fill = type)) +
  geom_ribbon(aes(ymin = lower * 100, ymax = upper * 100), alpha = 0.2, lwd = 0) +
  geom_path() +
  ylab("% change in activity") +
  xlab("Percentile") +
  scale_color_manual(values = PNWColors::pnw_palette(name="Starfish",n=7,type="discrete")[c(6,2)]) +
  scale_fill_manual(values = PNWColors::pnw_palette(name="Starfish",n=7,type="discrete")[c(6,2)]) +
  pretty_plot() +
  my_theme

# Save the plot
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(2))
cowplot::save_plot(
  paste0(outdir, "fig5/fig5b_satmut_nuc_sub_effects.pdf"),
  plt_combined,
  base_height = 2,
  base_width = 2.2,
  device = cairo_pdf
)

