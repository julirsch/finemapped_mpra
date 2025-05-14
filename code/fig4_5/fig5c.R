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

source(paste0(code_path,"code/Preprocess/utils.R"))

# Load motifs
load(paste0(code_path,"data/preprocess/motifs_processed_20230225.RData"))
motifs_pfml <- c(hocomoco_pfml, jaspar_pfml, cisbp_pfml, vierstra_pfml)

# Read in MPRA data
mpra <- vroom::vroom(paste0(code_path,"data/preprocess/core_mpra.txt.gz")) %>%
  distinct() %>% 
  dplyr::mutate(cohort = ifelse(is.na(cohort), "control", cohort))

# Files
final_sequences <- paste0(code_path,"data/preprocess/final.sequences.txt")
satmut_metadata <- paste0(code_path,"data/SatMut_MPRA/satmut_metadata_v4.txt")
satmut_final_128 <- paste0(code_path,"data/SatMut_MPRA/satmut_final_128elements.txt")

hg19_to_hg38_chain_file<-paste0(code_path,"annotations/hg19ToHg38.over.chain")
mammalian241_file<-paste0(code_path,"annotations/241-mammalian-2020v2.bigWig")


# Cell type colors
ct_colors <- c(
  pnw_palette('Bay',8,type='continuous')[5],
  pnw_palette('Bay',8,type='continuous')[3],
  pnw_palette('Starfish')[5],
  pnw_palette('Bay',8,type='continuous')[2],
  pnw_palette('Bay',8,type='continuous')[8]
)
names(ct_colors) <- c("A549", "HCT116", "HEPG2", "K562", "SKNSH")

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
files <- list.files(paste0(code_path,"data/SatMut_MPRA"), pattern = "*_20230216.out", full.names = T)
files <-  files[!grepl("emVAR", files)]
celltypes <- c("HEPG2", "K562")
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

# Add emVar annoations
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

# Look at concordance with original MPRA allelic skews
satmut_emVar_var <- satmut %>% 
  group_by(sat_ref_parent, cell_type) %>% 
  dplyr::filter(is_var1 == T, var1_emVar == T, allele %in% c("A", "R", is_haplo == T)) %>%
  distinct(var_emVar, cell_type, allele, log2FC_baseline) %>%
  dplyr::mutate(log2Skew_A = log2FC_baseline[allele == "A"],
                log2Skew_R = log2FC_baseline[allele == "R"],
                log2Skew_AR = log2Skew_A - log2Skew_R) %>%
  dplyr::filter(allele == "R") %>%
  distinct(var_emVar, cell_type, log2Skew_AR)
satmut_emVar_var1 <- satmut %>% 
  group_by(sat_ref_parent, cell_type) %>% 
  dplyr::filter(is_var1 == T, var1_emVar == T, allele %in% c("AR", "RR", is_haplo == T)) %>%
  distinct(var_emVar, cell_type, allele, log2FC_baseline) %>%
  dplyr::mutate(log2Skew_A = log2FC_baseline[allele == "AR"],
                log2Skew_R = log2FC_baseline[allele == "RR"],
                log2Skew_AR = log2Skew_A - log2Skew_R) %>%
  dplyr::filter(allele == "RR") %>%
  distinct(var_emVar, cell_type, log2Skew_AR)
satmut_emVar_var2 <- satmut %>% 
  group_by(sat_ref_parent, cell_type) %>% 
  dplyr::filter(is_var2 == T, var2_emVar == T, allele %in% c("RA", "RR", is_haplo == T)) %>%
  distinct(var_emVar, cell_type, allele, log2FC_baseline) %>%
  dplyr::mutate(log2Skew_A = log2FC_baseline[allele == "RA"],
                log2Skew_R = log2FC_baseline[allele == "RR"],
                log2Skew_AR = log2Skew_A - log2Skew_R) %>%
  dplyr::filter(allele == "RR") %>%
  distinct(var_emVar, cell_type, log2Skew_AR)
satmut_emVar_var_tmp <- bind_rows(satmut_emVar_var, satmut_emVar_var1, satmut_emVar_var2)

# Compute "ground truth" using collapsed baselines for each set of oligos
satmut_emVar_var <- satmut %>% 
  group_by(sat_ref_parent, cell_type) %>% 
  dplyr::filter(is_var1 == T, var1_emVar == T, mut_base == alt, allele == "R") %>% 
  distinct(ID, var_emVar, cell_type, log2Skew, post_log2Skew)
satmut_emVar_var1 <- satmut %>% 
  group_by(sat_ref_parent, cell_type) %>% 
  dplyr::filter(is_var1 == T, var1_emVar == T, mut_base == alt, allele == "RR") %>% 
  distinct(ID, var_emVar, cell_type, log2Skew, post_log2Skew)
satmut_emVar_var2 <- satmut %>% 
  group_by(sat_ref_parent, cell_type) %>% 
  dplyr::filter(is_var2 == T, var2_emVar == T, mut_base == alt2, allele == "RR") %>% 
  distinct(ID, var_emVar, cell_type, log2Skew, post_log2Skew)
satmut_emVar_var_tmp2 <- bind_rows(satmut_emVar_var, satmut_emVar_var1, satmut_emVar_var2)

# Merge in original MPRA skews
mpra_skew <- mpra %>%
  dplyr::filter(cohort %ni% "control") %>%
  dplyr::mutate(filt = sqrt(mean_Plasmid_ref) + sqrt(mean_Plasmid_alt)) %>% 
  group_by(variant, cell_type) %>%
  dplyr::filter(filt == max(filt)) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>% 
  distinct(variant, cell_type, "og_log2Skew" = log2Skew, "og_log2Skew_SE" = Skew_SE)
mpra_plot1 <- satmut_emVar_var_tmp %>%
  left_join(satmut_emVar_var_tmp2) %>%
  ungroup() %>% 
  inner_join(mpra_skew,
             by = c("var_emVar" = "variant", "cell_type")) %>%
  left_join(emvars_ct %>%
              dplyr::mutate(emVar_ctmatch = T),
            by = c("var_emVar" = "variant", "cell_type"))

# Match with phyloP
in_gr <- satmut %>%
  ungroup() %>%
  distinct(chr, pos, ID) %>%
  na.omit() %>%
  dplyr::mutate(chr = paste0("chr", chr)) %>%
  makeGRangesFromDataFrame(.,
                           seqnames.field = "chr", 
                           start.field = "pos", 
                           end.field = "pos",
                           keep.extra.columns = T)

ch <- import.chain(hg19_to_hg38_chain_file)
in_gr <- liftOver(in_gr, ch)
in_gr <- unlist(in_gr)
in_gr <- in_gr[!(duplicated(in_gr$ID) | duplicated(in_gr$ID, fromLast = TRUE)),]

phylop_gr <- rtracklayer::import(mammalian241_file, format = "BigWig", selection = in_gr)
in_df <- in_gr %>%
  as_tibble() %>%
  distinct(ID, seqnames, start) %>%
  inner_join(phylop_df <- phylop_gr %>%
               as_tibble() %>%
               distinct(seqnames, start, score),
             by = c("seqnames", "start"))

satmut <- satmut %>%
  left_join(in_df %>%
              dplyr::select(ID, "phylop_241m" = score),
            by = "ID") %>%
  dplyr::mutate(conserved = ifelse(phylop_241m > 2.27, T, F)) %>%
  dplyr::mutate(phylop_241m_trunc = ifelse(phylop_241m < 0, 0, phylop_241m))

# Pick a single best element by highest m0 log2FC
mpra_plot <- satmut %>%
  group_by(sat_ref_parent) %>%
  dplyr::filter(log2FC_baseline == max(log2FC_baseline)) %>%
  dplyr::filter(oligomut != "m0") %>%
  ungroup()

# Per substitution
mpra_plot_sub <- mpra_plot %>%
  group_by(sat_ref_parent, pos, mut_base) %>%
  dplyr::mutate(post_log2Skew = max(abs(post_log2Skew), na.rm = T)) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup()

# Investigate emVar effects vs max effect
mpra_plot_subvar <- mpra_plot_sub %>%
  group_by(sat_ref_parent) %>%
  dplyr::mutate(rank = rank(-post_log2Skew)) %>%
  dplyr::filter(rank <= 5) %>%
  dplyr::summarize(max_effect = mean(1 - 2^(-1 * abs(post_log2Skew))),
                   log2FC_baseline = mean(log2FC_baseline)) %>%
  dplyr::select(sat_ref_parent, max_effect, log2FC_baseline) %>%
  inner_join(mpra_plot1 %>%
               dplyr::filter(emVar_ctmatch == T) %>%
               dplyr::rename("variant" = var_emVar) %>%
               group_by(sat_ref_parent, variant) %>%
               dplyr::filter(abs(post_log2Skew) == max(abs(post_log2Skew))) %>%
               dplyr::mutate(var_effect = 1 - 2^(-1 * abs(post_log2Skew))) %>%
               distinct(variant, sat_ref_parent, var_effect),
             by = "sat_ref_parent") %>%
  left_join(mpra_plot_sub %>%
              group_by(sat_ref_parent) %>%
              dplyr::summarize(per_conserved = sum(conserved) / length(conserved)),
            by = "sat_ref_parent")

# Add cohort labels
mpra_plot_subvar <- mpra_plot_subvar %>%
  ungroup() %>%
  left_join(mpra %>%
              dplyr::filter(cohort %in% c("UKBB", "BBJ"), pip > 0.1) %>%
              distinct(variant) %>%
              dplyr::mutate(in_UKBB = T),
            by = "variant") %>%
  dplyr::mutate(in_UKBB = ifelse(is.na(in_UKBB), F, T)) %>%
  left_join(mpra %>%
              dplyr::filter(cohort %in% c("GTEx"), pip > 0.1) %>%
              distinct(variant) %>%
              dplyr::mutate(in_GTEx = T),
            by = "variant")  %>%
  dplyr::mutate(in_GTEx = ifelse(is.na(in_GTEx), F, T)) %>%
  dplyr::mutate(cohort = case_when(in_GTEx == T & in_UKBB == T ~ "UKBB/GTEx",
                                   in_UKBB == T ~ "UKBB",
                                   in_GTEx == T ~ "GTEx",
                                   T ~ "Neither"))

mpra_plot_subvar <- mpra_plot_subvar %>% 
  dplyr::mutate(var_effect = ifelse(var_effect < max_effect, var_effect, max_effect))

# Create the plot
p1 <- mpra_plot_subvar %>% 
  arrange(cohort) %>%
  ggplot(aes(x = -1 * max_effect, y = -1 * var_effect, color = cohort, size = log2FC_baseline)) +
  geom_abline(aes(slope = 1, intercept = 0), color = "grey") +
  geom_point() +
  scale_size(range = c(0, 2)) +
  xlab("Largest single nucleotide effect") +
  ylab("Trait-associated variant effect") +
  pretty_plot() +
  my_theme +
  xlim(-0.96, -0.4) +
  ylim(-0.96, 0)

plt_combined <- p1 + plot_layout(nrow = 1, heights = c(2))
cowplot::save_plot(
  paste0(outdir, "fig5/fig5c_satmut_var_es.pdf"),
  plt_combined,
  base_height = 2,
  base_width = 2.2,
  device = cairo_pdf
)

