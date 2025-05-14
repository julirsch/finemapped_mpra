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
code_path<-""

source(paste0(code_path,"code/Preprocess/utils.R"))

# Load motifs
load(paste0(code_path,"data/preprocess/motifs_processed_20230225.RData"))
motifs_pfml <- c(hocomoco_pfml, jaspar_pfml, cisbp_pfml, vierstra_pfml)

# Read in MPRA data
mpra <- vroom::vroom(paste0(code_path,"data/preprocess/core_mpra.txt.gz")) %>%
  distinct() %>% 
  dplyr::mutate(cohort = ifelse(is.na(cohort), "control", cohort))

#Files
final_sequences<-paste0(code_path,"data/preprocess/final.sequences.txt")
satmut_metadata<-paste0(code_path,"data/SatMut_MPRA/satmut_metadata_v4.txt")
satmut_final_128<-paste0(code_path,"data/SatMut_MPRA/satmut_final_128elements.txt")

satmut_rc_gb<-paste0(code_path,"data/preprocess/satmut_20230413_with_gaussian_blocks_v8.txt.gz")
satmut_rc_motif<-paste0(code_path,"data/preprocess/satmut_motif_calls_v8.txt.gz")
satmut_rc_topmatch<-paste0(code_path,"data/preprocess/block_top_matches_v8_50.txt.gz")
satmut_rc_submatch<-paste0(code_path,"data/preprocess/sub_block_top_matches_v8_4814.txt.gz")
satmut_rc_subsubmatch<-paste0(code_path,"data/preprocess/sub_sub_block_top_matches_v8_4814.txt.gz")

motifbeta_CRE_file<-paste0(code_path,"data/preprocess/motifbeta_CRE_120k_3cts.txt.gz")

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
    # legend.position = "none",
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
#files <-  files[!grepl("emVAR", files)]
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



# Read in and munge gaussian activity blocks
satmut2 <- read_delim(satmut_rc_gb, delim = "\t") %>%
  dplyr::mutate(chr = ifelse(is.na(chr), "X", chr)) %>%
  dplyr::select(ID, sat_ref_parent, sat_ref, cell_type, mut_pos, gaussian_block, block_direction) %>%
  left_join(satmut %>%
              distinct(ID, var_emVar))

#Removed RST
#blocks <- read_delim(satmut_rc_motif, delim = "\t") %>%
#  dplyr::filter(match_pearson > 0.75) %>%
#  mutate(motif_name = gsub(" -$", "", motif_name),
#         motif_name = gsub(" \\+$", "", motif_name),
#         motif = gsub(".*_", "", motif_name)) %>%   
#  inner_join(motifs_meta,
#             by = "motif") %>%
#  dplyr::select(-motif_name)

blocks2 <- read_delim(satmut_rc_topmatch, delim = "\t") %>%
  dplyr::mutate(motif_length = abs(motif_stop - motif_start)) %>%
  dplyr::filter(combined_pearson > 0.75) %>%
  mutate(motif_name = gsub(" -$", "", motif_name),
         motif_name = gsub(" \\+$", "", motif_name),
         motif = gsub(".*_", "", motif_name)) %>%   
  inner_join(motifs_meta,
             by = "motif") %>%
  dplyr::select(-motif_name)

blocks_sub2 <- read_delim(satmut_rc_submatch, delim = "\t") %>%
  dplyr::mutate(motif_length = abs(motif_stop - motif_start)) %>%
  dplyr::filter(combined_pearson > 0.75) %>%
  mutate(motif_name = gsub(" -$", "", motif_name),
         motif_name = gsub(" \\+$", "", motif_name),
         motif = gsub(".*_", "", motif_name)) %>%   
  inner_join(motifs_meta,
             by = "motif") %>%
  dplyr::select(-motif_name) %>%
  dplyr::select(-sub_block_id, -sub_block_start, -sub_block_stop)

blocks_subsub2 <- read_delim(satmut_rc_subsubmatch, delim = "\t") %>%
  dplyr::mutate(motif_length = abs(motif_stop - motif_start)) %>%
  dplyr::filter(combined_pearson > 0.75) %>%
  mutate(motif_name = gsub(" -$", "", motif_name),
         motif_name = gsub(" \\+$", "", motif_name),
         motif = gsub(".*_", "", motif_name)) %>%   
  inner_join(motifs_meta,
             by = "motif") %>%
  dplyr::select(-motif_name) %>%
  dplyr::select(-sub_block_id, -sub_block_start, -sub_block_stop, -sub_sub_block_id, -sub_sub_block_start, -sub_sub_block_stop)

# Add GB motif matches to SatMut
blocks2_gr <- bind_rows(blocks2, blocks_sub2, blocks_subsub2) %>%
  dplyr::mutate(cell_type = ifelse(cell_type == "HepG2", "HEPG2", cell_type)) %>%
  dplyr::rename("sat_ref" = sequence_id) %>%
  makeGRangesFromDataFrame(.,
                           seqnames.field = "sat_ref", 
                           start.field = "motif_start", 
                           end.field = "motif_stop",
                           keep.extra.columns = T)

satmut2_gr <- satmut2 %>%
  makeGRangesFromDataFrame(.,
                           seqnames.field = "sat_ref", 
                           start.field = "mut_pos", 
                           end.field = "mut_pos",
                           keep.extra.columns = T)

idx <- findOverlaps(satmut2_gr, blocks2_gr)
satmut2_motif_df <- bind_cols(satmut2_gr[idx@from] %>%
                                as_tibble() %>%
                                dplyr::select(ID, sat_ref_parent, "sat_ref" = seqnames, cell_type, var_emVar, mut_pos = "start", gaussian_block, block_direction),
                              blocks2_gr[idx@to] %>%
                                as_tibble(.name_repait = "minimal") %>%
                                dplyr::rename("cell_type2" = cell_type)) %>%
  dplyr::filter(cell_type == cell_type2) %>%
  dplyr::select(-width, -strand, -cell_type2) %>%
  distinct()

satmut_tmp <- satmut2 %>%
  dplyr::select(ID, cell_type, mut_pos, gaussian_block, block_direction) %>%
  dplyr::mutate(in_gb = ifelse(gaussian_block > 0, 1, 0)) %>%
  left_join(satmut2_motif_df %>%
              dplyr::select(ID, cell_type, mut_pos) %>%
              distinct() %>%
              dplyr::mutate(in_gbmotif = 1),
            by = c("ID", "cell_type", "mut_pos")) %>%
  dplyr::mutate(in_gbmotif = ifelse(is.na(in_gbmotif), 0, 1)) %>%
  dplyr::mutate(in_gbandmotif = ifelse(in_gb + in_gbmotif > 0, 1, 0))

satmut <- satmut %>%
  left_join(satmut_tmp,
            by = c("ID", "cell_type", "mut_pos"))

# Proportion GB activating repressing
gb_dir_out <- satmut %>%
  dplyr::filter(gaussian_block != 0) %>%
  dplyr::distinct(sat_ref, cell_type, gaussian_block, block_direction) %>%
  dplyr::count(block_direction) %>%
  dplyr::mutate(perc = n / sum(n))

# Read in motif betas from OLS
beta_df <- vroom::vroom(motifbeta_CRE_file)

# Motifs in each block
satmut_motif_enrich_df <- satmut2_motif_df %>%
  dplyr::select(sat_ref, cell_type, gaussian_block, block_direction, motif) %>%
  filter(gaussian_block > 0) %>%
  distinct(sat_ref, cell_type, gaussian_block, block_direction, motif) %>%
  count(cell_type, block_direction, motif, name = "num") %>% 
  dplyr::left_join(gb_dir_out %>%
                     dplyr::select(block_direction, "denom" = n)) %>%
  mutate(block_direction2 = ifelse(block_direction == -1, "R", "A")) %>%
  pivot_wider(id_cols = c("cell_type", "motif"), names_from = c("block_direction2"), values_from = c("num", "denom")) %>%
  collapse::replace_na() %>%
  mutate(denom_A = max(denom_A),
         denom_R = max(denom_R)) %>%
  group_by(cell_type, motif) %>% 
  mutate(num_A = num_A + 1,
         num_R = num_R + 1,
         denom_A = denom_A + 1,
         denom_R = denom_R + 1,
         prop_A = (num_A / denom_A),
         prop_R = (num_R / denom_R),
         enrich = (num_A / denom_A) / (num_R / denom_R),
         enrich_lower = fisher.test(matrix(c(num_A, denom_A, num_R, denom_R), nrow = 2))$conf.int[1],
         enrich_upper = fisher.test(matrix(c(num_A, denom_A, num_R, denom_R), nrow = 2))$conf.int[2],
         pval = fisher.test(matrix(c(num_A, denom_A, num_R, denom_R), nrow = 2))$p.value)

# Merge with motif betas
satmut_motif_enrich_tmp_df <- satmut_motif_enrich_df %>% 
  ungroup() %>% 
  dplyr::mutate(motif = gsub(".*_", "", motif)) %>% 
  inner_join(beta_df %>%
               dplyr::filter(cell_type %in% c("K562_adj", "HEPG2_adj")) %>%
               dplyr::mutate(cell_type = gsub("_adj", "", cell_type)),
             by = c("motif" = "motif", "cell_type")) %>%
  dplyr::filter(prop_elem > 0.001) %>%
  left_join(motifs_meta, by = "motif")

# Plot final figure
set.seed(1234)
p1 <- satmut_motif_enrich_tmp_df %>%
  arrange(-beta) %>%
  ggplot(aes(x = log2(enrich), y = beta, color = cell_type, size = -log10(pval))) +
  geom_point() +
  scale_color_manual(values = ct_colors) + 
  scale_size(range = c(0, 2)) +
  pretty_plot() +
  my_theme +
  ylab("Activity contribution (120k CREs)") +
  xlab("Motif enrichment (active vs repressive AB)") +
  coord_cartesian(y = c(-0.6, 1.4), x = c(-6, 6)) +
  ggpubr::stat_cor(aes(group = "1"), method = "spearman")

# Save the plot
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(2))
cowplot::save_plot(
  paste0(outdir, "fig4/fig4c_satmut_abmotifs_v_beta.pdf"),
  plt_combined,
  base_height = 2,
  base_width = 2.2,
  device = cairo_pdf
)

