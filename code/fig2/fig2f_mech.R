# Load libraries
library(tidyverse)
library(vroom)
library(GenomicRanges)
library(rtracklayer)
library(plyranges)
library(matrixStats)
library(annotables)
library(BuenColors)
library(PNWColors)
library(cowplot)
library(patchwork)
library(ggExtra)
library(ggrastr)
source("code/utils.R")
source("code/theme.R")
options(stringsAsFactors = FALSE)

# Read in processed MPRA results
mpra_df <- vroom("data/preprocess/core_mpra.txt.gz")
mpra_meta_df <- vroom("data/preprocess/mpra_meta.txt.gz")
mpra_df <- mpra_df %>%
  left_join(mpra_meta_df,
            by = c("variant", "cohort")) %>%
  ungroup() 

# Read in dataset for 4 categories of variants
mech_df <- vroom::vroom("/mnt/sdb/gtex_mpra/data/mech_df.txt.gz") 

# Redefine low-PIP categories for consistency
mech_df <- mech_df %>%
  dplyr::filter(!(pip > 0.02 & category %in% c("LowPIP_emVar", "LowPIP_only")))

# Filter to CRE variants
mpra_df <- mpra_df %>%
  dplyr::filter(CRE > 0, cs_id > 0)

# Read in motif results
load("data/motifs/motifs_processed_20230225.RData")
#motifdisrupt_df <- vroom::vroom("/mnt/sdb/gtex_mpra/data/fimo/motifdisrupt.txt.gz")
#motifdisrupt_df <- vroom::vroom("/mnt/sdb/gwas_eqtl_mpra/data/annotations/motif/motifdisrupt.txt.gz")
motifdisrupt_df <- vroom::vroom("data/annotations/motif/motifdisrupt.txt.gz")
#variant_chip_lite_df <- vroom::vroom("/mnt/sdb/gwas_eqtl_mpra/data/annotations/motif/variant_chip_lite.txt.gz")
#motifdisrupt_chip_df <- vroom::vroom("/mnt/sdb/gwas_eqtl_mpra/data/annotations/motif/motifdisrupt_chip_jaspar_hocmoco_only.txt.gz")
motifdisrupt_chip_df <- vroom::vroom("data/annotations/motif/motifdisrupt_chip_jaspar_hocmoco_only.txt.gz")
#motifcorr_filt_df <- vroom::vroom("/mnt/sdb/gtex_mpra/data/motif_corr_filt_jaspar_hocomoco.txt.gz")
#motifcorr_filt_df <- vroom::vroom("/mnt/sdb/gwas_eqtl_mpra/data/annotations/motif/motif_corr_filt_jaspar_hocomoco.txt.gz")
gc()

# Add motif metadata information
motifdisrupt_df <- motifdisrupt_df %>%
  inner_join(motifs_meta %>%
               na.omit() %>% 
               dplyr::select(-motif_db),
             by = "motif") %>%
  dplyr::filter(motif_db %in% c("jaspar", "hocomoco"))

# Filter to CRE variants and make make motif calls
motifdisrupt_df <- motifdisrupt_df %>% 
  dplyr::filter(variant %in% unique(mpra_df$variant)) %>%  
  ungroup() %>% 
  dplyr::mutate(pcta_allele1 = ifelse(pcta_allele1 < 0, 0, pcta_allele1),
                pcta_allele2 = ifelse(pcta_allele2 < 0, 0, pcta_allele2),
                pct = pmax(pct_allele1, pct_allele2),
                pcta = pmax(pcta_allele1, pcta_allele2),
                motif_diff = pct_allele2 - pct_allele1,
                motif_diffa = pcta_allele2 - pcta_allele1,
                pval = pmin(pval),
                widtha = abs(start - stop)) %>%
  dplyr::select(variant, motif_db, motif, TF, width, icscore, pct, pcta, pval, motif_diff, motif_diffa) %>%
  distinct() %>%
  ungroup() %>%
  left_join(motifdisrupt_chip_df,
            by = c("variant", "motif")) %>%
  filter(!is.na(TF)) %>%
  dplyr::mutate(motif_call = ifelse(width >= 8 & pval < 10^-4 & abs(motif_diff) > 0.1 & pct > 0.75 & pcta > 0.75 & icscore > 12, T, F)) %>%
  dplyr::mutate(motif_call_nodiff = ifelse(width >= 8 & pval < 10^-4 & pct > 0.75 & pcta > 0.75 & icscore > 12, T, F)) %>%
  dplyr::mutate(motif_weakcall = ifelse(width >= 8 & pval < 10^-3 & abs(motif_diff) > 0.05 & pct > 0.75 & icscore > 8, T, F)) %>%
  dplyr::mutate(motif_weakcall_nodiff = ifelse(width >= 8 & pval < 10^-3 & pct > 0.75 & icscore > 8, T, F))

# Merge MPRA with motif disruption data
mpra_motif_df <- mpra_df %>% 
  distinct(variant, log2Skew, cell_type, emVar, emVar_any) %>%
  left_join(motifdisrupt_df %>%
              dplyr::filter(motif_db %in% c("jaspar", "hocomoco")) %>%
              distinct(variant, motif, TF, antigen, motif_call, motif_call_nodiff, motif_weakcall, motif_weakcall_nodiff, motif_diff))

# Add motif information to subset of variants for Fig. 2f
mech_motif_df <- mech_df %>% 
  distinct(variant, variant_hg38, category, footprints) %>%
  left_join(motifdisrupt_df %>%
              distinct(variant, motif, TF, antigen, pct, motif_call, motif_call_nodiff, motif_weakcall, motif_weakcall_nodiff, motif_diff),
            by = "variant")

# Read in motif matches for flanking motifs
#motifmatch_df <- vroom::vroom("/mnt/sdb/gtex_mpra/data/fimo/motifmatch.txt.gz")
#motifmatch_chip_df <- vroom::vroom("/mnt/sdb/gtex_mpra/data/motifmatch_chip.txt.gz")
motifmatch_df <- vroom::vroom("data/fimo/motifmatch.txt.gz")
motifmatch_chip_df <- vroom::vroom("data/motifmatch_chip.txt.gz")
motifmatch_df <- motifmatch_df %>%
  inner_join(motifs_meta %>%
               na.omit() %>% 
               dplyr::select(-motif_db),
             by = "motif") %>%
  left_join(motifmatch_chip_df) %>%
  dplyr::mutate(motif_call = ifelse(width >= 8 & pval < 10^-4 & pct > 0.75 & pcta > 0.75 & icscore > 12, T, F)) %>%
  ungroup()

# Read in Enformer
#enf_chip <- vroom::vroom("/mnt/sdb/gtex_mpra/data/annotations/enformer/enformer_CHIP_mech_satmut.txt.gz", delim = "\t", col_names = T) %>%
#  tidytable::as_tidytable()
#enf_dhs <- vroom::vroom("/mnt/sdb/gtex_mpra/data/annotations/enformer/enformer_DHS_mech_satmut.txt.gz", delim = "\t", col_names = T) %>%
#  tidytable::as_tidytable()
enf_chip <- vroom::vroom("data/annotations/enformer/enformer_CHIP_mech_satmut.txt.gz", delim = "\t", col_names = T) %>%
  tidytable::as_tidytable()
enf_dhs <- vroom::vroom("data/annotations/enformer/enformer_DHS_mech_satmut.txt.gz", delim = "\t", col_names = T) %>%
  tidytable::as_tidytable()

# Get SS thresholds (L2 norm threshold for Enformer)
ss_chip <- matrixStats::rowSums2(enf_chip[,-1]^2 %>% as.matrix())
ss_dhs <- matrixStats::rowSums2(enf_dhs[,-1]^2 %>% as.matrix())
ss_df <- bind_cols("variant" = enf_dhs$variant, "ss_chip" = ss_chip, "ss_dhs" = ss_dhs)
ss_thresh <- ss_df %>%
  inner_join(mech_df %>%
               distinct(variant, category)) %>%
  group_by(category) %>%
  dplyr::summarize(ss_chip_thresh = quantile(ss_chip, 0.95),
                   ss_dhs_thresh = quantile(ss_dhs, 0.95))
ss_chip_thresh <- ss_thresh %>% 
  dplyr::filter(category == "LowPIP_only") %>% 
  pull(ss_chip_thresh)
ss_dhs_thresh <- ss_thresh %>% 
  dplyr::filter(category == "LowPIP_only") %>% 
  pull(ss_dhs_thresh)

# Add Enformer calls
mech_motif_df <- mech_motif_df %>%
  left_join(ss_df) %>%
  dplyr::mutate(enf_chip = ifelse(ss_chip > ss_chip_thresh, 1, 0),
                enf_dhs = ifelse(ss_dhs > ss_dhs_thresh, 1, 0))

# Get variants breaking motifs +/- occupied
tmp_motif_call <- mech_motif_df %>%
  dplyr::filter(motif_call == T) %>%
  distinct(variant)
tmp_motif_call_occupied <- mech_motif_df %>%
  dplyr::filter(motif_call == T, !is.na(antigen)) %>%
  distinct(variant)
tmp_motif_weakcall <- mech_motif_df %>%
  dplyr::filter(motif_weakcall == T) %>%
  distinct(variant)
tmp_motif_call_nodiff <- mech_motif_df %>%
  dplyr::filter(motif_call_nodiff == T) %>%
  distinct(variant)
tmp_motif_flank_occupied <- motifmatch_df %>%
  tidytable::filter(motif_call == T, !is.na(antigen)) %>%
  tidytable::group_by(variant, antigen) %>%
  tidytable::summarize(flank_motif2 = max(flank_motif)) %>%
  tidytable::filter(flank_motif2 == 0) %>%
  tidytable::filter(variant %ni% tmp_motif_call_occupied$variant) %>%
  distinct(variant)

# Function for 95% CIs
m_enrich <- function(df, attr, intype) {
  t <- df %>% distinct(variant, category, !! sym(attr)) %>%
    dplyr::count(category, !! sym(attr)) %>% 
    #na.omit() %>%
    group_by(category) %>% 
    dplyr::mutate(tot = sum(n)) %>% 
    dplyr::filter(!! sym(attr) == T) %>% 
    dplyr::select(-!! sym(attr)) %>% 
    dplyr::mutate(prop = n / tot, 
                  lower = binom.confint(n, tot, method = 'wilson')$lower,
                  upper = binom.confint(n, tot, method = 'wilson')$upper,
                  pvalue = binom.test(n, tot)$p.value)
  t %>% 
    dplyr::mutate(attribute = attr, type = intype)
}

# Compute proportions for each category
motif_plot <- bind_rows(m_enrich(mech_motif_df %>% dplyr::mutate(annot = ifelse(variant %in% tmp_motif_call$variant, 1, 0)), "annot", "Breaks motif"),
                        m_enrich(mech_motif_df %>% dplyr::mutate(annot = ifelse(variant %in% tmp_motif_call_occupied$variant, 1, 0)), "annot", "Breaks occupied motif"),
                        m_enrich(mech_motif_df %>% dplyr::mutate(annot = ifelse(variant %in% tmp_motif_flank_occupied$variant, 1, 0)), "annot", "In occupied flanking motif"),
                        m_enrich(mech_motif_df %>% dplyr::mutate(annot = ifelse(footprints > 0, 1, 0)), "annot", "Footprints"),
                        m_enrich(mech_motif_df %>% dplyr::filter(!is.na(enf_chip), !is.na(enf_dhs)) %>% dplyr::mutate(annot = ifelse(enf_chip + enf_dhs > 0, 1, 0)), "annot", "Enformer"),
                        m_enrich(mech_motif_df %>% dplyr::filter(!is.na(enf_chip), !is.na(enf_dhs)) %>% dplyr::mutate(annot = ifelse((enf_chip + enf_dhs > 0) & (variant %ni% tmp_motif_call$variant), 1, 0)), "annot", "Enformer no motif"))
m_enrich(mech_motif_df %>% dplyr::mutate(annot = ifelse(variant %in% tmp_motif_call_nodiff$variant, 1, 0)), "annot", "In occupied motif")
prop.test(c(1039, 2415), c(1530, 14694))$p.value
prop.test(c(892, 2415), c(1742, 14694))$p.value

# Write out high PIP emVar mechanisms
mech_motif_df %>%
  dplyr::filter(motif_call == T, !is.na(antigen)) %>%
  dplyr::distinct(variant, category, motif, TF, pct, motif_diff) %>%
  vroom::vroom_write("tables/stable16.txt.gz")

# Prepare to plot
motif_plot <- motif_plot %>%
  ungroup() %>%
  dplyr::mutate(category = factor(category, levels=c("LowPIP_only", "HighPIP_only", "LowPIP_emVar", "HighPIP_emVar")),
                type = factor(type, levels=c("Breaks motif", "Breaks occupied motif", "Footprints", "In occupied flanking motif", "Enformer", "Enformer no motif")))

# Plot proportion in each category (Fig. 2f)
p1 <- ggplot(motif_plot, aes(x = type, y = prop, fill = category)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                size = 0.5,
                width = 0,
                col = "black",
                position = position_dodge(0.9, preserve = "single")) +
  scale_color_manual(values = c(BuenColors::jdb_palette('brewer_spectra')[8],BuenColors::jdb_palette('brewer_spectra')[2],
                                BuenColors::jdb_palette('brewer_spectra')[1], 
                                '#40366F'),
                     aesthetics = c("color","fill")) +
  pretty_plot()  + theme(text = element_text(size = 20),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank())+
  theme( plot.background = element_blank() ,
         panel.grid.major = element_blank() ,
         panel.grid.minor = element_blank() ,
         panel.border = element_blank() ,
         panel.background = element_blank(),
         legend.position = 'None') +
  theme(axis.line = element_line(color = 'black'))
p1 

# Save plot (Fig 2f)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(2))
cowplot::save_plot(
  paste0(outdir, "mpra_motif_overlap.pdf"),
  plt_combined,
  base_height = 2,
  base_width = 3,
  device = cairo_pdf
)

# # Number of TF / motifs disrupted by high PIP emVars
# mech_motif_df %>%
#   dplyr::filter(category == "HighPIP_emVar", motif_call == T) %>%
#   distinct(TF) %>%
#   dplyr::summarize(n = length(TF))
# 
# # TF motif proportion plot
# keep_tfs <- plot_in3 %>%
#   dplyr::filter(fdr < 0.05) %>%
#   distinct(motif, TF)
# 
# b <- mech_motifs_df %>%
#   dplyr::filter(category == "HighPIP_emVar", motif_call == T, !is.na(antigen)) %>%
#   inner_join(keep_tfs) %>%
#   distinct(variant, TF) %>%
#   dplyr::count(variant) %>%
#   dplyr::rename(n_vars = n)
# a <- mech_motifs_df %>%
#   dplyr::filter(category == "HighPIP_emVar", motif_call == T, !is.na(antigen)) %>%
#   inner_join(keep_tfs) %>%
#   distinct(variant, TF) %>%
#   dplyr::count(TF) %>%
#   dplyr::mutate(n = n / sum(n))
# 
# p1 <- mech_motifs_df %>%
#   dplyr::filter(category == "HighPIP_emVar", motif_call == T, !is.na(antigen)) %>%
#   inner_join(keep_tfs) %>%
#   distinct(variant, TF) %>%
#   left_join(a,
#             by = "TF") %>%
#   left_join(b,
#             by = "variant") %>%
#   arrange(-n) %>%
#   distinct(TF, n) %>%
#   arrange(-n) %>%
#   print(n = 200)
# ggplot(aes(x = reorder(TF, -n), alpha = cut(n_vars, c(0, 1, 5, 10)))) +
#   geom_bar(aes(y = ..count../sum(..count..)), fill = PNWColors::pnw_palette('Sunset2')[1]) +
#   #pretty_plot() +
#   my_theme +
#   theme(legend.position = "none", axis.text.x = element_text(angle=90, hjust=1))
# p1
# 
# plt_combined <- p1 + plot_layout(nrow = 1, heights = c(2))
# cowplot::save_plot(
#   paste0(outdir, "mpra_motif_proportion.pdf"),
#   plt_combined,
#   base_height = 2,
#   base_width = 4,
#   device = cairo_pdf
# )

  


