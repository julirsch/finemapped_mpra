# Load libraries
library(tidyverse)
library(vroom)
library(BuenColors)
library(binom)
library(cowplot)
library(patchwork)
library(rtracklayer)
source("code/utils.R")
source("code/github_ready/prc.functions.R")
options(stringsAsFactors = FALSE)

# Read in processed MPRA results
mpra_df <- vroom("data/preprocess/core_mpra.txt.gz") %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0))

# Read in ChIP-seq data
mpra_chip_sums <- vroom('/mnt/sdb/gwas_eqtl_mpra/data/preprocess/chip_mpra.txt.gz') 
# binarize chip data
tmp <- mpra_chip_sums %>% dplyr::select(starts_with('TF')) 
tmp[,1:1009] <- (tmp[, 1:1009] > 0) + 0
mpra_chip_sums['chip_sums2'] = rowSums(tmp[ , 1:1009], na.rm=TRUE)
mpra_chip_sums <- mpra_chip_sums %>% dplyr::select(variant, chip_sums, chip_sums2)

# Read in E2G
files <- list.files("/mnt/sdb/gtex_mpra/data/e2g/", ".bed.gz", full.names = T)[-c(45, 78)]
e2g_gr <- vroom::vroom(files) %>%
  distinct("chr" = `#chr`, start, end, "gene" = TargetGene, "ENSGID" = TargetGeneEnsemblID, "cell_type" = CellType, isSelfPromoter, "rE2G" = Score, "TSS_dist" = distanceToTSS.Feature, "ABC" = ABCScoreDNaseOnlyAvgHicTrack2.Feature) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chr", 
                           start.field =  "start", 
                           end.field = "end",
                           keep.extra.columns = T)

# Liftover MPRA
variant_gr <- mpra_df %>%
  dplyr::filter(cohort != 'control') %>% 
  distinct(variant) %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  na.omit() %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "position", 
                           end.field = "position",
                           keep.extra.columns = T)
ch <- import.chain("/mnt/sdb/ukbb-finemapping/data/liftover/hg19ToHg38.over.chain")
variant_gr <- liftOver(variant_gr, ch)
variant_gr <- unlist(variant_gr)
variant_gr <- variant_gr[!(duplicated(variant_gr$variant) | duplicated(variant_gr$variant, fromLast = TRUE)),]

# Find overlaps
idx <- findOverlaps(variant_gr, e2g_gr)
e2g_filt_df <- bind_cols(variant_gr[idx@from] %>% 
                           as_tibble() %>%
                           dplyr::select(variant),
                         e2g_gr[idx@to] %>%
                           as_tibble()) %>%
  distinct()

# Add additional annotation columns
eqtl_mpra_df <- mpra_df %>%
  dplyr::filter(cohort == 'GTEx' | grepl('tissue',type)) %>%
  add_annots("CRE_yes_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 1, CRE == 0), .) %>%
  add_annots("CRE_no_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 0, CRE == 1), .) %>%
  add_annots("CRE_yes_active_any_yes", rlang::exprs(active_any == 1, active_any == 0, CRE == 1, CRE == 0), .) %>%
  left_join(.,mpra_chip_sums, by='variant') %>%
  dplyr::mutate(chip_thresh = ifelse(chip_sums > 10, 1, 0)) %>%
  dplyr::mutate(chip2_thresh = ifelse(chip_sums2 > 10, 1, 0)) %>%
  dplyr::mutate(e2g = ifelse(variant %in% unique(e2g_filt_df$variant), T, F)) %>%
  add_annots("e2g_yes_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, e2g == 1, e2g== 0), .) %>%
  add_annots("e2g_no_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, e2g == 0, e2g == 1), .) 

traits_mpra_df <- mpra_df %>% 
  dplyr::filter(cohort %in% c('UKBB',"BBJ") | type %in% c('loc_CS','loc_PIP10','annot_CS','annot_PIP10','null_CS','null_PIP10')) %>% 
  add_annots("CRE_yes_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 1, CRE == 0), .) %>%
  add_annots("CRE_no_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 0, CRE == 1), .) %>%
  add_annots("CRE_yes_active_any_yes", rlang::exprs(active_any == 1, active_any == 0, CRE == 1, CRE == 0), .) %>%
  left_join(.,mpra_chip_sums, by='variant') %>% 
  dplyr::mutate(chip_thresh = ifelse(chip_sums > 10, 1, 0)) %>%
  dplyr::mutate(chip2_thresh = ifelse(chip_sums2 > 10, 1, 0)) %>%
  dplyr::mutate(e2g = ifelse(variant %in% unique(e2g_filt_df$variant), T, F)) %>%
  add_annots("e2g_yes_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, e2g == 1, e2g== 0), .) %>%
  add_annots("e2g_no_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, e2g == 0, e2g == 1), .) 

# Make df of causal variants for PrC
eqtl_prset_df <- make_prset(eqtl_mpra_df, "eQTLs")
traits_prset_df <-  make_prset(traits_mpra_df, "complex_traits")

# Count number of causal variants in data frame
bind_rows(eqtl_prset_df %>% 
            dplyr::filter(causal == T) %>% 
            dplyr::select(variant),
          traits_prset_df %>% 
            dplyr::filter(causal == T) %>% 
            dplyr::select(variant)) %>%
  distinct(variant)

# Calculate precision and recall
annotlist <- c("emVar_any", "CRE", 
               "active_any","CRE_yes_active_any_yes",
               "promoter","chip_thresh","chip2_thresh",
               "CRE_yes_emVar_any_yes", "CRE_no_emVar_any_yes",
               "e2g", "e2g_yes_emVar_any_yes", "e2g_no_emVar_any_yes")
eqtl_out <- lapply(annotlist, function(x) {prec_rec(x, eqtl_prset_df)}) %>%
  bind_rows() %>%
  dplyr::mutate(cohort = "eQTL", .before = annot)
traits_out <- lapply(annotlist, function(x) {prec_rec(x, traits_prset_df)}) %>%
  bind_rows() %>%
  dplyr::mutate(cohort = "Complex traits", .before = annot)
out <- bind_rows(eqtl_out, traits_out)

out$annot = factor(out$annot, levels = c('active_any',"chip_thresh","chip2_thresh", "CRE","CRE_yes_emVar_any_yes","emVar_any","CRE_no_emVar_any_yes","promoter","CRE_yes_active_any_yes","e2g", "e2g_yes_emVar_any_yes", "e2g_no_emVar_any_yes"))

out %>%
  vroom::vroom_write(paste0('/mnt/sdb/gwas_eqtl_mpra/reviews/tables/', Sys.Date(), '.updated.prc.expanded.txt'))

# Plot extended PrC
pr <- ggplot(data = out, aes(x = recall, y = precision, color=annot, shape = cohort)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin = prec_lower,ymax = prec_upper),width=0)+
  geom_errorbarh(aes(xmin = rec_lower,xmax = rec_upper),height=0)+ 
  BuenColors::pretty_plot(fontsize = 20)+
  theme(aspect.ratio=1,
        panel.border = element_blank(),
        axis.line = element_line())+ 
  scale_color_manual(values = BuenColors::jdb_palette("corona")[c(4,6,11,2,3,1,5,7,10,12,9,16)]) +
  geom_hline(yintercept = 0.5, color="grey50", linetype='dashed')


pr
# Calculate sensitivity and specificity
eqtl_out <- lapply(annotlist, function(x) {sens_spec(x, eqtl_prset_df)}) %>%
  bind_rows() %>%
  dplyr::mutate(cohort = "eQTL", .before = annot)
traits_out <- lapply(annotlist, function(x) {sens_spec(x, traits_prset_df)}) %>%
  bind_rows() %>%
  dplyr::mutate(cohort = "Complex traits", .before = annot)
out <- bind_rows(eqtl_out, traits_out)

out$annot = factor(out$annot, levels = c('active_any',"chip_thresh","chip2_thresh", "CRE","CRE_yes_emVar_any_yes","emVar_any","CRE_no_emVar_any_yes","promoter","CRE_yes_active_any_yes","e2g", "e2g_yes_emVar_any_yes", "e2g_no_emVar_any_yes"))

# Write out table (part of supp table 10)
out %>%
  vroom::vroom_write(paste0('tables/stabl10.updated.ss.expanded.txt'))

# plot sensitivity and specificity (Supplementary Fig 9)
ss <- ggplot(data = out, aes(x = sensitivity, y = specificity, color=annot, shape = cohort)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin = spec_lower,ymax = spec_upper),width=0)+
  geom_errorbarh(aes(xmin = sens_lower,xmax = sens_upper),height=0)+ 
  BuenColors::pretty_plot(fontsize = 20)+
  theme(aspect.ratio=1,
        panel.border = element_blank(),
        axis.line = element_line())+ 
  scale_color_manual(values = BuenColors::jdb_palette("corona")[c(4,6,11,2,3,1,5,7,10,12,9,16)]) +
  geom_hline(yintercept = 0.5, color="grey50", linetype='dashed')



ss <- ss +
  theme(legend.position = "None")

leg <- get_legend(pr)
pr <- pr +
  theme(legend.position = "None")

# Convert to a ggplot and print
as_ggplot(leg)

# Save plot (Supplementary Fig 9)
plt_combined <- ggarrange(pr, ss, ncol = 2, labels = c("a","b"))

plt_combined
cowplot::save_plot(paste0('figures/s9ab.pdf'),
                   plt_combined,
                   base_height = 6,
                   base_width = 12,
                   device = cairo_pdf
)

cowplot::save_plot(paste0('figures/s9ab.legend.pdf'),
                   leg,
                   base_height = 12,
                   base_width = 12,
                   device = cairo_pdf
)
