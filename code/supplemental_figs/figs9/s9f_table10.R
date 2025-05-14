
# Load libraries
library(tidyverse)
library(vroom)
library(BuenColors)
library(binom)
library(cowplot)
library(patchwork)
library(PNWColors)
library(ggpubr)
library(rtracklayer)
source("code/utils.R")
source("code/github_ready/prc.functions.R")
options(stringsAsFactors = FALSE)


# Read in processed MPRA results
mpra_df <- vroom("data/preprocess/core_mpra.txt.gz") %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0)) %>%
  dplyr::filter(cohort == 'UKBB')

# Get positive set variants, don't filter down yet to a single row
positives_df <- mpra_df %>%
  dplyr::filter(pip > 0.9, 
                pchisq(chisq_marginal, 1, log.p = TRUE, lower.tail = F) / -log(10) > -log10(5 * 10^-8), 
                type %in% c("CS", "PIP10"),
                consequence %ni% c("synonymous", "missense", "LoF")) %>% 
  dplyr::group_by(variant) %>%
  dplyr::mutate(active_any = ifelse(any(active==T),T, F)) %>%
  dplyr::mutate(pip = max(pip, na.rm = T)) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  dplyr::mutate(Skew_logPadj = ifelse(is.na(Skew_logPadj), 0, Skew_logPadj)) %>%
  arrange(-active,-Skew_logPadj) %>% 
  ungroup() %>%
  dplyr::mutate(causal = T) %>%
  distinct(variant, cell_type, active, emVar, causal)

# Get negative set variants
negatives_df <- mpra_df %>%
  dplyr::filter(!is.na(pip)) %>%
  group_by(variant) %>%
  dplyr::mutate(active_any = ifelse(any(active==T),T, F)) %>%
  dplyr::mutate(pip = max(pip, na.rm = T)) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  dplyr::mutate(Skew_logPadj = ifelse(is.na(Skew_logPadj), 0, Skew_logPadj)) %>%
  arrange(-active,-Skew_logPadj) %>% 
  ungroup() %>%
  dplyr::filter(pip < 0.01, 
                type == "CS", 
                consequence %ni% c("synonymous", "missense", "LoF"), 
                variant %ni% positives_df$variant) %>% 
  ungroup() %>%
  dplyr::mutate(causal = F) %>%
  distinct(variant, cell_type, active, emVar, causal)

# Bind together
mpra_mc <- bind_rows(positives_df, negatives_df)

# Let's look at cell type specific correlations
DHS_Meuleman.meta <- read_delim("/mnt/sdb/gwas_eqtl_mpra/data/annotations/meuleman/DHS_Index_and_Vocabulary_metadata.tsv", delim = "\t", col_names = T)
DHS_Meuleman.meta <- DHS_Meuleman.meta[-734,]

# Read in DHS quantitative
in1.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/meuleman/meuleman_en_gr.rds")
in2.gr <- in1.gr[,0]

# Convert variants to granges
variant.gr <- mpra_mc %>%
  distinct(variant) %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  na.omit() %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "position", 
                           end.field = "position",
                           keep.extra.columns = T)

# Overlap variants with DHS
idx <- findOverlaps(variant.gr, in2.gr)
dhs_vars <- bind_cols(variant.gr[idx@from] %>%
                        as_tibble() %>%
                        dplyr::select(variant),
                      in1.gr[idx@to] %>%
                        as_tibble(.name_repait = "minimal") %>%
                        dplyr::select(-seqnames, -start, -end, -width, -strand))

# Altius columns corresponding to tissues of interest
altius <- DHS_Meuleman.meta %>%
  dplyr::rename("tissue" = `Biosample name`,
                "celltype" = `Altius Aggregation ID`,
                "system" = System) %>%
  dplyr::select(tissue, celltype) %>%
  dplyr::mutate(tissue = ifelse(tissue == "HepG2", "HEPG2", tissue)) %>% 
  dplyr::filter(tissue %in% unique(mpra_df$cell_type))  %>%
  .$celltype

# Create starting df
indf <- mpra_mc %>%
  distinct(variant,cell_type) %>%
  left_join(dhs_vars,
            by = "variant") %>%
  dplyr::select(variant, cell_type, any_of(altius)) %>%
  distinct() 


# Get matching tissues
out <- indf %>%
  ungroup() %>%
  pivot_longer(cols = !c(variant, cell_type)) %>%
  left_join(DHS_Meuleman.meta %>%
              dplyr::rename("tissue" = `Biosample name`,
                            "celltype" = `Altius Aggregation ID`,
                            "system" = System) %>%
              dplyr::select(tissue, celltype) %>%
              dplyr::mutate(tissue = ifelse(tissue == "HepG2", "HEPG2", tissue)),
            by = c("name" = "celltype")) %>%
  ungroup() %>%
  dplyr::filter(tissue %in% cell_type) %>%
  dplyr::group_by(variant, cell_type, tissue) %>%
  dplyr::filter(value == max(value)) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup()

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
variant_gr <- mpra_mc %>%
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

# Keep for MPRA tissues
e2g_filt_cts_df <- e2g_filt_df %>%
  dplyr::filter(cell_type %in% c("K562", "HepG2", "A549", "SK-N-DZ", "SK-N-MC")) %>%
  dplyr::mutate(cell_type = ifelse(cell_type == "HepG2", "HEPG2", cell_type),
                cell_type = ifelse(cell_type %in% c("SK-N-DZ", "SK-N-MC"), "SKNSH", cell_type))

# add in cell-type matched e2g and dhs
mpra_overall <- mpra_mc %>% 
  left_join(out %>%
              dplyr::filter(cell_type == tissue, value > 0) %>%
              distinct(variant, cell_type, value),
            by = c('variant','cell_type')) %>%
  replace(is.na(.),0) %>%
  dplyr::rename(is.CRE = value) %>%
  dplyr::mutate(is.CRE = ifelse(is.CRE > 0, 1,0)) %>% 
  left_join(e2g_filt_cts_df %>%
              distinct(variant, cell_type) %>%
              dplyr::mutate(is.e2g = 1),
            by = c('variant','cell_type')) %>%
  replace(is.na(.),0) 

mpra_overall

prset_df <- mpra_overall %>%
  add_annots("CRE_yes_emVar_yes", rlang::exprs(emVar == 1, emVar == 0, is.CRE == 1, is.CRE == 0), .) %>%
  add_annots("e2g_yes_emVar_yes", rlang::exprs(emVar == 1, emVar == 0, is.e2g == 1, is.e2g == 0), .) 

# cell-type matching PrC function
prec_rec_ct <- function(in_annot, df) {
  
  # Get 2x2 table
  out_df <- df %>%
    dplyr::select(!! sym(in_annot), causal, cell_type) %>%
    dplyr::mutate(annot = case_when(!! sym(in_annot) > 0 ~ 1,
                                    !! sym(in_annot) == 0 ~ 0)) %>%
    dplyr::filter(!is.na(!! sym(in_annot))) %>%
    dplyr::count(cell_type, causal, annot)
  
  # Compute precision, recall, and 95% CIs
  out_df <- out_df %>%
    group_by(cell_type) %>% 
    dplyr::summarize(TP = sum(n[causal == T & annot == T]),
                     FP = sum(n[causal == F & annot == T]),
                     FN = sum(n[causal == T & annot == F]),
                     TN = sum(n[causal == F & annot == F]),
                     total_pos = sum(n[annot == T]),
                     total = sum(n)) 
  
  # there are more negatives:
  out_df <- out_df  %>%
    dplyr::mutate(FPR = FP / (FP + TN),
                  N = TP + FN,
                  FP_new = FPR * N,
                  TN_new = N - FP_new) %>%
    dplyr::mutate(FP = FP_new,
                  TN = TN_new) 
  # end new
  out_df <- out_df %>%
    dplyr::mutate(precision = TP / (TP + FP),
                  recall = TP / (TP + FN),
                  prec_upper = binom.confint(TP, (TP + FP), method = "wilson")$upper,
                  prec_lower = binom.confint(TP, (TP + FP), method = "wilson")$lower,
                  rec_upper = binom.confint(TP, (TP + FN), method = "wilson")$upper,
                  rec_lower = binom.confint(TP, (TP + FN), method = "wilson")$lower) %>% 
    dplyr::mutate(annot = in_annot, .before = TP)
  return(out_df)
}

# get cell type matched precision and recall
annotlist <- c("emVar", "is.CRE", "is.e2g", "CRE_yes_emVar_yes", "e2g_yes_emVar_yes")
traits_out <- lapply(annotlist, function(x) {prec_rec_ct(x, prset_df)}) %>%
  bind_rows() %>%
  dplyr::mutate(cohort = "Complex Traits", .before = annot)

# plot
clrs = c(
  pnw_palette('Bay',8,type='continuous')[5],
  #pnw_palette('Bay',8,type='continuous')[3],
  pnw_palette('Starfish')[5],
  pnw_palette('Bay',8,type='continuous')[1],
  pnw_palette('Bay',8,type='continuous')[8],
  'black'
)


names(clrs) = c("A549", "HEPG2", "K562","SKNSH","any")


p1 <- ggplot(data = traits_out %>% dplyr::filter(cohort == 'Complex Traits',
                                                 annot %in% c("is.e2g","e2g_yes_emVar_yes",
                                                              "is.CRE", "CRE_yes_emVar_yes")), aes(x = recall, y = precision, color=cell_type, shape = annot)) +
  scale_shape_manual(values=c(1, 5, 19, 18))+
  geom_point(size=5) +
  xlim(0,0.15)+
  geom_errorbar(aes(ymin = prec_lower,ymax = prec_upper),width=0)+
  geom_errorbarh(aes(xmin = rec_lower,xmax = rec_upper),height=0)+ 
  BuenColors::pretty_plot(fontsize = 20)+
  theme(#legend.position = 'none', 
    aspect.ratio=1,
    panel.border = element_blank(),
    axis.line = element_line())+ 
  scale_color_manual(values = clrs) +
  geom_hline(yintercept = 0.5, color="grey50", linetype='dashed')+
  ggtitle('Complex Traits')

p1

# Save off part of supplementary table 10
traits_out %>% dplyr::filter(cohort == 'Complex Traits',
                             annot %in% c("is.e2g","e2g_yes_emVar_yes",
                                          "is.CRE", "CRE_yes_emVar_yes")) %>%
  vroom::vroom_write(paste0('tables/s9f.txt.gz'))

# Convert to a ggplot and save
leg <- get_legend(p1)
as_ggplot(leg)

cowplot::save_plot(paste0('figures/s9f.pdf'),
                   p1+theme(legend.position = 'None'),
                   base_height = 6,
                   base_width = 6,
                   device = cairo_pdf
)

cowplot::save_plot(paste0('figures/s9f.legend.pdf'),
                   leg,
                   base_height = 6,
                   base_width = 6,
                   device = cairo_pdf
)
