# Setup output
outdir <- "figures/fig3/"

# Load in libraries
source("code/utils.R")
library(vroom)
library(tidyverse)
library(BuenColors)
library(ggplot2)
library(dplyr)
library(binom)
library(PNWColors)
library(stringr)
library(rtracklayer)
library(patchwork)
options(stringsAsFactors = FALSE)


# Set up plotting theme
my_theme <-
  BuenColors::pretty_plot(fontsize = 12) +
  BuenColors::L_border() +
  theme(
    plot.background = element_blank(),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"),
    plot.title = element_text(hjust = 4e-3, margin = margin(b = -12)),
    # legend.position = "none",
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.title = element_text(margin = margin(0.0, 0.0, 0.0, 0.0)),
    legend.background = element_blank(),
    legend.key.size = unit(0.2, "cm")
  )

# Read in MPRA data
# UKBB / BBJ
traits_mpra_df <- readRDS("/mnt/sdb/gtex_mpra/final_mpra_data/processed/traits_mpra_paired_final20230117.rds") %>%
  dplyr::filter(cohort %in% c("UKBB", "GTEx", "control")) %>%
  dplyr::filter(type %in% c("CS", "loc_CS", "null_CS", "annot_CS")) %>%
  dplyr::filter(cell_type != "GM12878") %>%
  dplyr::rename(log2Skew = Log2Skew) %>%
  dplyr::mutate(position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  ungroup() %>%
  tidyr::replace_na(list(Skew_logPadj = 0)) %>%
  dplyr::mutate(active = ifelse(logPadj_BF >= -log10(0.01) & abs(log2FC) >= 1, T, F),
                emVar = ifelse(active & Skew_logPadj >= -log10(0.1) & !is.na(Skew_logPadj) & abs(log2Skew) >= 0, T, F)) %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(active_any = ifelse(any(active), T, F),
                emVar_any = ifelse(any(emVar), T, F),
                emVar_all = ifelse(all(emVar),T,F)) %>%
  dplyr::mutate(cs_uid = paste0(cohort, ";", trait, ";", region, ";", cs_id))

# GTEx
eqtl_mpra_df <- readRDS("/mnt/sdb/gtex_mpra/final_mpra_data/processed/gtex_mpra_paired_final20230117.rds") %>%
  dplyr::mutate(cohort = ifelse(is.na(cohort), "control", cohort)) %>%
  dplyr::filter(cohort %in% c("GTEx", "control")) %>%
  dplyr::filter(type %in% c("3tissue_CS", "3tissue_locctrl", "3tissue_annotctrl")) %>%
  dplyr::rename(log2Skew = Log2Skew) %>%
  dplyr::mutate(position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  ungroup() %>%
  tidyr::replace_na(list(Skew_logPadj = 0)) %>%
  dplyr::mutate(active = ifelse(logPadj_BF >= -log10(0.01) & abs(log2FC) >= 1, T, F),
                emVar = ifelse(active & Skew_logPadj >= -log10(0.1) & !is.na(Skew_logPadj) & abs(log2Skew) >= 0, T, F)) %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(active_any = ifelse(any(active), T, F),
                emVar_any = ifelse(any(emVar), T, F),
                emVar_all = ifelse(all(emVar),T,F)) %>%
  dplyr::mutate(cs_uid = paste0(cohort, ";", tissue, ";", gene, ";", cs_id))


# Annotate control variants
# Read in general annotations
coding.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/Coding_UCSC.bed", format="bed")
UTR5.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/UTR_5_UCSC.bed", format="bed")
UTR3.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/UTR_3_UCSC.bed", format="bed")
conserved.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/Conserved_LindbladToh.bed", format="bed")
DHS_Roadmap.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/DHS_peaks_Trynka.bed", format="bed")
H3K27ac.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/H3K27ac_PGC2.bed", format="bed")
promoter.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/Promoter_UCSC.bed", format="bed")
intron.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/Intron_UCSC.bed", format="bed")
DHS.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[1]]
K27_1.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[2]]
K27_2.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[3]]

# Read in specific chromatin accesibility annotations
DHS_Meuleman.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/meuleman/meuleman_en_gr.rds")
DHS_Domcke.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/domcke/domcke_en_gr.rds")
DHS_Corces.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/corces/corces_en_gr.rds")
DHS_Buen.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/buen/buen_en_gr.rds")
DHS_Calderon.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/calderon/calderon_en_gr.rds")
CA_DHS.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/ChIPAtlas.dnase.rds")
CA_DHS.q500.gr <- CA_DHS.gr[CA_DHS.gr$qval > 500] %>%
  GRsize(., shift = 150) %>%
  GenomicRanges::reduce()

# Create annotation GR list
annotations.gr.list <- list(coding.gr, UTR5.gr, UTR3.gr, conserved.gr, promoter.gr, intron.gr, DHS.gr, K27_1.gr, K27_2.gr, DHS_Roadmap.gr, DHS_Meuleman.gr, DHS_Domcke.gr, DHS_Corces.gr, DHS_Buen.gr, DHS_Calderon.gr, CA_DHS.q500.gr)
names(annotations.gr.list) <- c("coding", "UTR5", "UTR3", "conserved", "promoter", "intron", "DHS", "K27ac_1", "K27ac_2", "DHS_Roadmap", "DHS_Meuleman", "DHS_Domcke", "DHS_Corces", "DHS_Buen", "DHS_Calderon", "DHS_ChipAtlas")

# Annotate GTEx variants
variant.gr <- bind_rows(traits_mpra_df,
                        eqtl_mpra_df) %>%
  distinct(variant, chromosome, position) %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "position", 
                           end.field = "position",
                           keep.extra.columns = T)
annot.df <- data.frame(lapply(1:length(annotations.gr.list), function(x) {countOverlaps(variant.gr, annotations.gr.list[[x]])}))
names(annot.df) <- names(annotations.gr.list)
var.df <- bind_cols(variant.gr %>% as_tibble() %>% dplyr::select(variant), annot.df) %>%
  dplyr::mutate(CRE_Roadmap = ifelse(DHS_Roadmap >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_Meuleman = ifelse(DHS_Meuleman >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_Domcke = ifelse(DHS_Domcke >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_Corces = ifelse(DHS_Corces >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_Buen = ifelse(DHS_Buen >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_Calderon = ifelse(DHS_Calderon >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_ChipAtlas = ifelse(DHS_ChipAtlas >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE = ifelse(CRE_Meuleman + CRE_Domcke + CRE_Corces + CRE_Buen + CRE_Calderon + CRE_ChipAtlas >= 1, 1, 0)) %>%
  dplyr::select(variant, CRE2 = CRE)

eqtl_mpra_df <- eqtl_mpra_df %>%
  left_join(var.df,
            by = "variant") %>%
  dplyr::mutate(CRE = CRE2) %>%
  dplyr::select(-CRE2) 
eqtl_mpra_df <- eqtl_mpra_df %>%
  dplyr::mutate(CRE_emVar = ifelse(CRE > 0 & emVar == 1, T, F))
traits_mpra_df <- traits_mpra_df %>%
  left_join(var.df,
            by = "variant") %>%
  dplyr::mutate(CRE = CRE2) %>%
  dplyr::select(-CRE2)
traits_mpra_df <- traits_mpra_df %>%
  dplyr::mutate(CRE_emVar = ifelse(CRE > 0 & emVar == 1, T, F)) 
traits_mpra_df <- traits_mpra_df %>%
  dplyr::mutate(CRE_active = ifelse(CRE > 0 & active == 1, T, F)) 
eqtl_mpra_df <- eqtl_mpra_df %>%
  dplyr::mutate(CRE_active = ifelse(CRE > 0 & active == 1, T, F)) 

# Add annotation of near variants
cre_variants <- bind_rows(traits_mpra_df, eqtl_mpra_df) %>%
  filter(CRE == 1) %>%
  pull(variant) %>%
  unique()

# Step 1: Filter the GRanges object
filtered_gr <- variant.gr[variant.gr$variant %in% cre_variants]

# Step 2: Flank the filtered GRanges object
flanked_gr <- flank(filtered_gr, width = 150, both = TRUE)

# Step 3: Reduce the flanked GRanges object
CRE.gr <- GenomicRanges::reduce(flanked_gr)

max(width(CRE.gr))

idx <- findOverlaps(variant.gr, CRE.gr) 
CRE_df <- cbind(variant.gr[idx@from]$variant, paste0("CRE_", idx@to)) 
colnames(CRE_df) <- c("variant", "CRE_id")
CRE_df <- CRE_df %>%
  as_tibble()


# Compute enrichment of annotation in CS over background
mcv_enrich <- function(df, bset, fset, annot, cssize, r2){
  
  # Get background hits
  background <- df %>%
    ungroup() %>%
    dplyr::filter(type %in% bset) %>% 
    group_by(cs_uid) %>%
    dplyr::filter(any(consequence2 %ni% c("LoF", "missense", "synonymous"))) %>%
    ungroup() %>%
    dplyr::mutate(keep = ifelse(type %in% c("CS", "3tissue_CS"), (pip < 0.01) & (cs_id > 0), T)) %>%
    dplyr::filter(keep == T) %>%
    #dplyr::filter(pip < 0.01, cs_id > 0) %>%
    dplyr::mutate(annotation = ifelse(!! sym(annot) == T, T, F)) %>%
    dplyr::select(variant, library, cell_type, annotation) %>%
    distinct() %>%
    dplyr::count(library, cell_type, annotation) %>%
    group_by(library, cell_type) %>%
    dplyr::summarize(b_denom = sum(n),
                     b_num = sum(n[annotation == T]))
  
  # Get CS hits
  foreground_tmp <- df %>%
    ungroup() %>%
    dplyr::filter(type %in% fset) %>% 
    group_by(cs_uid) %>%
    dplyr::filter(any(consequence2 %ni% c("LoF", "missense", "synonymous"))) %>%
    ungroup() %>%
    #dplyr::filter(type %in% c("3tissue_locctrl")) %>%
    #dplyr::filter(type %in% c("3tissue_CS"), cs_n < 10, cs_min_r2 > 0.8, pip > 0.01, cs_id > 0) %>%
    dplyr::filter(cs_min_r2 > r2, cs_id > 0) %>%
    dplyr::mutate(annotation = ifelse(!! sym(annot) == T, T, F)) %>%
    dplyr::select(cs_uid, variant, pip, library, cell_type, annotation) %>%
    distinct() %>%
    dplyr::group_by(cs_uid, library, cell_type) %>%
    dplyr::mutate(sum_pip = sum(pip),
                  vrank = rank(pip, ties.method = 'first')) %>%
    ungroup() %>% 
    dplyr::filter(sum_pip > 0.9) %>%
    dplyr::filter(pip > 0.05, vrank <= cssize)
  foreground <- foreground_tmp %>%
    dplyr::group_by(cs_uid, library, cell_type) %>%
    dplyr::mutate(cs_n = length(variant),
                  has_annotation = sum(annotation) >= 1) %>%
    ungroup() %>%
    dplyr::filter(cs_n > 1, has_annotation == T) %>%
    #left_join(CRE_df,
    #          by = "variant") %>%
    #dplyr::mutate(CRE_cs_id = paste0(CRE_id, "_", cs_uid)) %>%
    group_by(library, cell_type) %>%
    dplyr::summarize(f_cs = length(unique(cs_uid[has_annotation == 1])),
                     f_denom = sum(has_annotation == 1),
                     f_num = sum(has_annotation == 1 & annotation == 1))#,
  #f_num2 = length(unique(CRE_cs_id[has_annotation == 1 & annotation == 1])),
  #f_num3 = sum(duplicated(CRE_cs_id[has_annotation == 1 & annotation == 1])) + f_cs)
  # Compute RRs
  pre_input <- foreground %>%
    left_join(background,
              by = c("library", "cell_type")) %>%
    group_by(library, cell_type) %>%
    dplyr::mutate(a = f_num - 1 * f_cs, 
                  b = f_denom - f_num,
                  c = b_num,
                  d = b_denom - b_num) 
  
  input_RR <- pre_input %>%
    metafor::escalc(measure = "RR", ai = a, bi = b, ci = c, di = d, data = .)
  input_RD <- pre_input %>%
    metafor::escalc(measure = "RD", ai = a, bi = b, ci = c, di = d, data = .)
  
  # Meta analyze (RR)
  input_RR %>%
    metafor::rma(yi, vi, data = ., digits = 5, control = list(stepadj = 0.5, maxiter = 1000)) %>%
    broom::tidy(conf.int = T,
                conf.level = 0.95,
                exponentiate = T) %>%
    dplyr::mutate(cohort = case_when(!! (fset) %in% c("CS") ~ "GWAS",
                                     !! (fset) %in% c("3tissue_CS") ~ "eQTL"),
                  annotation = annot,
                  control = case_when(!! (bset) %in% c("3tissue_annotctrl", "annot_CS") ~ "annotation",
                                      !! (bset) %in% c("3tissue_locctrl", "loc_CS") ~ "location",
                                      !! (bset) %in% c("3tissue_CS", "CS") ~ "low pip"),
                  cs_size = cssize,
                  r_2 = r2) %>%
    dplyr::select(cohort, annotation, control, cs_size, r_2, estimate, conf.low, conf.high, p.value) %>%
    bind_cols(input_RD %>%
                metafor::rma(yi, vi, data = ., digits = 5, control = list(stepadj = 0.5, maxiter = 1000)) %>%
                broom::tidy(conf.int = T,
                            conf.level = 0.95,
                            exponentiate = F) %>%
                bind_cols(pre_input %>%
                            ungroup() %>%
                            dplyr::summarize(n_vars = sum(b),
                                             n_CSs = sum(f_cs))) %>%
                dplyr::mutate(excess = n_vars * estimate, eprop = excess / n_CSs, eprop.low = n_vars * conf.low / n_CSs, eprop.high = n_vars * conf.high / n_CSs) %>%
                dplyr::select(n_CSs, excess, eprop, eprop.low, eprop.high)) %>%
    bind_cols(df %>%
                dplyr::filter(type %in% c("3tissue_CS", "CS")) %>% 
                group_by(cs_uid) %>%
                dplyr::filter(any(consequence2 %ni% c("LoF", "missense", "synonymous"))) %>%
                ungroup() %>%
                dplyr::filter(cs_min_r2 > r2, cs_id > 0) %>%
                dplyr::select(cs_uid, variant, pip, CRE, emVar_any) %>%
                distinct() %>%
                dplyr::mutate(CRE_emVar = ifelse(CRE + emVar_any == 2, T, F)) %>%
                dplyr::group_by(cs_uid) %>%
                dplyr::mutate(sum_pip = sum(pip),
                              vrank = rank(pip, ties.method = 'first')) %>%
                ungroup() %>% 
                dplyr::filter(sum_pip > 0.9) %>%
                dplyr::filter(pip > 0.05, vrank <= cssize) %>%
                dplyr::group_by(cs_uid) %>%
                dplyr::mutate(cs_n = length(variant),
                              has_annotation = sum(CRE_emVar) > 1) %>%
                ungroup() %>%
                dplyr::filter(cs_n > 1, has_annotation == T) %>%
                distinct(cs_uid, variant, CRE_emVar) %>%
                left_join(CRE_df,
                          by = "variant") %>%
                dplyr::mutate(CRE_cs_id = paste0(CRE_id, "_", cs_uid)) %>%
                dplyr::filter(CRE_emVar == T) %>%
                #group_by(cs_uid) %>%
                #dplyr::summarize(num_CREs = length(unique(CRE_id)), num_CRE_emVars = length(CRE_id)) %>%
                #dplyr::count(num_CRE_emVars - num_CREs > 0)
                distinct(variant, CRE_id) %>%
                dplyr::count(num_CREs = length(unique(CRE_id)), num_CRE_emVars = length(CRE_id)) %>%
                dplyr::mutate(prop_uCRE = num_CREs / num_CRE_emVars) %>%
                dplyr::select(prop_uCRE))
}


# Read in processed MPRA results
mpra_df <- vroom("/mnt/sdb/gwas_eqtl_mpra/data/preprocess/core_mpra.txt.gz")

# read in credible set metadata
files <- list.files("/mnt/sdb/gtex_mpra/data/merged_CSs/", pattern = "*.txt.gz", full.names = T)
csm_full <- vroom::vroom(files, col_names = T)

mpra_df <- mpra_df %>%
  left_join(csm_full %>%
              dplyr::rename("cs_uid" = "id") %>%
              distinct())

# r2 <- 0.8
# blah <- mpra_df %>% 
#   dplyr::filter(cohort == "GTEx", type %in% c("3tissue_CS", "3tissue_PIP10")) %>% 
#   dplyr::group_by(cs_uid) %>% 
#   dplyr::filter(any(consequence2 %ni% c("LoF", "missense", "synonymous"))) %>% 
#   ungroup() %>%
#   dplyr::filter(cs_min_r2 > r2, cs_id > 0, !is.na(cs_uid), !is.na(csm_id)) %>%
#   dplyr::select(tissue, gene, cs_id, cs_uid, csm_id, variant, pip, cell_type, emVar) %>% 
#   distinct() %>%
#   dplyr::group_by(cell_type, cs_uid) %>% 
#   dplyr::mutate(sum_pip = sum(pip),
#                 vrank = rank(pip, ties.method = 'first')) %>%
#   ungroup() %>% 
#   dplyr::filter(sum_pip > 0.9) %>%
#   dplyr::filter(pip > 0.01, vrank <= cssize) %>%
#   dplyr::group_by(cell_type, tissue, gene, cs_id) %>% 
#   dplyr::mutate(n_vars = sum(pip > 0.01)) %>%
#   dplyr::filter(n_vars >= 2) %>%
#   ungroup()
# blah %>%
#   distinct(csm_id)
# blah %>%
#   distinct(variant)
# blah %>%
#   distinct(csm_id, variant) %>%
#   group_by(csm_id) %>%
#   dplyr::summarize(cs_n = length(csm_id)) %>%
#   ungroup() %>%
#   dplyr::summarize(median = median(cs_n))
# blah %>%
#   distinct(csm_id, variant, cell_type, emVar) %>%
#   dplyr::group_by(cell_type, csm_id) %>% 
#   #dplyr::filter(length(csm_id) <= 30) %>%
#   dplyr::summarize(mcvs = sum(emVar)) %>% 
#   ungroup() %>% 
#   #dplyr::filter(mcvs >= 1) %>%
#   group_by(cell_type) %>%
#   dplyr::summarize(mcvs_n = sum(mcvs >= 2),
#                    cs_n = length(unique(csm_id))) %>%
#   ungroup() %>% 
#   dplyr::mutate(prop = mcvs_n / cs_n)

mcvs_df <- bind_rows(mcv_enrich(traits_mpra_df, "CS", "CS", "CRE_emVar", 5, 0.9),
                         mcv_enrich(traits_mpra_df, "CS", "CS", "CRE", 5, 0.9),
                         mcv_enrich(traits_mpra_df, "CS", "CS", "emVar", 5, 0.9),
                         mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE_emVar", 5, 0.9),
                         mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE", 5, 0.9),
                         mcv_enrich(traits_mpra_df, "loc_CS", "CS", "emVar", 5, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE_emVar", 5, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE", 5, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "emVar", 5, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE_emVar", 5, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE", 5, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "emVar", 5, 0.9),
                         mcv_enrich(traits_mpra_df, "CS", "CS", "CRE_emVar", 10, 0.9),
                         mcv_enrich(traits_mpra_df, "CS", "CS", "CRE", 10, 0.9),
                         mcv_enrich(traits_mpra_df, "CS", "CS", "emVar", 10, 0.9),
                         mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE_emVar", 10, 0.9),
                         mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE", 10, 0.9),
                         mcv_enrich(traits_mpra_df, "loc_CS", "CS", "emVar", 10, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE_emVar", 10, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE", 10, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "emVar", 10, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE_emVar", 10, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE", 10, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "emVar", 10, 0.9),
                         mcv_enrich(traits_mpra_df, "CS", "CS", "CRE_emVar", 15, 0.9),
                         mcv_enrich(traits_mpra_df, "CS", "CS", "CRE", 15, 0.9),
                         mcv_enrich(traits_mpra_df, "CS", "CS", "emVar", 15, 0.9),
                         mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE_emVar", 15, 0.9),
                         mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE", 15, 0.9),
                         mcv_enrich(traits_mpra_df, "loc_CS", "CS", "emVar", 15, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE_emVar", 15, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE", 15, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "emVar", 15, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE_emVar", 15, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE", 15, 0.9),
                         mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "emVar", 15, 0.9),
                         mcv_enrich(traits_mpra_df, "CS", "CS", "CRE_emVar", 5, 0.8),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "CRE", 5, 0.8),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "emVar", 5, 0.8),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE_emVar", 5, 0.8),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE", 5, 0.8),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "emVar", 5, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE_emVar", 5, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE", 5, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "emVar", 5, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE_emVar", 5, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE", 5, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "emVar", 5, 0.8),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "CRE_emVar", 10, 0.8),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "CRE", 10, 0.8),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "emVar", 10, 0.8),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE_emVar", 10, 0.8),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE", 10, 0.8),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "emVar", 10, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE_emVar", 10, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE", 10, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "emVar", 10, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE_emVar", 10, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE", 10, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "emVar", 10, 0.8),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "CRE_emVar", 15, 0.8),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "CRE", 15, 0.8),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "emVar", 15, 0.8),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE_emVar", 15, 0.8),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE", 15, 0.8),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "emVar", 15, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE_emVar", 15, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE", 15, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "emVar", 15, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE_emVar", 15, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE", 15, 0.8),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "emVar", 15, 0.8),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "CRE_emVar", 5, 0.99),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "CRE", 5, 0.99),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "emVar", 5, 0.99),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE_emVar", 5, 0.99),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE", 5, 0.99),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "emVar", 5, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE_emVar", 5, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE", 5, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "emVar", 5, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE_emVar", 5, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE", 5, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "emVar", 5, 0.99),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "CRE_emVar", 10, 0.99),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "CRE", 10, 0.99),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "emVar", 10, 0.99),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE_emVar", 10, 0.99),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE", 10, 0.99),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "emVar", 10, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE_emVar", 10, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE", 10, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "emVar", 10, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE_emVar", 10, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE", 10, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "emVar", 10, 0.99),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "CRE_emVar", 15, 0.99),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "CRE", 15, 0.99),
                       mcv_enrich(traits_mpra_df, "CS", "CS", "emVar", 15, 0.99),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE_emVar", 15, 0.99),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "CRE", 15, 0.99),
                       mcv_enrich(traits_mpra_df, "loc_CS", "CS", "emVar", 15, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE_emVar", 15, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "CRE", 15, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_CS", "3tissue_CS", "emVar", 15, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE_emVar", 15, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "CRE", 15, 0.99),
                       mcv_enrich(eqtl_mpra_df, "3tissue_locctrl", "3tissue_CS", "emVar", 15, 0.99),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE_emVar", 5, 0.8),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE", 5, 0.8),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "emVar", 5, 0.8),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE_emVar", 5, 0.8),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE", 5, 0.8),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "emVar", 5, 0.8),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE_emVar", 10, 0.8),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE", 10, 0.8),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "emVar", 10, 0.8),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE_emVar", 10, 0.8),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE", 10, 0.8),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "emVar", 10, 0.8),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE_emVar", 15, 0.8),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE", 15, 0.8),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "emVar", 15, 0.8),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE_emVar", 15, 0.8),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE", 15, 0.8),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "emVar", 15, 0.8),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE_emVar", 5, 0.9),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE", 5, 0.9),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "emVar", 5, 0.9),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE_emVar", 5, 0.9),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE", 5, 0.9),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "emVar", 5, 0.9),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE_emVar", 10, 0.9),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE", 10, 0.9),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "emVar", 10, 0.9),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE_emVar", 10, 0.9),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE", 10, 0.9),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "emVar", 10, 0.9),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE_emVar", 15, 0.9),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE", 15, 0.9),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "emVar", 15, 0.9),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE_emVar", 15, 0.9),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE", 15, 0.9),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "emVar", 15, 0.9),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE_emVar", 5, 0.99),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE", 5, 0.99),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "emVar", 5, 0.99),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE_emVar", 5, 0.99),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE", 5, 0.99),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "emVar", 5, 0.99),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE_emVar", 10, 0.99),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE", 10, 0.99),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "emVar", 10, 0.99),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE_emVar", 10, 0.99),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE", 10, 0.99),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "emVar", 10, 0.99),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE_emVar", 15, 0.99),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "CRE", 15, 0.99),
                     mcv_enrich(traits_mpra_df, "annot_CS", "CS", "emVar", 15, 0.99),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE_emVar", 15, 0.99),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "CRE", 15, 0.99),
                     mcv_enrich(eqtl_mpra_df, "3tissue_annotctrl", "3tissue_CS", "emVar", 15, 0.99))


#test <- vroom::vroom('/mnt/sdb/gwas_eqtl_mpra/tables/stable18_mcvs_v3.txt.gz')

mcvs_df %>%
  group_by(cohort, control, cs_size, r_2) %>%
  #dplyr::filter(annotation != "emVar") %>%
  dplyr::mutate(diff = log(estimate[annotation == "CRE_emVar"]) - log(estimate[annotation == "CRE"]),
                SE1 = (log(conf.high[annotation == "CRE_emVar"]) - log(conf.low[annotation == "CRE_emVar"])) / 2 / 1.96,
                SE2 = (log(conf.high[annotation == "CRE"]) - log(conf.low[annotation == "CRE"])) / 2 / 1.96,
                SE = sqrt(SE1^2 + SE2^2),
                diff_pval = pnorm(abs(diff) / SE, lower.tail = F) * 2,
                lower = exp(diff  - 1.96 * SE) - 1, 
                upper = exp(diff  + 1.96 * SE) - 1,
                diff = exp(diff) - 1) %>%
  dplyr::select(-excess, -prop_uCRE, -SE1, -SE2, -SE) %>%
  dplyr::mutate(diff = ifelse(annotation == "CRE_emVar", diff, NA),
                diff_pval = ifelse(annotation == "CRE_emVar", diff_pval, NA),
                lower = ifelse(annotation == "CRE_emVar", lower, NA),
                upper = ifelse(annotation == "CRE_emVar", upper, NA)) %>%
  vroom::vroom_write("/mnt/sdb/gwas_eqtl_mpra/tables/stable21_mcvs.txt.gz")

#mcvs_df <- vroom::vroom('/mnt/sdb/gwas_eqtl_mpra/tables/stable18_mcvs_v3.txt.gz')
  
mcvs_df %>% 
  dplyr::filter(annotation == "CRE_emVar", control == "low pip") %>%
  #dplyr::filter(p.value < 0.05 & estimate > 1) %>%
  arrange(eprop) %>% 
  print(n = 100)

mcvs_df %>% 
  dplyr::filter(annotation == "CRE_emVar", control == "low pip") %>%
  #dplyr::filter(p.value < 0.05 & estimate > 1) %>%
  arrange(estimate) %>%
  print(n = 100)

mcvs_df %>% 
  dplyr::filter(annotation == "CRE", control == "low pip") %>%
  #dplyr::filter(p.value < 0.05 & estimate > 1) %>%
  arrange(estimate) %>% 
  print(n = 100)

mcvs_df %>% 
  dplyr::filter(annotation == "emVar", control == "low pip") %>%
  #dplyr::filter(estimate > 1, p.value < 0.05) %>%
  arrange(estimate) %>% 
  print(n = 100)

mpra_df %>%
  dplyr::filter(cohort == "UKBB", trait == "ALP", cs_id == 1, region == "chr3:121565778-124565778") %>%
  distinct(variant, pip, emVar_any, CRE)

traits_mpra_df %>%
  dplyr::distinct(variant, trait, region, cs_id, CRE, emVar_any) %>%
  dplyr::group_by(trait, region, cs_id) %>%
  dplyr::filter(length(cs_id) < 20) %>%
  dplyr::filter(CRE == 1 & emVar_any == 1) %>%
  dplyr::filter(length(CRE) > 4)


p1 <- mcvs_df %>%
  dplyr::mutate(annotation = factor(annotation, c("CRE_emVar", "CRE", "emVar")),
                control = factor(control, c("low pip", "annotation", "location"))) %>%
  dplyr::filter(cs_size == 5, r_2 == 0.9) %>%
  #dplyr::filter(cs_size == 5, r_2 == 0.9, control) %>%
  ggplot(aes(x = annotation, y = estimate , fill = control)) +
  geom_hline(aes(yintercept = 1)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                size = 0.5,
                width = 0,
                col = "black",
                position = position_dodge(0.9, preserve = "single")) +
  scale_color_manual(values = PNWColors::pnw_palette('Shuksan')[c(2,3,4)],
                     aesthetics = c("color","fill")) +
  pretty_plot() +
  facet_wrap(~cohort, nrow = 2) +
  my_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle=90, hjust=1))
p1

plt_combined <- p1 + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot(
  paste0(outdir, "fig3f_mcvs.pdf"),
  plt_combined,
  base_height = 3,
  base_width = 2,
  device = cairo_pdf
)

# Get background hits
background <- eqtl_mpra_df %>%
  ungroup() %>%
  dplyr::filter(type %in% c("3tissue_CS")) %>%
  dplyr::filter(cell_type != "GM12878") %>%
  #dplyr::filter(type %in% c("3tissue_locctrl")) %>%
  group_by(tissue, gene, cs_id) %>%
  dplyr::filter(any(consequence2 %ni% c("LoF", "missense", "synonymous"))) %>%
  ungroup() %>%
  dplyr::filter(pip < 0.01, cs_id > 0) %>%
  dplyr::mutate(CRE_emVar = ifelse(emVar == T, T ,F)) %>%
  dplyr::select(variant, library, cell_type, CRE_emVar) %>%
  distinct() %>%
  dplyr::count(library, cell_type, CRE_emVar) %>%
  group_by(library, cell_type) %>%
  dplyr::summarize(b_denom = sum(n),
                   b_num = sum(n[CRE_emVar == T])) %>%
  dplyr::mutate(prop = b_num / b_denom)

# Get CS hits
foreground <- eqtl_mpra_df %>%
  dplyr::filter(cell_type != "GM12878") %>%
  dplyr::filter(type %in% c("3tissue_CS")) %>%
  ungroup() %>%
  group_by(tissue, gene, cs_id) %>%
  dplyr::filter(any(consequence2 %ni% c("LoF", "missense", "synonymous"))) %>%
  ungroup() %>%
  #dplyr::filter(type %in% c("3tissue_locctrl")) %>%
  #dplyr::filter(type %in% c("3tissue_CS"), cs_n < 10, cs_min_r2 > 0.8, pip > 0.01, cs_id > 0) %>%
  dplyr::filter(cs_n <= 10, cs_n != 0, cs_min_r2 > 0.6, pip > 0.01, cs_id > 0) %>%
  dplyr::mutate(CRE_emVar = ifelse(emVar == T, T, F)) %>%
  dplyr::select(cs_uid, variant, pip, library, cell_type, CRE_emVar) %>%
  distinct() %>%
  dplyr::group_by(cs_uid, library, cell_type) %>%
  dplyr::filter(sum(pip) > 0.8) %>%
  dplyr::mutate(has_CRE_emVar = sum(CRE_emVar) >= 1) %>%
  ungroup() %>%
  group_by(library, cell_type) %>%
  dplyr::summarize(f_cs = length(unique(cs_uid[has_CRE_emVar == 1])),
                   f_denom = sum(has_CRE_emVar == 1),
                   f_num = sum(has_CRE_emVar == 1 & CRE_emVar == 1))

# Compute ORs
input <- foreground %>%
  left_join(background,
            by = c("library", "cell_type")) %>%
  group_by(library, cell_type) %>%
  dplyr::mutate(a = f_num - 1 * f_cs, 
                b = f_denom - 1 * f_cs,
                c = b_num,
                d = b_denom) %>%
  metafor::escalc(measure = "OR", ai = a, bi = b, ci = c, di = d, data = .)

# Meta analyze
input %>%
  metafor::rma(yi, vi, data = .) %>%
  broom::tidy(conf.int = T,
              conf.level = 0.95,
              exponentiate = T)


