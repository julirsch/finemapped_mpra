# set temp storage
unixtools:: set.tempdir("/mnt/sdb/rstudio/sirajl/")

# set working directory
setwd('/mnt/sdb/gwas_eqtl_mpra')

# source files
source("code/utils.R")
source("code/github_ready/Preprocess/read_in_functions.R")

# load libraries
library(DESeq2)
library(plyranges)
library(tidyverse)
library(GenomicRanges)
library(BuenColors)
library(rtracklayer)
library(reshape2)
library(ComplexHeatmap)
library(binom)
options(stringsAsFactors = FALSE)
library("magrittr")
library(ggExtra)
library(dplyr)
library(ggplot2)
library(TFBSTools)
library(motifmatchr)
library(universalmotif)
library(motifbreakR)
data('hocomoco')

# Read in data
cohort = 'gtex'
paired = 'paired'

#path <- "/mnt/sdb/gtex_mpra/final_mpra_data/processed/"
path <- "/mnt/sdb/gwas_eqtl_mpra/data/preprocess/gridsearch/"
print(paste0(path,cohort,"_mpra_",paired,"_final20230117.rds"))
mpra <- readRDS(paste0(path,'gtex',"_mpra_",'paired',"_final20230117.rds")) %>% 
  dplyr::filter(best_library == 1) %>%
  dplyr::filter(type != "other_test") %>%
  dplyr::rename(log2Skew = Log2Skew)  %>% 
  dplyr::mutate(cs_uid = case_when(cohort == 'control' ~ NA_character_,
                                   cs_id == -1 ~ NA_character_,
                                   TRUE ~ paste0(cohort, ";", tissue, ";", gene, ";", cs_id))) 

cohort = 'traits'
paired = 'paired'

#path <- "/mnt/sdb/gtex_mpra/final_mpra_data/processed/"
path <- "/mnt/sdb/gwas_eqtl_mpra/data/preprocess/gridsearch/"
print(paste0(path,cohort,"_mpra_",paired,"_final20230117.rds"))
traits_mpra_df <- readRDS(paste0(path,cohort,"_mpra_",paired,"_final20230117.rds")) %>% 
  dplyr::filter(best_library == 1) %>%
  dplyr::filter(type != "other_test") %>% 
  dplyr::filter(cell_type != "GM12878") %>%
  dplyr::rename(log2Skew = Log2Skew) %>% 
  dplyr::mutate(cs_uid = case_when(cohort == 'control' ~ NA_character_,
                                   cs_id == -1 ~ NA_character_,
                                   TRUE ~ paste0(cohort, ";", trait, ";", region, ";", cs_id)))

big_mpra <- bind_rows(traits_mpra_df, mpra) %>%
  ungroup()
rm(traits_mpra_df); rm(mpra)

# get positions
big_mpra <- big_mpra %>%
  dplyr::mutate(position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  ungroup()

# annotate with GTEx system info
#GTEx.info <- read_delim("/mnt/sdb/gtex_mpra/data/GTEx_49tissues.txt", delim = "\t", col_names = T)
GTEx.info <- read_delim("/mnt/sdb/gwas_eqtl_mpra/data/gtex/GTEx_49tissues.txt", delim = "\t", col_names = T)
big_mpra <- big_mpra %>% 
  left_join(GTEx.info %>%
              dplyr::select(tissue, system),
            by = "tissue")
rm(GTEx.info)


# add in CRE and consequence annotations for all variants
# drop columns
big_mpra <- big_mpra %>% 
  dplyr::select(-c('consequence','consequence2'),
                -c("coding", "UTR5", "UTR3", "conserved", "promoter", "intron", "DHS", "K27ac_1", "K27ac_2", "DHS_Roadmap", "DHS_Meuleman", "DHS_Domcke", "DHS_Corces", "DHS_Buen", "DHS_Calderon", "DHS_ChipAtlas"),
                -c('CRE','CRE_Roadmap','CRE_Meuleman','CRE_Domcke','CRE_Corces','CRE_Buen','CRE_Calderon','CRE_ChipAtlas'))

# Annotate control variants
# Read in general annotations
#coding.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/Coding_UCSC.bed", format="bed")
coding.gr <- import("/mnt/sdb/gwas_eqtl_mpra/data/annotations/ldsc/Coding_UCSC.bed", format="bed")
# UTR5.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/UTR_5_UCSC.bed", format="bed")
UTR5.gr <- import("/mnt/sdb/gwas_eqtl_mpra/data/annotations/ldsc/UTR_5_UCSC.bed", format="bed")
# UTR3.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/UTR_3_UCSC.bed", format="bed")
UTR3.gr <- import("/mnt/sdb/gwas_eqtl_mpra/data/annotations/ldsc/UTR_3_UCSC.bed", format="bed")
# conserved.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/Conserved_LindbladToh.bed", format="bed")
conserved.gr <- import("/mnt/sdb/gwas_eqtl_mpra/data/annotations/ldsc/Conserved_LindbladToh.bed", format="bed")
# DHS_Roadmap.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/DHS_peaks_Trynka.bed", format="bed")
DHS_Roadmap.gr <- import("/mnt/sdb/gwas_eqtl_mpra/data/annotations/ldsc/DHS_peaks_Trynka.bed", format="bed")
# H3K27ac.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/H3K27ac_PGC2.bed", format="bed")
H3K27ac.gr <- import("/mnt/sdb/gwas_eqtl_mpra/data/annotations/ldsc/H3K27ac_PGC2.bed", format="bed")
# promoter.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/Promoter_UCSC.bed", format="bed")
promoter.gr <- import("/mnt/sdb/gwas_eqtl_mpra/data/annotations/ldsc/Promoter_UCSC.bed", format="bed")
# intron.gr <- import("/mnt/sdb/gtex_mpra/data/annotations/ldsc/Intron_UCSC.bed", format="bed")
intron.gr <- import("/mnt/sdb/gwas_eqtl_mpra/data/annotations/ldsc/Intron_UCSC.bed", format="bed")
# DHS.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[1]]
# K27_1.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[2]]
# K27_2.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[3]]
DHS.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[1]]
K27_1.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[2]]
K27_2.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[3]]

# Read in specific chromatin accesibility annotations
# DHS_Meuleman.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/meuleman/meuleman_en_gr.rds")
# DHS_Domcke.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/domcke/domcke_en_gr.rds")
# DHS_Corces.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/corces/corces_en_gr.rds")
# DHS_Buen.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/buen/buen_en_gr.rds")
# DHS_Calderon.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/calderon/calderon_en_gr.rds")
# CA_DHS.gr <- readRDS("/mnt/sdb/gtex_mpra/data/annotations/ChIPAtlas.dnase.rds")

DHS_Meuleman.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/meuleman/meuleman_en_gr.rds")
DHS_Domcke.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/domcke/domcke_en_gr.rds")
DHS_Corces.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/corces/corces_en_gr.rds")
DHS_Buen.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/buen/buen_en_gr.rds")
DHS_Calderon.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/calderon/calderon_en_gr.rds")
CA_DHS.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/ChIPAtlas.dnase.rds")
CA_DHS.q500.gr <- CA_DHS.gr[CA_DHS.gr$qval > 500] %>%
  GRsize(., shift = 150) %>%
  GenomicRanges::reduce()

# Create annotation GR list
annotations.gr.list <- list(coding.gr, UTR5.gr, UTR3.gr, conserved.gr, promoter.gr, intron.gr, DHS.gr, K27_1.gr, K27_2.gr, DHS_Roadmap.gr, DHS_Meuleman.gr, DHS_Domcke.gr, DHS_Corces.gr, DHS_Buen.gr, DHS_Calderon.gr, CA_DHS.q500.gr)
names(annotations.gr.list) <- c("coding", "UTR5", "UTR3", "conserved", "promoter", "intron", "DHS", "K27ac_1", "K27ac_2", "DHS_Roadmap", "DHS_Meuleman", "DHS_Domcke", "DHS_Corces", "DHS_Buen", "DHS_Calderon", "DHS_ChipAtlas")



# Annotate GTEx variants
variant.gr <- big_mpra %>%
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
                CRE = ifelse(CRE_Meuleman + CRE_Domcke + CRE_Corces + CRE_Buen + CRE_Calderon + CRE_ChipAtlas >= 1, 1, 0))

big_mpra <- big_mpra %>%
  left_join(var.df,
            by = "variant")


# add consequence annotations
big_mpra <- big_mpra %>%
  dplyr::mutate(consequence = case_when(lof %in% c("HC", "OS") ~ "LoF",
                                        lof == "LC" ~ "missense",
                                        vep.most_severe == "missense_variant" ~ "missense",
                                        vep.most_severe == "synonymous_variant" ~ "synonymous",
                                        vep.most_severe == "3_prime_UTR_variant" ~ "UTR3",
                                        vep.most_severe == "5_prime_UTR_variant" ~ "UTR5",
                                        TRUE ~ "non-genic")) %>%
  dplyr::mutate(consequence2 = case_when(consequence == "non-genic" & promoter >= 1 ~ "promoter",
                                         consequence == "non-genic" & CRE >= 1 ~ "CRE",
                                         T ~ consequence))


rm(coding.gr); rm(UTR5.gr); rm(UTR3.gr); rm(conserved.gr); rm(DHS_Roadmap.gr); rm(H3K27ac.gr); rm(promoter.gr); rm(intron.gr); rm(DHS.gr); rm(K27_1.gr); rm(K27_2.gr)
rm(DHS_Meuleman.gr); rm(DHS_Domcke.gr); rm(DHS_Corces.gr); rm(DHS_Buen.gr); rm(DHS_Calderon.gr); rm(CA_DHS.gr); rm(CA_DHS.q500.gr)
rm(annotations.gr.list); rm(variant.gr); rm(var.df); rm(annot.df)
gc()

# Define emVars and emVar intersecting columns
big_mpra <- def_emVars(big_mpra)
big_mpra <- def_intersecting(big_mpra, "")

# add in footprints
variant.gr <- get_variants(big_mpra)

big_mpra <- add_footprints(big_mpra, variant.gr)
big_mpra <- add_footprints_5(big_mpra, variant.gr)

big_mpra <- def_intersecting_2(big_mpra)

big_mpra <- add_chip(big_mpra, variant.gr)
big_mpra <- add_sei(big_mpra, variant.gr)

big_mpra <- big_mpra %>% ungroup() %>%  dplyr::mutate(chip_sums = rowSums(dplyr::select(., starts_with("TF")))) %>% group_by(variant)

# get gr format of variants
variant.gr <- big_mpra %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "position", 
                           end.field = "position",
                           keep.extra.columns = T)

# add in enformer annotations
# Read in processed ENFORMER scores
#/mnt/sdb/gtex_mpra/data/annotations/enformer/enformer_chip_z_ss.rds
enf_chip <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/enformer/enformer_chip_z_ss.rds") %>%
  dplyr::rename("enf_chip_ss" = enf_ss,
                "enf_chip_ss_bestannot" = enf_ss_bestannot) %>%
  dplyr::filter(variant %in% big_mpra$variant)
#/mnt/sdb/gtex_mpra/data/annotations/enformer/enformer_dhs_z_ss.rds
enf_dhs <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/enformer/enformer_dhs_z_ss.rds") %>%
  dplyr::filter(variant %in% big_mpra$variant)
#/mnt/sdb/gtex_mpra/data/annotations/enformer/enformer_chip_z_new_ss.rds
enf_chip_new <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/enformer/enformer_chip_z_new_ss.rds") %>%
  dplyr::rename("enf_chip_gtex_ss" = enf_ss,
                "enf_chip_gtex_ss_bestannot" = enf_ss_bestannot) %>%
  dplyr::rename("variant_hg38" = variant) %>%
  dplyr::filter(variant_hg38 %in% big_mpra$variant_hg38)
#/mnt/sdb/gtex_mpra/data/annotations/enformer/enformer_dhs_z_new_ss.rds
enf_dhs_new <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/enformer/enformer_dhs_z_new_ss.rds") %>%
  dplyr::rename("enf_dhs_gtex_ss" = enf_dhs_ss,
                "enf_dhs_gtex_ss_bestannot" = enf_dhs_ss_bestannot)%>%
  dplyr::rename("variant_hg38" = variant) %>%
  dplyr::filter(variant_hg38 %in% big_mpra$variant_hg38)

# Merge enformer scores in
big_mpra <- big_mpra %>%
  left_join(bind_cols(enf_chip, enf_dhs %>% dplyr::select(-variant)),
            by = "variant") %>%
  left_join(bind_cols(enf_chip_new, enf_dhs_new %>% dplyr::select(-variant_hg38)),
            by = "variant_hg38")
# scores, 4 score columns - batch 1 for each score and batch 2 
# thresholding - get the quantile in control variants - what's the value that 90% (whatever thresh) of control variants are below


# read in motifs
# all motif matches and motifbreaking
#motifbreakr.df <- readRDS('/mnt/sdb/gtex_mpra/code/motifs/motifbreakr.filtered.ls.rds')
motifbreakr.df <- readRDS('/mnt/sdb/gwas_eqtl_mpra/annotations/motif/motifbreakr.filtered.ls.rds')
# occupied motif matches and motifbreaking
#ALL.MB_C.df.match <- readRDS( "/mnt/sdb/gtex_mpra/code/motifs/motifbreakr.occupied.filtered.ls.rds")
ALL.MB_C.df.match <- readRDS( "/mnt/sdb/gwas_eqtl_mpra/annotations/motif/motifbreakr.occupied.filtered.ls.rds")
ALL.MB_C.df.match <- ALL.MB_C.df.match %>% dplyr::mutate(motif_start_10 = motif_start - 10, motif_end_10 = motif_end + 10)

# label occupied motifbreaking variants
motifbreakr.occupied.breaks.gr <- ALL.MB_C.df.match %>% dplyr::filter(effect != 'neut') %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  dplyr::mutate(motif_absolute_start = position + motif_start, motif_absolute_end = position + motif_end) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "motif_absolute_start", 
                           end.field = "motif_absolute_end",
                           keep.extra.columns = T)
motifs <-data.frame(countOverlaps(variant.gr, motifbreakr.occupied.breaks.gr))
names(motifs) <- 'motifbreakr_occupied_overlap'
big_mpra <- bind_cols(big_mpra, motifs)

# label motifbreaking variants
motifbreakr.breaks.gr <- motifbreakr.df %>% dplyr::filter(effect != 'neut') %>% dplyr::select(-c(seqnames, start, end, width, strand, name, name.1)) %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  dplyr::mutate(motif_absolute_start = position + motif_start, motif_absolute_end = position + motif_end) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "motif_absolute_start", 
                           end.field = "motif_absolute_end",
                           keep.extra.columns = T)
motifs <-data.frame(countOverlaps(variant.gr, motifbreakr.breaks.gr))
names(motifs) <- 'motifbreakr_overlap'
big_mpra <- bind_cols(big_mpra, motifs)


# label occupied motifs (with and without a 10bp buffer on either side)
motifbreakr.occupied.gr <- ALL.MB_C.df.match %>% 
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  dplyr::mutate(motif_absolute_start = position + motif_start, motif_absolute_end = position + motif_end) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "motif_absolute_start", 
                           end.field = "motif_absolute_end",
                           keep.extra.columns = T)
motifbreakr.occupied.gr.10 <- ALL.MB_C.df.match %>% 
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  dplyr::mutate(motif_absolute_start = position + motif_start_10, motif_absolute_end = position + motif_end_10) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "motif_absolute_start", 
                           end.field = "motif_absolute_end",
                           keep.extra.columns = T)
names(motifbreakr.occupied.gr) <-motifbreakr.occupied.gr$group_name_fixed
names(motifbreakr.occupied.gr.10) <-motifbreakr.occupied.gr.10$group_name_fixed
m10 <- data.frame(countOverlaps(variant.gr,motifbreakr.occupied.gr.10))
names(m10) <- c("motif_match_10")
big_mpra <- bind_cols(big_mpra,m10)
m <- data.frame(countOverlaps(variant.gr, motifbreakr.occupied.gr))
names(m) <- c("motif_match")
big_mpra <- bind_cols(big_mpra,m)

# Label motifs (with and without a 10bp buffer on either side )
motifbreakr.gr <- motifbreakr.df %>% dplyr::select(-c(seqnames, start, end, width, strand, name, name.1)) %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  dplyr::mutate(motif_absolute_start = position + motif_start, motif_absolute_end = position + motif_end) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "motif_absolute_start", 
                           end.field = "motif_absolute_end",
                           keep.extra.columns = T)
motifbreakr.gr.10 <- motifbreakr.df %>% dplyr::select(-c(seqnames, start, end, width, strand, name, name.1)) %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  dplyr::mutate(motif_absolute_start = position + motif_start-10, motif_absolute_end = position + motif_end+10) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "motif_absolute_start", 
                           end.field = "motif_absolute_end",
                           keep.extra.columns = T)
names(motifbreakr.gr) <- motifbreakr.gr$group_name_fixed
names(motifbreakr.gr.10) <- motifbreakr.gr.10$group_name_fixed
m10 <- data.frame(countOverlaps(variant.gr, motifbreakr.gr.10))
names(m10) <- c("motif_10")
big_mpra <- bind_cols(big_mpra,m10)
m <- data.frame(countOverlaps(variant.gr, motifbreakr.gr))
names(m) <- c("motif")
big_mpra <- bind_cols(big_mpra,m)

# remove other files
rm(ALL.MB_C.df.match)
rm(m); rm(m10); rm(motifs)
rm(motifbreakr.breaks.gr); rm(motifbreakr.df); rm(motifbreakr.gr); rm(motifbreakr.gr.10)
rm(motifbreakr.occupied.breaks.gr); rm(motifbreakr.occupied.gr); rm(motifbreakr.occupied.gr.10)
rm(enf_chip); rm(enf_chip_new); rm(enf_dhs); rm(enf_dhs_new)
rm(variant.gr); rm(ch)

gc()

big_mpra <- big_mpra %>% ungroup() %>% distinct()

big_mpra %>% distinct(cohort)

# subsets to save off
enf_subset <- big_mpra %>% 
  ungroup() %>%
  dplyr::select(variant,starts_with('enf_')) %>% 
  distinct()
chip_subset <- big_mpra %>% 
  ungroup() %>% 
  dplyr::select(variant,starts_with('TF'),'chip_sums') %>% 
  distinct()
sei_subset <- big_mpra %>% 
  ungroup() %>% 
  dplyr::select(variant,starts_with('SEI')) %>% 
  distinct()
motif_subset <- big_mpra %>% 
  ungroup() %>% 
  dplyr::select('variant','motifbreakr_overlap','motifbreakr_occupied_overlap','motif_match_10','motif_match','motif_10','motif') %>% 
  distinct()

core_mpra <- big_mpra %>%
  ungroup() %>% 
  dplyr::select(-starts_with('enf_'), -starts_with('TF'),-starts_with('SEI'),
                -c('motifbreakr_overlap','motifbreakr_occupied_overlap','motif_match_10','motif_match','motif_10','motif','chip_sums')) %>%
  distinct()


# save off subsets
enf_subset %>% vroom::vroom_write(file = "./enf_mpra.txt.gz", delim = "\t", col_names = T)
chip_subset %>% vroom::vroom_write(file = "./chip_mpra.txt.gz", delim = "\t", col_names = T)
sei_subset %>% vroom::vroom_write(file = "./sei_mpra.txt.gz", delim = "\t", col_names = T)
motif_subset %>% vroom::vroom_write(file = "./motif_mpra.txt.gz", delim = "\t", col_names = T)
core_mpra %>% vroom::vroom_write(file = "./core_mpra.txt.gz", delim = "\t", col_names = T)
big_mpra %>% vroom::vroom_write(file = "./big_mpra.txt.gz", delim = "\t", col_names = T)

