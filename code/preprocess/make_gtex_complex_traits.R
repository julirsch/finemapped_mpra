unixtools:: set.tempdir("/mnt/sdb/rstudio/sirajl/")

# set working directory
setwd('/mnt/sdb/gwas_eqtl_mpra/')
source("code/github_ready/Preprocess/read_in_functions.R")
source("code/utils.R")

library(vroom)
library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(plyranges)
library(BuenColors)
library(rtracklayer)
library(reshape2)
library(ComplexHeatmap)
library(binom)
library(ggExtra)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE)

# read in all annotation files
# Read in general annotations
# 
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

# Read in data for GTEx
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
                                   TRUE ~ paste0(cohort, ";", tissue, ";", gene, ";", cs_id))) %>%
  dplyr::mutate(position = as.numeric(str_split_fixed(variant, ":", 4)[, 2]))
# annotate with GTEx system info
#GTEx.info <- read_delim("/mnt/sdb/gtex_mpra/data/GTEx_49tissues.txt", delim = "\t", col_names = T)
GTEx.info <- read_delim("/mnt/sdb/gwas_eqtl_mpra/data/gtex/GTEx_49tissues.txt", delim = "\t", col_names = T)
mpra <- mpra %>% 
  left_join(GTEx.info %>%
              dplyr::select(tissue, system),
            by = "tissue")
rm(GTEx.info)

# drop columns to be re-done for GTEx
mpra <- mpra %>% 
  dplyr::select(-c('consequence','consequence2'),
                -c("coding", "UTR5", "UTR3", "conserved", "promoter", "intron", "DHS", "K27ac_1", "K27ac_2", "DHS_Roadmap", "DHS_Meuleman", "DHS_Domcke", "DHS_Corces", "DHS_Buen", "DHS_Calderon", "DHS_ChipAtlas"),
                -c('CRE','CRE_Roadmap','CRE_Meuleman','CRE_Domcke','CRE_Corces','CRE_Buen','CRE_Calderon','CRE_ChipAtlas'))

# Annotate GTEx variants
variant.gr <- mpra %>%
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

mpra <- mpra %>%
  left_join(var.df,
            by = "variant")
rm(variant.gr); rm(var.df); rm(annot.df)
gc()

# add consequence annotations to GTEx data
mpra <- mpra %>%
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

mpra <- def_emVars(mpra)
mpra <- def_intersecting(mpra, "")
variant.gr <- get_variants(mpra)
mpra <- add_footprints(mpra, variant.gr)
mpra <- add_footprints_5(mpra, variant.gr)
mpra <- def_intersecting_2(mpra)
# save off
mpra %>% vroom::vroom_write(file = "./gtex.txt.gz", delim = "\t", col_names = T)



# repeat all of the above for complex traits data
# read in data
cohort = 'traits'
paired = 'paired'
path <- "/mnt/sdb/gwas_eqtl_mpra/data/preprocess/gridsearch/"
#path <- "/mnt/sdb/gtex_mpra/final_mpra_data/processed/"
print(paste0(path,cohort,"_mpra_",paired,"_final20230117.rds"))
traits_mpra_df <- readRDS(paste0(path,cohort,"_mpra_",paired,"_final20230117.rds")) %>% 
  dplyr::filter(best_library == 1) %>%
  dplyr::filter(type != "other_test") %>% 
  dplyr::filter(cell_type != "GM12878") %>%
  dplyr::rename(log2Skew = Log2Skew) %>% 
  dplyr::mutate(cs_uid = case_when(cohort == 'control' ~ NA_character_,
                                   cs_id == -1 ~ NA_character_,
                                   TRUE ~ paste0(cohort, ";", trait, ";", region, ";", cs_id))) %>%
  dplyr::mutate(position = as.numeric(str_split_fixed(variant, ":", 4)[, 2]))

traits_mpra_df <- traits_mpra_df %>% 
  dplyr::select(-c('consequence','consequence2'),
                -c("coding", "UTR5", "UTR3", "conserved", "promoter", "intron", "DHS", "K27ac_1", "K27ac_2", "DHS_Roadmap", "DHS_Meuleman", "DHS_Domcke", "DHS_Corces", "DHS_Buen", "DHS_Calderon", "DHS_ChipAtlas"),
                -c('CRE','CRE_Roadmap','CRE_Meuleman','CRE_Domcke','CRE_Corces','CRE_Buen','CRE_Calderon','CRE_ChipAtlas'))

# Annotate complex trait variants
variant.gr <- traits_mpra_df %>%
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

traits_mpra_df <- traits_mpra_df %>%
  left_join(var.df,
            by = "variant")
rm(variant.gr); rm(var.df); rm(annot.df)
gc()

# add consequence annotations to complex trait variants
traits_mpra_df <- traits_mpra_df %>%
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

traits_mpra_df <- def_emVars(traits_mpra_df)
traits_mpra_df <- def_intersecting(traits_mpra_df, "")
variant.gr <- get_variants(traits_mpra_df)
traits_mpra_df <- add_footprints(traits_mpra_df, variant.gr)
traits_mpra_df <- add_footprints_5(traits_mpra_df, variant.gr)
traits_mpra_df <- def_intersecting_2(traits_mpra_df)

# save off
traits_mpra_df %>% vroom::vroom_write(file = "./complex_traits.txt.gz", delim = "\t", col_names = T)


rm(coding.gr); rm(UTR5.gr); rm(UTR3.gr); rm(conserved.gr); rm(DHS_Roadmap.gr); rm(H3K27ac.gr); rm(promoter.gr); rm(intron.gr); rm(DHS.gr); rm(K27_1.gr); rm(K27_2.gr)
rm(DHS_Meuleman.gr); rm(DHS_Domcke.gr); rm(DHS_Corces.gr); rm(DHS_Buen.gr); rm(DHS_Calderon.gr); rm(CA_DHS.gr); rm(CA_DHS.q500.gr)
rm(annotations.gr.list);
gc()
