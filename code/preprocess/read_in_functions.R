### makes read in functions for the MPRA data
### adapted from /mnt/sdb/gtex_mpra/code/make_dfs/read_in_functions.R
### Dec 18, 2023

library(rtracklayer)
# Manage temporary files
unixtools::set.tempdir("/mnt/sdb/rstudio/sirajl/")
# Load libraries
setwd('/mnt/sdb/gwas_eqtl_mpra/')
##################################
### liftover implement
ch <- import.chain("/mnt/sdb/gtex_mpra/data/liftover/hg38ToHg19.over.chain")

####### make MPRA - including defining emVars
def_emVars <- function(df) {
  df <- df %>%
    # if there isn't a skew_logPadj value, replace it with zero
    tidyr::replace_na(list(Skew_logPadj = 0)) %>%
    # active elements (lenient) are those where the logPadj(BF) is < -log10(0.05) and the effect size magnitude is at least 0.5
    # within that, emVars (lenients) are variants within active elements whose 
    dplyr::mutate(active = ifelse(logPadj_BF >= -log10(0.01) & abs(log2FC) >= 1, T, F),
                  emVar = ifelse(active & Skew_logPadj >= -log10(0.1) & !is.na(Skew_logPadj) & abs(log2Skew) >= 0, T, F)) %>%
    dplyr::group_by(variant) %>%
    dplyr::mutate(active_any = ifelse(any(active), T, F),
                  emVar_any = ifelse(any(emVar), T, F),
                  emVar_all = ifelse(all(emVar),T,F))
  return(df)
}


## define columns we are interested in
#
def_intersecting <- function(df, kind){
  df <- df %>% ungroup() %>%
    dplyr::mutate(CRE_and_mpra_any = case_when(CRE == 1 & emVar_any == T ~T,
                                               CRE == 0 & emVar_any == F ~F,
                                               CRE == 0 & emVar_any == T ~F,
                                               CRE == 1 & emVar_any == F ~F,
                                               T ~ NA ),
                  CRE_and_active_any = case_when(CRE == 1 & active_any == T ~T,
                                                 CRE == 0 & active_any == F ~F,
                                                 CRE == 0 & active_any == T ~F,
                                                 CRE == 1 & active_any == F ~F,
                                                 T ~ NA ),
                  # footprints_mpra_any = case_when(footprints == 1 & emVar_any == T ~T,
                  #                                 footprints == 0 & emVar_any == F ~F,
                  #                                 footprints == 0 & emVar_any == T ~F,
                  #                                 footprints == 1 & emVar_any == F ~F,
                  #                                 T ~ NA ),
                  CRE_and_mpra_all = case_when(CRE == 1 & emVar_all == T ~T,
                                               CRE == 0 & emVar_all == F ~F,
                                               CRE == 0 & emVar_all == T ~F,
                                               CRE == 1 & emVar_all == F ~F,
                                               T ~ NA ),
                  # footprints_mpra_all = case_when(footprints == 1 & emVar_all == T ~T,
                  #                                 footprints == 0 & emVar_all == F ~F,
                  #                                 footprints == 0 & emVar_all == T ~F,
                  #                                 footprints == 1 & emVar_all == F ~F,
                  #                                 T ~ NA ),
                  non_promoter = ifelse((consequence == "non-genic" & promoter == 0), T, F))
  if(kind=='GTEx'){
    df <- df %>% dplyr::mutate(cell_type_CRE_MPRA = case_when(SKNSH_DHS == 1 & (emVar == T & cell_type== 'SKNSH') ~T,
                                                              SKNSH_DHS == 0 & (emVar == T & cell_type== 'SKNSH') ~F,
                                                              SKNSH_DHS == 1 & (emVar == F & cell_type== 'SKNSH') ~F,
                                                              SKNSH_DHS == 0 & (emVar == F & cell_type== 'SKNSH') ~F,
                                                              K562_DHS == 1 & (emVar == T & cell_type== 'K562') ~T,
                                                              K562_DHS == 0 & (emVar == T & cell_type== 'K562') ~F,
                                                              K562_DHS == 1 & (emVar == F & cell_type== 'K562') ~F,
                                                              K562_DHS == 0 & (emVar == F & cell_type== 'K562') ~F,
                                                              HCT116_DHS == 1 & (emVar == T & cell_type== 'HCT116') ~T,
                                                              HCT116_DHS == 0 & (emVar == T & cell_type== 'HCT116') ~F,
                                                              HCT116_DHS == 1 & (emVar == F & cell_type== 'HCT116') ~F,
                                                              HCT116_DHS == 0 & (emVar == F & cell_type== 'HCT116') ~F,
                                                              HEPG2_DHS == 1 & (emVar == T & cell_type== 'HepG2') ~T,
                                                              HEPG2_DHS == 0 & (emVar == T & cell_type== 'HepG2') ~F,
                                                              HEPG2_DHS == 1 & (emVar == F & cell_type== 'HepG2') ~F,
                                                              HEPG2_DHS == 0 & (emVar == F & cell_type== 'HepG2') ~F,
                                                              T ~ NA),
                               
                               cell_type_CRE = case_when(SKNSH_DHS == 1 | HEPG2_DHS == 1 | K562_DHS == 1 | HCT116_DHS ==1 ~ T,
                                                         (SKNSH_DHS == 0 & HEPG2_DHS == 0) & (K562_DHS == 0 & HCT116_DHS == 0) ~F,
                                                         T ~ NA ))}
  df$gsa <- abs(df$log2Skew) * df$Skew_logPadj
  df <- df %>% group_by(variant)
  return(df)
}


#### here - add in footprint and ChIP-seq annotations
###########################################
### make variant names granges object
get_variants <- function(df){
  variant.gr <- df %>%
    dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                  position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
    makeGRangesFromDataFrame(., 
                             seqnames.field = "chromosome", 
                             start.field =  "position", 
                             end.field = "position",
                             keep.extra.columns = T)
  return(variant.gr)
}

# Add footprint annotations
add_footprints <- function(df, variant.gr){
  footprints.gr <- read_delim("/mnt/sdb/ukbb-finemapping/data/annotations/meuleman/consensus_footprints_and_collapsed_motifs_hg38.bed.gz", delim = "\t",
                              col_names = c("seqname", "start", "end", "identifier", "mean_signal", "num_samples", "num_fps", "width", "summit_pos", "core_start", "core_end", "motif_clusters")) %>%
    dplyr::select(seqname, start, end, identifier, mean_signal, motif_clusters) %>% 
    makeGRangesFromDataFrame(., 
                             seqnames.field = "seqname", 
                             start.field =  "start", 
                             end.field = "end", 
                             keep.extra.columns = T) 
  ### lift over to hg19
  footprints.hg19.gr <- liftOver(footprints.gr, ch)
  rm(footprints.gr)
  footprints.hg19.gr <- unlist(footprints.hg19.gr)
  footprints.hg19.gr <- footprints.hg19.gr[!is.na(footprints.hg19.gr$identifier)]
  footprints.hg19.gr <- footprints.hg19.gr[!(duplicated(footprints.hg19.gr$identifier) | duplicated(footprints.hg19.gr$identifier, fromLast = TRUE)),]
  names(footprints.hg19.gr) <- footprints.hg19.gr$identifier
  fp <- data.frame(countOverlaps(variant.gr, footprints.hg19.gr))
  names(fp) <- c("footprints")
  df <- bind_cols(df,fp)
  ### make footprints binary
  df <- df %>% dplyr::mutate(fp = case_when(footprints > 0 ~ TRUE, 
                                            footprints == 0 ~ FALSE, 
                                            T ~ NA))
  return(df)
}

add_footprints_5 <- function(df, variant.gr){
  footprints.gr <- read_delim("/mnt/sdb/ukbb-finemapping/data/annotations/meuleman/consensus_footprints_and_collapsed_motifs_hg38.bed.gz", delim = "\t",
                              col_names = c("seqname", "start", "end", "identifier", "mean_signal", "num_samples", "num_fps", "width", "summit_pos", "core_start", "core_end", "motif_clusters")) %>%
    dplyr::select(seqname, start, end, identifier, mean_signal, motif_clusters) %>% 
    makeGRangesFromDataFrame(., 
                             seqnames.field = "seqname", 
                             start.field =  "start", 
                             end.field = "end", 
                             keep.extra.columns = T) 
  ### lift over to hg19
  footprints.hg19.gr <- liftOver(footprints.gr, ch)
  rm(footprints.gr)
  footprints.hg19.gr <- unlist(footprints.hg19.gr)
  footprints.hg19.gr <- footprints.hg19.gr[!is.na(footprints.hg19.gr$identifier)]
  footprints.hg19.gr <- footprints.hg19.gr[!(duplicated(footprints.hg19.gr$identifier) | duplicated(footprints.hg19.gr$identifier, fromLast = TRUE)),]
  footprints.hg19.10.gr <- footprints.hg19.gr %>% anchor_center() %>% mutate(width=width + 10)
  names(footprints.hg19.10.gr) <- footprints.hg19.10.gr$identifier
  fp <- data.frame(countOverlaps(variant.gr, footprints.hg19.10.gr))
  names(fp) <- c("footprints_padded_5")
  df <- bind_cols(df,fp)
  ### make footprints binary
  df <- df %>% dplyr::mutate(fp_padded_5 = case_when(footprints_padded_5 > 0 ~ TRUE, 
                                                     footprints_padded_5 == 0 ~ FALSE, 
                                                     T ~ NA))
  return(df)
}


### find overlaps

### add in combination columns
def_intersecting_2 <- function(df){
  df <- df %>% ungroup() %>% 
    dplyr::mutate(footprints_mpra_any = case_when(fp == 1 & emVar_any == T ~T,
                                                  fp == 0 & emVar_any == F ~F,
                                                  fp == 0 & emVar_any == T ~F,
                                                  fp == 1 & emVar_any == F ~F,
                                                  T ~ NA ),
                  footprints_mpra_all = case_when(fp == 1 & emVar_all == T ~T,
                                                  fp == 0 & emVar_all == F ~F,
                                                  fp == 0 & emVar_all == T ~F,
                                                  fp == 1 & emVar_all == F ~F,
                                                  T ~ NA ),
                  no_CRE_yes_mpra_any = case_when(CRE == 1 & emVar_any == T ~F,
                                                  CRE == 0 & emVar_any == F ~F,
                                                  CRE == 0 & emVar_any == T ~T,
                                                  CRE == 1 & emVar_any == F ~F,
                                                  T ~ NA ),
                  no_CRE_yes_mpra_all = case_when(CRE == 1 & emVar_all == T ~F,
                                                  CRE == 0 & emVar_all == F ~F,
                                                  CRE == 0 & emVar_all == T ~T,
                                                  CRE == 1 & emVar_all == F ~F,
                                                  T ~ NA ),
                  no_footprints_yes_mpra_any = case_when(fp == 1 & emVar_any == T ~F,
                                                         fp == 0 & emVar_any == F ~F,
                                                         fp == 0 & emVar_any == T ~T,
                                                         fp == 1 & emVar_any == F ~F,
                                                         T ~ NA ),
                  no_footprints_yes_mpra_all = case_when(fp == 1 & emVar_all == T ~F,
                                                         fp == 0 & emVar_all == F ~F,
                                                         fp == 0 & emVar_all == T ~T,
                                                         fp == 1 & emVar_all == F ~F,
                                                         T ~ NA ),
                  yes_CRE_no_mpra_any = case_when(CRE == 1 & emVar_any == T ~F,
                                                  CRE == 0 & emVar_any == F ~F,
                                                  CRE == 0 & emVar_any == T ~F,
                                                  CRE == 1 & emVar_any == F ~T,
                                                  T ~ NA ),
                  no_CRE_yes_active_any = case_when(CRE == 1 & active_any == T ~F,
                                                    CRE == 0 & active_any == F ~F,
                                                    CRE == 0 & active_any == T ~T,
                                                    CRE == 1 & active_any == F ~F,
                                                    T ~ NA ),
                  promoters_mpra_any = case_when(promoter == 1 & emVar_any == T ~T,
                                                 promoter == 0 & emVar_any == F ~F,
                                                 promoter == 0 & emVar_any == T ~F,
                                                 promoter == 1 & emVar_any == F ~F,
                                                 T ~ NA ),
                  promoters_mpra_all = case_when(promoter == 1 & emVar_all == T ~T,
                                                 promoter == 0 & emVar_all == F ~F,
                                                 promoter == 0 & emVar_all == T ~F,
                                                 promoter == 1 & emVar_all == F ~F,
                                                 T ~ NA ))
  df <- df %>% group_by(variant)
  return(df)
}


#### annotate with ChIP
add_chip <- function(df,variant.gr){
  chip.gr <- readRDS("/mnt/sdb/ukbb-finemapping/data/chipatlas/ChIPAtlas.ChIP.rds")
  print('loaded chip in')
  #qval neg log10 FDR
  # merged peaks of ChIP-seq experiments
  ## can filter to one antigen and loop through, or can filter to a cell type class and loop through
  chip.gr.list <- split(chip.gr,as.factor(chip.gr$antigen))
  print('split into gr list')
  ##### for the entire chIp-seq df - annotate
  chip.df <- data.frame(lapply(1:length(chip.gr.list), function(x) {countOverlaps(variant.gr, chip.gr.list[[x]])})) ## and then filter to greater than or equal to 1
  print('found overlaps')
  ## this can be slow - take this ChIP GRanges object and turn into a list of granges per TF
  names(chip.df) <- paste("TF" ,names(chip.gr.list),sep=".")
  print('renamed')
  chip.names <- names(chip.df)
  print('renamed again')
  #names(chip.df) <- names(chip.gr.list)
  
  test_chip.df <- data.frame(lapply(names(chip.df), function(x) {ifelse(chip.df[[x]]==0,0,1)}))
  print('making some adjustment here')
  
  names(test_chip.df) <- names(chip.df)
  print('renaming')
  
  df <- bind_cols(df,test_chip.df)
  print('finished combining')
  return(df)
}


#### read in SEI
add_sei <- function(df, variant.gr){
  sei.names <- read_delim("/mnt/sdb/ukbb-finemapping/data/sei/resources/cnames.tsv", delim = "\t", col_names = T)
  sei.gr <-  import("/mnt/sdb/gtex_mpra/data/annotations/sei/sorted.hg38.tiling.bed.ipca_randomized_300.labels.merged.bed", format="bed")
  print('imported sei')
  
  
  sei.hg19.gr <- liftOver(sei.gr, ch)
  print('lifted sei over')
  #rm(sei.gr)
  sei.hg19.gr <- unlist(sei.hg19.gr)
  sei.hg19.gr <- sei.hg19.gr[!is.na(sei.hg19.gr$name)]
  sei.hg19.gr <- makeGRangesFromDataFrame(left_join(transform(as_tibble(sei.hg19.gr), name=as.numeric(name)), sei.names, by= c("name"="index")), keep.extra.columns = TRUE)
  #names(sei.hg19.gr) <- sei.hg19.gr$name.y
  sei.hg19.gr <- sei.hg19.gr[!is.na(sei.hg19.gr$name.y)]
  print('fixed up sei')
  
  sei.hg19.gr.list <- split(sei.hg19.gr,as.factor(sei.hg19.gr$name.y))
  print('split sei')
  
  ### get the overlaps with variant.gr (made earlier)
  sei.df <- data.frame(lapply(1:length(sei.hg19.gr.list), function(x) {countOverlaps(variant.gr, sei.hg19.gr.list[[x]])})) ## and then filter to greater than or equal to 1
  ## attach to mpra
  #names(sei.df) <- names(sei.hg19.gr.list)
  names(sei.df) <- paste("SEI",names(sei.hg19.gr.list) %>% str_replace_all('[ -]', '.') %>% str_replace_all('/','or'),sep='.')
  print('calculated overlaps and renamed')
  rm(sei.gr)
  rm(sei.hg19.gr)
  rm(sei.hg19.gr.list)
  df <- bind_cols(df,sei.df)
  return(df)
}


