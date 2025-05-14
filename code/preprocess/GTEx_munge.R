# Load libraries
library(tidyverse)
#library(Vennerable)
library(GenomicRanges)
library(BuenColors)
library(rtracklayer)
#library(regioneR)
library(reshape2)
#library(MatchIt)
library(ComplexHeatmap)
library(binom)
#library(meta)
setwd('/mnt/sdb/gwas_eqtl_mpra/')
source("code/utils.R")
options(stringsAsFactors = FALSE)
library("magrittr")


###############
# (1) Read in GTEx fine-mapping
###############

# Read in GTEx fine-mapping variant info
GTEx.info <- read_delim("../data/GTEx_49tissues.txt", delim = "\t", col_names = T)
GTEx.CS.cols <- read_delim("../release/GTEx_49tissues_release1.cols", delim = "\t", col_names = F)
GTEx.CS.df <- read_delim("../release/GTEx_49tissues_release1.tsv.gz", delim = "\t", col_names = GTEx.CS.cols$X1) %>%
  dplyr::mutate("position" = start) %>%
  dplyr::select(-start, -end) %>%
  dplyr::mutate(variant = gsub("_", ":", variant))

# Read in GTEx variant fine-mapping variant info (restricted to MPRA)
GTEx.test.df <- read_delim("../mpra_design/GTEx/GTEx_test_variants.df.txt.gz", delim = " ", col_names = T) %>%
  dplyr::mutate("position" = start) %>%
  dplyr::select(-start, -end) %>%
  dplyr::mutate(variant = gsub("_", ":", variant))

# Filter to variants we tested by MPRA
GTEx.CS.df <- GTEx.CS.df %>%
  dplyr::filter(variant %in% GTEx.test.df$variant)

# Add CS and type info to main dataset
GTEx.CS.df <- GTEx.CS.df %>%
  left_join(GTEx.test.df %>%
              dplyr::select(variant, method, tissue, gene, cs_id, cs_log10bf, cs_avg_r2, cs_min_r2, "cs_n" = n, type),
            by = c("variant", "method", "tissue", "gene", "cs_id"))

# Filter to just SuSiE for now
GTEx.test.df <- GTEx.test.df %>%
  dplyr::filter(method == "SUSIE")
GTEx.CS.df <- GTEx.CS.df %>%
  dplyr::filter(method == "SUSIE")

# Read in control variant data
ctrl_loc <- read_ctrl_loc()
ctrl_annot <- read_ctrl_annot()

###############
# (3) Annotate variants
###############

# Read in VEP annotations
GTEx.vep <- read_delim("../release/GTEX_49tissues_release1.vep.most_severe.tsv.gz", delim = "\t",
                       col_types = cols(.default = "c")) %>%
  distinct() %>%
  dplyr::rename("vep.most_severe" = most_severe,
                "vep.gene_most_severe" = gene_most_severe) %>%
  dplyr::select(variant, vep.most_severe, vep.gene_most_severe, lof)

# Add VEP annotations
GTEx.CS.df <- GTEx.CS.df %>%
  dplyr::mutate(variant = gsub("_", ":", variant)) %>%
  left_join(GTEx.vep,
            by = "variant")

# Read in general annotations
coding.gr <- import("../data/annotations/ldsc/Coding_UCSC.bed", format="bed")
UTR5.gr <- import("../data/annotations/ldsc/UTR_5_UCSC.bed", format="bed")
UTR3.gr <- import("../data/annotations/ldsc/UTR_3_UCSC.bed", format="bed")
conserved.gr <- import("../data/annotations/ldsc/Conserved_LindbladToh.bed", format="bed")
DHS_Roadmap.gr <- import("../data/annotations/ldsc/DHS_peaks_Trynka.bed", format="bed")
H3K27ac.gr <- import("../data/annotations/ldsc/H3K27ac_PGC2.bed", format="bed")
promoter.gr <- import("../data/annotations/ldsc/Promoter_UCSC.bed", format="bed")
intron.gr <- import("../data/annotations/ldsc/Intron_UCSC.bed", format="bed")
DHS.gr <- readRDS("../data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[1]]
K27_1.gr <- readRDS("../data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[2]]
K27_2.gr <- readRDS("../data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[3]]

# Read in specific chromatin accesibility annotations
DHS_Meuleman.gr <- readRDS("../data/annotations/meuleman/meuleman_en_gr.rds")
DHS_Domcke.gr <- readRDS("../data/annotations/domcke/domcke_en_gr.rds")
DHS_Corces.gr <- readRDS("../data/annotations/corces/corces_en_gr.rds")
DHS_Buen.gr <- readRDS("../data/annotations/buen/buen_en_gr.rds")
DHS_Calderon.gr <- readRDS("../data/annotations/calderon/calderon_en_gr.rds")
CA_DHS.gr <- readRDS("../data/annotations/ChIPAtlas.dnase.rds")
CA_DHS.q500.gr <- CA_DHS.gr[CA_DHS.gr$qval > 500] %>%
  GRsize(., shift = 150) %>%
  GenomicRanges::reduce()

# Create annotation GR list
annotations.gr.list <- list(coding.gr, UTR5.gr, UTR3.gr, conserved.gr, promoter.gr, intron.gr, DHS.gr, K27_1.gr, K27_2.gr, DHS_Roadmap.gr, DHS_Meuleman.gr, DHS_Domcke.gr, DHS_Corces.gr, DHS_Buen.gr, DHS_Calderon.gr, CA_DHS.q500.gr)
names(annotations.gr.list) <- c("coding", "UTR5", "UTR3", "conserved", "promoter", "intron", "DHS", "K27ac_1", "K27ac_2", "DHS_Roadmap", "DHS_Meuleman", "DHS_Domcke", "DHS_Corces", "DHS_Buen", "DHS_Calderon", "DHS_ChipAtlas")

# Annotate GTEx variants
variant.gr <- GTEx.CS.df %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "position", 
                           end.field = "position",
                           keep.extra.columns = T)
annot.df <- data.frame(lapply(1:length(annotations.gr.list), function(x) {countOverlaps(variant.gr, annotations.gr.list[[x]])}))
names(annot.df) <- names(annotations.gr.list)
GTEx.CS.df <- bind_cols(GTEx.CS.df, annot.df) %>%
  dplyr::mutate(CRE_Roadmap = ifelse(DHS_Roadmap >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_Meuleman = ifelse(DHS_Meuleman >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_Domcke = ifelse(DHS_Domcke >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_Corces = ifelse(DHS_Corces >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_Buen = ifelse(DHS_Buen >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_Calderon = ifelse(DHS_Calderon >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE_ChipAtlas = ifelse(DHS_ChipAtlas >= 1 & (K27ac_1 + K27ac_2 >= 1), 1, 0),
                CRE = ifelse(CRE_Meuleman + CRE_Domcke + CRE_Corces + CRE_Buen + CRE_Calderon + CRE_ChipAtlas >= 1, 1, 0)) %>%
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

###############
# (2) Read in MPRA results
###############

# Set up parameters
arrays <- c("OL41", "OL41B", "OL41_42", "OL42")
# For Method 1
fullfiles <- list.files(path = "../final_mpra_data/emvars/", pattern = "*paired*", full.names = T)
files <- list.files(path = "../final_mpra_data/emvars/", pattern = "*paired*")
# For Method 2
fullfiles <- list.files(path = "../final_mpra_data/emvars/", pattern = "*glm_20*", full.names = T)
files <- list.files(path = "../final_mpra_data/emvars/", pattern = "*glm_20*")
libraries <- gsub("_emVAR_.*", "", files) %>%
  sub("_([^_]*)$", "", .)
cell_types <- gsub("_emVAR_.*", "", files) %>%
  gsub(".*_", "", .)

# Read in MPRA data
mpra_input <- lapply(1:length(files),
                     function(x) {
                       read_mpra(fullfiles[x], cell_types[x], libraries[x])
                     }) %>%
  bind_rows() %>%
  dplyr::filter(library %in% arrays)

# Filter to best covered variants
mpra_input <- mpra_input %>%
  dplyr::mutate(cell_type = ifelse(cell_type == "SK.N.SH", "SKNSH", cell_type),
                cell_type = ifelse(cell_type == "HepG2", "HEPG2", cell_type)) %>%
  left_join(mpra_dnarna,
            by = c("variant", "cell_type", "library")) %>%
  dplyr::mutate(read_thresh = ifelse(mean_Plasmid_alt > 30 & mean_Plasmid_ref > 30 & mean_RNA_alt > 0 & mean_RNA_ref > 0, 1, 0)) %>%
  dplyr::group_by(ID, cell_type) %>%
  dplyr::mutate(best_library = ifelse((read_thresh == 1) & ((dplyr::n() == 1) | (sqrt(mean_Plasmid_alt) + sqrt(mean_Plasmid_ref)) == max(sqrt(mean_Plasmid_alt) + sqrt(mean_Plasmid_ref))), 1 , 0)) %>%
  ungroup()

# Read in MPRA metadata
mpra_meta <- data.table::fread("../final_mpra_data/attributes/OL41-42.attributes", sep = "\t") 
#names(mpra_meta) <- strsplit(names(mpra_meta)[1], "\\s+")[[1]]

# Add raw counts
files <- list.files(path = "../final_mpra_data/counts/", pattern = "*_counts.out")
libraries <- gsub("_20.*_counts.out", "", files)
keep <- libraries %in% arrays
files <- files[keep]
libraries <- libraries[keep]
counts <- lapply(1:length(files), 
                 function(x) {
                   cnames <- c("ID", colnames(read_delim(paste0("../final_mpra_data/counts/", files[x]), delim = "\t", n_max = 0)))
                   read_delim(paste0("../final_mpra_data/counts/", files[x]), delim = "\t", col_names = cnames, skip = 1) %>%   
                     dplyr::mutate(library = libraries[x]) %>%
                     pivot_longer(cols = -c(ID, library)) %>%
                     dplyr::mutate(name = gsub("_.*", "", name)) %>%
                     dplyr::group_by(ID, name, library) %>%
                     dplyr::summarize(mean = mean(value))
                 }) %>%
  bind_rows() %>%
  ungroup() %>%
  dplyr::mutate(ID = gsub("\\(|\\)", "", ID)) %>%
  separate_rows(ID, sep = ";")
counts <- counts %>%
  dplyr::mutate(name = ifelse(name == "plasmid", "Plasmid", toupper(name)))

# Filter to Plasmid counts
mpra_dna <- mpra_meta %>%
  dplyr::filter(window == "center",
                strand == "fwd",
                haplotype == "ref") %>%
  dplyr::select("ID", "variant" = SNP, allele) %>%
  distinct() %>%
  left_join(counts %>%
              dplyr::filter(name == "Plasmid") %>%
              distinct(),
            by = "ID") %>%
  dplyr::filter(!is.na(mean)) %>%
  dplyr::select(-ID) %>%
  pivot_wider(names_from = c("name", "allele"), 
              values_from = c("mean"), 
              names_prefix = c("mean_"))

# Filter to RNA counts  
mpra_rna <- mpra_meta %>%
  dplyr::filter(window == "center",
                strand == "fwd",
                haplotype == "ref") %>%
  dplyr::select("ID", "variant" = SNP, allele) %>%
  distinct() %>%
  left_join(counts %>%
              dplyr::filter(name != "Plasmid"),
            by = "ID") %>%
  dplyr::filter(!is.na(mean)) %>%
  dplyr::select(-ID) %>%
  pivot_wider(names_from = c("allele"), 
              values_from = c("mean"), 
              names_prefix = c("mean_RNA_"))

mpra_dnarna <- mpra_dna %>%
  left_join(mpra_rna,
            by = c("variant", "library")) %>%
  dplyr::rename("cell_type" = name) %>%
  dplyr::mutate(variant = paste0("chr", variant))

# Annotate with fine-mapping data
mpra_test <- GTEx.CS.df %>%
  left_join(mpra_input,
            by = "variant") %>%
  left_join(GTEx.info %>%
              dplyr::select(tissue, group),
            by = "tissue") %>%
  dplyr::mutate(type = ifelse(is.na(type), "other_test", type))

# Optional: remove variants not successfully measured by MPRA
mpra_test <- mpra_test %>%
  dplyr::filter(!is.na(ID))

# Annotate with control data
ctls <- bind_rows(ctrl_loc, ctrl_annot) %>%
  distinct() %>%
  dplyr::filter(indexVar %in% mpra_test$variant)
mpra_ctrls <- ctls %>%
  inner_join(mpra_input,
             by = "variant") %>%
  dplyr::mutate(cohort = "control")

# Combine into final dataset
mpra <- bind_rows(mpra_test, mpra_ctrls)
mpra <- mpra %>% 
  dplyr::filter(!is.na(cohort))
# Method 1
saveRDS(mpra, "../final_mpra_data/processed/gtex_mpra_paired_final20230117.rds")
gtex_mpra_df <- readRDS("../final_mpra_data/processed/gtex_mpra_paired_final20230117.rds")
