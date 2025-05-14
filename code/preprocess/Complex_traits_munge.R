# Load libraries
library(tidyverse)
library(GenomicRanges)
library(BuenColors)
library(rtracklayer)
library(reshape2)
library(ComplexHeatmap)
library(binom)
library(magrittr)
library(MESS)
setwd('/mnt/sdb/gwas_eqtl_mpra/')
source("code/utils.R")
options(stringsAsFactors = FALSE)

###############
# (1) Read in complex trait fine-mapping
###############

# Read in traits
array_layout <- read_delim("/mnt/sdb/gtex_mpra/mpra_design/ukbb_arrays.txt", delim = "\t", col_names = T)

# Read in BBJ fine-mapping
BBJ_df <- list.files(path = "/mnt/sdb/encode-finemapping/results/BBJ/combined", pattern = ".SuSiE.snp.gz", full.names = T) %>% 
  lapply(function(x) {read_delim(x, delim = " ")} %>% 
           dplyr::mutate(chromosome = as.character(chromosome),
                         variant = paste0("chr", chromosome, ":", position, ":", allele1, ":", allele2))) %>%
  bind_rows() %>%
  dplyr::filter(trait %in% array_layout$trait)

# Annotate BBJ with VEP
VEP <- read_delim("/mnt/sdb/encode-finemapping/data/annotations/VEP/bbj.ukbb.hg19.all.vep.tsv.gz", delim = "\t") %>%
  dplyr::mutate(variant = paste0("chr", v_hg19)) %>%
  dplyr::filter(variant %in% BBJ_df$variant)
BBJ_df <- BBJ_df %>%
  left_join(VEP %>%
              dplyr::select(variant, most_severe, gene_most_severe, impact, lof),
            by = "variant")
BBJ_df <- BBJ_df %>%
  dplyr::mutate(pip = prob,
                chisq_marginal = qchisq(p, df = 1, lower.tail = F),
                chromosome = paste0("chr", chromosome)) %>%
  dplyr::select(-flip, -p) %>%
  dplyr::rename(cohort = pop, beta_marginal = beta, se_marginal = se, susie.beta_posterior = susie_beta,  susie.sd_posterior = susie_sd,
                susie.pip = prob, vep.most_severe = most_severe, vep.gene_most_severe = gene_most_severe, vep.impact = impact,
                susie.cs_id = cs, vep.lof = lof)
saveRDS(BBJ_df, "/mnt/sdb/gtex_mpra/release/BBJ_processed.rds")

# Read in UKBB fine-mapping variant info
UKBB_df <- data.table::fread(cmd="cat /mnt/sdb/gtex_mpra/release/UKBB_96traits_release1.average.tsv.bgz | gunzip -dc", data.table=F) %>%
  as_tibble() %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2]))

# Add CS and type info to UKBB
# Filter to just SuSiE
UKBB_cs <- data.table::fread(cmd="cat /mnt/sdb/gtex_mpra/release/UKBB_96traits_release1_SuSiE_cred.tsv.gz | gunzip -dc", data.table=F) %>%
  as_tibble() %>%
  dplyr::mutate(id = paste(trait, region, cs_id, sep = ";"))
UKBB_df <- UKBB_df %>%
  dplyr::filter(!is.na(susie.pip)) %>%
  dplyr::mutate(pip = ifelse(is.na(pip), susie.pip, pip)) %>%
  dplyr::mutate(method = "SUSIE") %>%
  left_join(UKBB_cs %>%
              dplyr::rename("cs_n" = n),
            by = c("trait", "region", "susie.cs_id" = "cs_id"))
saveRDS(UKBB_df, "../release/UKBB_processed.rds")

# Combine UKBB and BBJ
# BBJ missing: minor_allele, cs_log10bf, cs_avg_r2, cs_min_r2, cs_n
BBJ_df <- readRDS("../release/BBJ_processed.rds")
UKBB_df <- readRDS("../release/UKBB_processed.rds")
traits_df <- bind_rows(UKBB_df, BBJ_df) %>%
  dplyr::mutate(method = "SUSIE") %>%
  dplyr::rename(minor_allele = minorallele, cs_id = susie.cs_id, beta_posterior = susie.beta_posterior, 
                sd_posterior = susie.sd_posterior, lof = vep.lof) %>%
  dplyr::select(chromosome, variant, allele1, allele2, minor_allele, cohort, method, trait, region, maf, 
                beta_marginal, se_marginal, chisq_marginal, pip, cs_id, beta_posterior, sd_posterior, position, cs_log10bf,
                cs_avg_r2, cs_min_r2, cs_n, vep.most_severe, vep.gene_most_severe, lof)

###############
# (3) Annotate variants
###############

# Change DHS size
GRsize <- function(gr, shift = 250) {
  midpt <- (start(gr) + floor(width(gr)/2))
  start(gr) <- as.integer(midpt - shift)
  end(gr) <- as.integer(midpt + shift)
  return(gr)
}

# Read in general annotations
coding.gr <- rtracklayer::import("../data/annotations/ldsc/Coding_UCSC.bed", format="bed")
UTR5.gr <- rtracklayer::import("../data/annotations/ldsc/UTR_5_UCSC.bed", format="bed")
UTR3.gr <- rtracklayer::import("../data/annotations/ldsc/UTR_3_UCSC.bed", format="bed")
conserved.gr <- rtracklayer::import("../data/annotations/ldsc/Conserved_LindbladToh.bed", format="bed")
DHS_Roadmap.gr <- rtracklayer::import("../data/annotations/ldsc/DHS_peaks_Trynka.bed", format="bed")
H3K27ac.gr <- rtracklayer::import("../data/annotations/ldsc/H3K27ac_PGC2.bed", format="bed")
promoter.gr <- rtracklayer::import("../data/annotations/ldsc/Promoter_UCSC.bed", format="bed")
intron.gr <- rtracklayer::import("../data/annotations/ldsc/Intron_UCSC.bed", format="bed")
DHS.gr <- readRDS("../data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[1]]
K27_1.gr <- readRDS("../data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[2]]
K27_2.gr <- readRDS("../data/annotations/DHSmerged_ROADMAPK27ac_CAK27ac.rds")[[3]]

# Read in specific chromatin accesibility annotations
DHS_Meuleman.gr <- readRDS("../data/annotations/meuleman/meuleman_en_gr.rds")
DHS_Domcke.gr <- readRDS("../data/annotations/domcke/domcke_en_gr.rds")
DHS_Corces.gr <- readRDS("../data/annotations/corces/corces_en_gr.rds")
DHS_Buen.gr <- readRDS("../data/annotations/buen/buen_en_gr.rds")
DHS_Calderon.gr <- readRDS("../data/annotations/calderon/calderon_en_gr.rds")
CA_DHS.gr <- readRDS("../../ukbb-finemapping/data/chipatlas/ChIPAtlas.dnase.rds")
CA_DHS.q500.gr <- CA_DHS.gr[CA_DHS.gr$qval > 500] %>%
  GRsize(., shift = 150) %>%
  GenomicRanges::reduce()

# Create annotation GR list
annotations.gr.list <- list(coding.gr, UTR5.gr, UTR3.gr, conserved.gr, promoter.gr, intron.gr, DHS.gr, K27_1.gr, K27_2.gr, DHS_Roadmap.gr, DHS_Meuleman.gr, DHS_Domcke.gr, DHS_Corces.gr, DHS_Buen.gr, DHS_Calderon.gr, CA_DHS.q500.gr)
names(annotations.gr.list) <- c("coding", "UTR5", "UTR3", "conserved", "promoter", "intron", "DHS", "K27ac_1", "K27ac_2", "DHS_Roadmap", "DHS_Meuleman", "DHS_Domcke", "DHS_Corces", "DHS_Buen", "DHS_Calderon", "DHS_ChipAtlas")

# Annotate UKBB variants
variant.gr <- traits_df %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "position", 
                           end.field = "position",
                           keep.extra.columns = T)
annot.df <- data.frame(lapply(1:length(annotations.gr.list), function(x) {countOverlaps(variant.gr, annotations.gr.list[[x]])}))
names(annot.df) <- names(annotations.gr.list)
traits_df <- bind_cols(traits_df, annot.df) %>%
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

saveRDS(traits_df, "../release/ukbb_bbj.rds")
traits_df <- readRDS("../release/ukbb_bbj.rds")

###############
# (2) Read in MPRA variant information
###############
cell_types <- c("K562", "HepG2", "GM12878", "SKNSH", "A549")
array_layout <- read_delim("../mpra_design/ukbb_arrays.txt", delim = "\t", col_names = T)

# Define function for reading in variants
# Read in test variants
read_ukbb_bbj_test <- function(trait_names) {
  pmap_dfr(
    crossing(trait_name = trait_names),
    function(trait_name) {
      df1 <- read_delim(paste0("../../encode-finemapping/mpra_design/", trait_name, "/", trait_name, "_test_variants.df.txt.gz"), delim = " ", col_names = T) %>%
        dplyr::mutate("variant" = v)
      df2 <- df1 %>%
        dplyr::filter(prob > 0.1) %>%
        dplyr::mutate(type = "PIP10")
      df3 <- df1 %>%
        dplyr::filter(cs > 0) %>%
        dplyr::mutate(type = "CS")
      bind_rows(df2, df3) %>%
        dplyr::select("cohort" = pop, trait, region, variant, "cs_id" = cs, method, type)
    }
  )
}

# Read in null GWAS control variants
read_ctrl_null <- function(trait_names) {
  pmap_dfr(
    crossing(trait_name = trait_names),
    function(trait_name) {
      read_delim(paste0("../../encode-finemapping/mpra_design/", trait_name, "/", trait_name, "_nullGWAScontrol_variants.df.txt.gz"), delim = " ", col_names = T) %>%
        dplyr::mutate(match_type = ifelse(prob > 0.1, "null_PIP10", "null_CS")) %>%
        dplyr::select("variant" = v, match_type) %>%
        dplyr::mutate(trait = trait_name)
    }
  )
}

# Read in location control variants 
read_ctrl_loc <- function(trait_names) {
  pmap_dfr(
    crossing(trait_name = trait_names),
    function(trait_name) {
      read_delim(paste0("../../encode-finemapping/mpra_design/", trait_name, "/", trait_name, "_locationcontrol_variants.df.txt.gz"), delim = " ", col_names = T) %>%
        dplyr::select("variant" = v, "match_type" = match, "match_variant" = highPPvar_v) %>%
        dplyr::mutate(match_type = ifelse(match_type == "PP10", "loc_PIP10", "loc_CS")) %>%
        dplyr::mutate(trait = trait_name)
    }
  )
}

# Read in annoation control variants
read_ctrl_annot <- function(trait_names) {
  pmap_dfr(
    crossing(trait_name = trait_names),
    function(trait_name) {
      read_delim(paste0("../../encode-finemapping/mpra_design/", trait_name, "/", trait_name, "_annotcontrol_variants.df.txt.gz"), delim = " ", col_names = T) %>%
        dplyr::select("variant" = v, "match_type" = match) %>%
        dplyr::mutate(match_type = ifelse(match_type == "PP10", "annot_PIP10", "annot_CS")) %>%
        dplyr::mutate(trait = trait_name)
    }
  )
}

# Read in test variant data
test_df <- read_ukbb_bbj_test(array_layout$trait) %>%
  dplyr::filter(method == "SUSIE") %>%
  dplyr::select(-method)
traits_df <- traits_df %>%
  left_join(test_df,
            by = c("cohort", "trait", "region", "cs_id", "variant"))

# Read in control variant data
ctrl_null <- read_ctrl_null(array_layout$trait)
ctrl_loc <- read_ctrl_loc(array_layout$trait)
ctrl_annot <- read_ctrl_annot(array_layout$trait)

###############
# (2) Read in MPRA results
###############

# Read in MPRA metadata
mpra_meta <- list.files(path = "../final_mpra_data/attributes/", pattern = "UKBB.attributes", full.names = T) %>% 
  lapply(function(x) {read_delim(x, delim = "\t")} %>% 
           dplyr::mutate(chr = as.character(chr))) %>%
  bind_rows()

# Add raw counts
files <- list.files(path = "../final_mpra_data/counts/", pattern = "*_counts.out")
libraries <- gsub("_.*_counts.out", "", files)
keep <- libraries %in% unique(array_layout$library)
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

# Set up parameters
# For Method 1
fullfiles <- list.files(path = "../final_mpra_data/emvars/", pattern = "*paired*", full.names = T)
files <- list.files(path = "../final_mpra_data/emvars/", pattern = "*paired*")
# For Method 2
fullfiles <- list.files(path = "../final_mpra_data/emvars/", pattern = "*glm_20*", full.names = T)
files <- list.files(path = "../final_mpra_data/emvars/", pattern = "*glm_20*")
libraries <- gsub("_.*", "", files)
cell_types <- gsub("_emVAR_.*", "", files) %>%
  gsub(".*_", "", .)

# Read in MPRA data
mpra_input <- lapply(1:length(files),
                     function(x) {
                       read_mpra(fullfiles[x], cell_types[x], libraries[x])
                     }) %>%
  bind_rows() %>%
  dplyr::filter(library %in% unique(array_layout$library))

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

# Annotate with fine-mapping data
mpra_test <- traits_df %>%
  left_join(mpra_input,
            by = "variant") %>%
  dplyr::mutate(type = ifelse(is.na(type), "other_test", type))

# Optional: remove variants not successfully measured by MPRA
mpra_test <- mpra_test %>%
  dplyr::filter(!is.na(ID))

# Annotate with control data
ctrl_loc <- ctrl_loc %>%
  distinct() %>%
  dplyr::rename(indexVar = match_variant) %>%
  dplyr::filter(indexVar %in% mpra_test$variant)
ctrls <- bind_rows(ctrl_loc, ctrl_annot, ctrl_null) %>%
  dplyr::rename(type = match_type) %>%
  distinct()
mpra_ctrls <- ctrls %>%
  inner_join(mpra_input,
             by = c("variant")) %>%
  dplyr::mutate(cohort = "control")

# Combine into final dataset
mpra <- bind_rows(mpra_test, mpra_ctrls)
mpra <- mpra %>% 
  dplyr::filter(!is.na(cohort))
# Method 1
saveRDS(mpra, "../final_mpra_data/processed/traits_mpra_paired_final20230117.rds")
traits_mpra_df <- readRDS("../final_mpra_data/processed/traits_mpra_paired_final20230117.rds")
# Method 2
saveRDS(mpra, "../final_mpra_data/processed/traits_mpra_unpaired_final20230117.rds")
traits_mpra_df <- readRDS("../final_mpra_data/processed/traits_mpra_unpaired_final20230117.rds")
