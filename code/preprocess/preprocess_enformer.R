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
source("code/utils.R")
source("code/theme.R")
options(stringsAsFactors = FALSE)

# Read in TF meta data
targets <- read_delim("https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt", delim = "\t") %>%
  as_tibble()

# Get Enformer files (two separate runs)
files_gtex <- list.files("../data/annotations/enformer/raw/", pattern = "*.tsv.gz") %>% 
  str_remove_all("_refscore") %>% 
  unique()
files_ukbb <- list.files("../../ukbb-finemapping/data/enformer/", pattern = "*.vcf.gz") %>% 
  str_remove_all("_refscore") %>% 
  unique()

# Get variants
keep_vars <- mpra %>%
  dplyr::filter(CRE > 0) %>%
  distinct(variant, variant_hg38)

# Get control variants for normalization
highpip_df <- mpra %>%
  filter(pip > 0.9, CRE > 0, cs_id > 0) %>%
  group_by(cs_id) %>%
  filter(any(consequence %ni% c("synonymous","missense","LoF"))) %>%
  ungroup() %>% 
  mutate(category = case_when(pchisq(chisq_marginal, 1, log.p = TRUE, lower.tail = F) / -log(10) > -log10(5 * 10^-8) & type %in% c("CS", "PIP10") ~ "HighPIP",
                              abs(z) > qnorm(1 - 5*10^-8) & type %in% c("49tissue_PIP50", "3tissue_PIP10", "3tissue_CS") ~ "HighPIP")) %>%
  filter(!is.na(category)) %>%
  dplyr::group_by(variant) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::ungroup() %>%
  mutate(category = case_when(emVar_any == T ~ "HighPIP_emVar",
                              T ~ "HighPIP_only"))
lowpip_df <- mpra %>%
  filter(!is.na(pip)) %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(pip = max(pip, na.rm = T)) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::ungroup() %>%
  filter(pip < 0.05, CRE > 0, cs_id > 0) %>%
  group_by(cs_id) %>%
  filter(any(consequence %ni% c("synonymous","missense","LoF"))) %>%
  ungroup() %>% 
  mutate(category = case_when(type %in% c("CS", "PIP10") ~ "LowPIP",
                              type %in% c("49tissue_PIP50", "3tissue_PIP10", "3tissue_CS") ~ "LowPIP")) %>%
  filter(!is.na(category)) %>%
  filter(variant %ni% highpip_df$variant) %>%
  dplyr::group_by(variant) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::ungroup() %>%
  mutate(category = case_when(emVar_any == T ~ "LowPIP_emVar",
                              T ~ "LowPIP_only"))
mech_df <- bind_rows(highpip_df, lowpip_df)
rm(highpip_df); rm(lowpip_df); gc()
mech_df <- mech_df %>%
  dplyr::filter(nchar(allele1) == 1, nchar(allele2) == 1)

# Write out mechanism filtered data
mech_df %>% 
  vroom::vroom_write("/mnt/sdb/gtex_mpra/data/mech_df.txt.gz")

# Get variants to match up to Enformer labels
vars_hg19 <- mech_df %>%
  distinct(variant) %>%
  pull(1)
vars_hg38 <- mech_df %>%
  distinct(variant_hg38) %>%
  pull(1)
vars_hg19 <- c(vars_hg19, mpra %>% 
                 dplyr::filter(variant %in% satmut_emVar_df$variant) %>% 
                 distinct(variant) %>%
                 na.omit() %>% 
                 .$variant)
vars_hg38 <- c(vars_hg38, mpra %>% 
                 dplyr::filter(variant %in% satmut_emVar_df$variant) %>% 
                 distinct(variant_hg38) %>%
                 .$variant_hg38) %>%
  na.omit()
ctrl_vars_hg19 <- mech_df %>%
  dplyr::filter(category == "LowPIP_only") %>%
  distinct(variant) %>%
  pull(1)
vars_gtex_df <- bind_rows(mech_df %>% 
                       distinct(variant, variant_hg38) %>%
                       na.omit(),
                     mpra %>% 
                       dplyr::filter(variant %in% satmut_emVar_df$variant) %>% 
                       distinct(variant, variant_hg38) %>% 
                       na.omit()) %>%
  distinct() 

# Remove duplicates and match indices
enf_gtex_mat <- data.table::fread(paste0("../data/annotations/enformer/raw/", files_gtex[1]), 
                             data.table = F, header = F, nrows = 1) %>%
  as_tibble()
enf_ukbb_mat <- data.table::fread(paste0("../../ukbb-finemapping/data/enformer/", files_ukbb[1]), 
                                  data.table = F, header = F, nrows = 1) %>%
  as_tibble()
enf_gtex_mat <- enf_gtex_mat %>% unlist() %>% as.vector() %>% gsub("^\\d+_", "", .)
enf_ukbb_mat <- enf_ukbb_mat %>% unlist() %>% as.vector()
enf_ukbb_indx <- match(unique(enf_ukbb_mat), enf_ukbb_mat)
enf_gtex_indx <- match(unique(enf_gtex_mat), enf_gtex_mat)
table(enf_ukbb_mat[enf_ukbb_indx] == enf_gtex_mat[enf_gtex_indx])

# Get proper names
ids_df <- bind_cols("id" = enf_gtex_mat[6:5318], "name" = targets$description) %>%
  group_by(id) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup()

# Read in ChIP only for mech variants
data.table::setDTthreads(24)
tmp1_gtex <- lapply(files_gtex, function(x) {
  enf.mat <- data.table::fread(paste0("../data/annotations/enformer/raw/", x), select = enf_gtex_indx, data.table = T) %>% 
    tidytable::as_tidytable() %>% 
    tidytable::filter(id %in% vars_hg38) %>%
    tidytable::select(id, contains("CHIP")) %>%
    tidytable::select(!contains("CHIP:H2") & !contains("CHIP:H3") & !contains("CHIP:H4")) %>%
    tidytable::select(!contains("batch") & !contains("abcam") & !contains("active"))
}) %>%
  tidytable::bind_rows()
gc()
tmp1_ukbb <- lapply(files_ukbb, function(x) {
  enf.mat <- data.table::fread(paste0("../../ukbb-finemapping/data/enformer/", x), select = enf_ukbb_indx, data.table = T) %>% 
    tidytable::as_tidytable() %>% 
    tidytable::filter(id %in% vars_hg19) %>%
    tidytable::select(id, contains("CHIP")) %>%
    tidytable::select(!contains("CHIP:H2") & !contains("CHIP:H3") & !contains("CHIP:H4")) %>%
    tidytable::select(!contains("batch") & !contains("abcam") & !contains("active"))
}) %>%
  tidytable::bind_rows()
gc()

# Munge outputs
tmp1_gtex <- tmp1_gtex %>%
  na.omit() %>% 
  tidytable::rename("variant_hg38" = id) %>%
  tidytable::inner_join(vars_gtex_df) %>%
  tidytable::mutate(variant_hg38 = variant) %>%
  tidytable::select(-variant) %>%
  tidytable::rename("variant" = variant_hg38)
tmp1_ukbb <- tmp1_ukbb %>%
  tidytable::rename("variant" = id) 
names(tmp1_gtex) <- gsub("^\\d+_", "", colnames(tmp1_gtex))

# Get mean when variant tested in both sets
tmp1 <- tidytable::bind_rows(tmp1_ukbb, tmp1_gtex) 
tmp1 <- tmp1 %>%
  tidytable::group_by(variant) %>%
  tidytable::summarize(tidytable::across(where(is.numeric), ~ mean(.x, na.rm = T)))
gc()

# Get mean and variance to normalize
tmp1_ctrl <- tmp1 %>%
  dplyr::filter(variant %in% ctrl_vars_hg19)
mu <- Matrix::colMeans(tmp1_ctrl[,-1])
tmp1_ctrl_mat <- tmp1_ctrl[,-1] %>% as.matrix()
sd <- matrixStats::colVars(tmp1_ctrl_mat) %>% sqrt()
#tmp1_ctrl_mat <- t((t(tmp1_ctrl_mat) - mu) / sd)
#var(tmp1_ctrl_mat[,2])
saveRDS(bind_cols("mu" = mu, "sd" = sd), "../data/annotations/enformer/norm_constants_CHIP.RDS")

# Normalize and save
tmp1_mat <- tmp1[,-1] %>% as.matrix()
tmp1_mat <- t((t(tmp1_mat) - mu) / sd)
tmp1 <- tidytable::bind_cols("variant" = tmp1$variant, tmp1_mat %>% tidytable::as_tidytable())
vroom::vroom_write(tmp1, "/mnt/sdb/gtex_mpra/data/annotations/enformer/enformer_CHIP_mech_satmut.txt.gz")
gc()

# Read in accessible chromatin only for mech variants
data.table::setDTthreads(24)
tmp1_gtex <- lapply(files_gtex, function(x) {
  enf.mat <- data.table::fread(paste0("../data/annotations/enformer/raw/", x), select = enf_gtex_indx, data.table = T) %>% 
    tidytable::as_tidytable() %>% 
    tidytable::filter(id %in% vars_hg38) %>%
    tidytable::select(id, contains("DNASE") | contains("ATAC"))
}) %>%
  tidytable::bind_rows()
gc()
tmp1_ukbb <- lapply(files_ukbb, function(x) {
  enf.mat <- data.table::fread(paste0("../../ukbb-finemapping/data/enformer/", x), select = enf_ukbb_indx, data.table = T) %>% 
    tidytable::as_tidytable() %>% 
    tidytable::filter(id %in% vars_hg19) %>%
    tidytable::select(id, contains("DNASE") | contains("ATAC"))
}) %>%
  tidytable::bind_rows()
gc()

# Munge outputs
tmp1_gtex <- tmp1_gtex %>%
  na.omit() %>% 
  tidytable::rename("variant_hg38" = id) %>%
  tidytable::inner_join(vars_gtex_df) %>%
  tidytable::mutate(variant_hg38 = variant) %>%
  tidytable::select(-variant) %>%
  tidytable::rename("variant" = variant_hg38)
tmp1_ukbb <- tmp1_ukbb %>%
  tidytable::rename("variant" = id) 
names(tmp1_gtex) <- gsub("^\\d+_", "", colnames(tmp1_gtex))

# Get mean when variant tested in both sets
tmp1 <- tidytable::bind_rows(tmp1_ukbb, tmp1_gtex) 
tmp1 <- tmp1 %>%
  tidytable::group_by(variant) %>%
  tidytable::summarize(across(where(is.numeric), ~ mean(.x)))
gc()

# Get mean and variance to normalize
tmp1_ctrl <- tmp1 %>%
  dplyr::filter(variant %in% ctrl_vars_hg19)
mu <- Matrix::colMeans(tmp1_ctrl[,-1])
tmp1_ctrl_mat <- tmp1_ctrl[,-1] %>% as.matrix()
sd <- matrixStats::colVars(tmp1_ctrl_mat) %>% sqrt()
tmp1_ctrl_mat <- t((t(tmp1_ctrl_mat) - mu) / sd)
var(tmp1_ctrl_mat[,2])
saveRDS(bind_cols("mu" = mu, "sd" = sd), "../data/annotations/enformer/norm_constants_DHS.RDS")

# Normalize and save
tmp1_mat <- tmp1[,-1] %>% as.matrix()
tmp1_mat <- t((t(tmp1_mat) - mu) / sd)
tmp1 <- tidytable::bind_cols("variant" = tmp1$variant, tmp1_mat %>% tidytable::as_tidytable())
vroom::vroom_write(tmp1, "/mnt/sdb/gtex_mpra/data/annotations/enformer/enformer_DHS_mech_satmut.txt.gz")
gc()

# Read in
enf_chip <- vroom::vroom("/mnt/sdb/gtex_mpra/data/annotations/enformer/enformer_CHIP_mech_satmut.txt.gz", delim = "\t", col_names = T) %>%
  tidytable::as_tidytable()
enf_dhs <- vroom::vroom("/mnt/sdb/gtex_mpra/data/annotations/enformer/enformer_DHS_mech_satmut.txt.gz", delim = "\t", col_names = T) %>%
  tidytable::as_tidytable()

# Fix missing SatMut variants
enf_chip_satmut <- tmp1_gtex %>%
  tidytable::filter(variant %in% satmut_emVar_df$variant) %>%
  tidytable::filter(variant %ni% enf_chip$variant)%>% 
  tidytable::as_tidytable() 
enf_dhs_satmut <- tmp2_gtex %>%
  tidytable::filter(variant %in% satmut_emVar_df$variant) %>%
  tidytable::filter(variant %ni% enf_dhs$variant)
saveRDS(list(enf_chip_satmut, enf_dhs_satmut), "/mnt/sdb/gtex_mpra/data/annotations/enformer/enformer_satmut_tmp.rds")
enf_chip <- tidytable::bind_rows(enf_chip, enf_chip_satmut)
enf_dhs <- tidytable::bind_rows(enf_dhs, enf_dhs_satmut)

# Get SS thresholds
ss_chip <- matrixStats::rowSums2(enf_chip[,-1]^2 %>% as.matrix())
ss_dhs <- matrixStats::rowSums2(enf_dhs[,-1]^2 %>% as.matrix())

# Get enformer thresholds
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

# SatMut
enf_chip %>%
  tidytable::filter(variant %in% satmut_emVar_df$variant)
satmut_emVar_df <- satmut_emVar_df %>%
  left_join(ss_df) %>%
  dplyr::mutate(enf_chip = ifelse(ss_chip > ss_chip_thresh, 1, 0),
                enf_dhs = ifelse(ss_dhs > ss_dhs_thresh, 1, 0)) %>%
  dplyr::select(-ss_chip, -ss_dhs) %>%
  left_join(enf_chip_thresh_df %>%
              tidytable::filter(variant %in% satmut_emVar_df$variant, call == T) %>%
              left_join(ids_df %>%
                          dplyr::rename("desc" = "name"),
                        by = c("name" = "id")) %>%
              tidytable::mutate(desc = gsub("CHIP:", "", desc),
                                desc = gsub(":.*", "", desc),
                                desc = gsub(".*FLAG-", "", desc),
                                desc = gsub(".*eGFP-", "", desc)) %>%
              tidytable::group_by(variant) %>% 
              tidytable::arrange(desc) %>% 
              tidytable::summarize(enf_chip_ids = paste0(desc, collapse = ";"))) %>%
  dplyr::mutate(enf_chip_ids = ifelse(enf_chip == T, enf_chip_ids, NA))
tmp %>%
  dplyr::select(variant, enf_chip, enf_chip_ids) %>% 
  dplyr::count(enf_chip)

# Add back in
mech_df <- mech_df %>%
  left_join(ss_df) %>%
  dplyr::mutate(enf_chip = ifelse(ss_chip > ss_chip_thresh, 1, 0),
                enf_dhs = ifelse(ss_dhs > ss_dhs_thresh, 1, 0)) 

### Not run after this
# Individual feautres passing thresholds
enf_chip_thresh_df <- enf_chip %>%
  tidytable::inner_join(mech_df %>%
                          distinct(variant, category) %>%
                          bind_rows(bind_cols("variant" = satmut_emVar_df$variant, category = "satmut"))) %>%
  #tidytable::mutate(across(where(is.numeric), ~ abs(.x))) %>%
  tidytable::pivot_longer(!c(variant, category)) %>%
  tidytable::group_by(name) %>%
  tidytable::mutate(thresh = quantile(abs(value[category == "LowPIP_only"]), 0.95),
                    call = value > abs(thresh))
enf_chip_thresh_df %>%
  tidytable::filter(value > thresh) %>%
  tidytable::select(variant, name)

enf_dhs_thresh_df <- enf_dhs %>%
  tidytable::inner_join(mech_df %>%
                          distinct(variant, category) %>%
                          bind_rows(bind_cols("variant" = satmut_emVar_df$variant, category = "satmut"))) %>%
  #tidytable::mutate(across(where(is.numeric), ~ abs(.x))) %>%
  tidytable::pivot_longer(!c(variant, category)) %>%
  tidytable::group_by(name) %>%
  tidytable::mutate(thresh = quantile(value[category == "LowPIP_only"], 0.95)) %>%
  tidytable::filter(value > thresh) %>%
  tidytable::select(variant, name)

enf_chip_var <- mech_df %>% dplyr::filter(enf_chip == T) %>% .$variant
enf_chip_thresh_ss_df <- enf_chip_thresh_df %>%
  tidytable::filter(variant %in% ctrl_vars_hg19) %>%
  tidytable::ungroup() 
#enf_chip_thresh_ss_df %>% 
enf_chip_thresh_plot_df <- enf_chip_thresh_df %>%
  tidytable::filter(category == "HighPIP_emVar") %>%
  #tidytable::mutate(denom = length(variant)) %>%
  tidytable::group_by(name, call) %>%
  tidytable::count() %>% 
  tidytable::group_by(name) %>%
  tidytable::mutate(denom = sum(n),
                       prop = n / denom) %>%
  tidytable::filter(call == T) %>%
  tidytable::group_by(name) %>%
  tidytable::mutate(lower = binom.test(n, denom, 0.05)$conf.int[1],
                    upper = binom.test(n, denom, 0.05)$conf.int[2],
                    pval = binom.test(n, denom, 0.05)$p.value) %>%
  left_join(ids_df %>%
              dplyr::rename("desc" = "name"),
            by = c("name" = "id")) %>%
  tidytable::mutate(desc = gsub("CHIP:", "", desc),
                    desc = gsub(":.*", "", desc),
                    desc = gsub(".*FLAG-", "", desc),
                    desc = gsub(".*eGFP-", "", desc)) %>%
  arrange(-n)



