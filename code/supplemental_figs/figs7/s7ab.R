# Load libraries
library(tidyverse)
library(vroom)
library(BuenColors)
library(binom)
library(cowplot)
library(patchwork)
library(GenomicRanges)
library(rtracklayer)
library(ggpubr)
source("code/utils.R")

# Read in processed MPRA results and filter to GTEx
mpra <- vroom("/mnt/sdb/gwas_eqtl_mpra/data/preprocess/core_mpra.txt.gz") %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0)) %>%
  dplyr::filter(cohort=='GTEx')

# Get positive set variants, don't filter down yet to a single row
positives_df <- mpra %>%
  dplyr::filter(pip > 0.9, 
                abs(z) > qnorm(1 - 5 * 10^-8), 
                type %in% c("49tissue_PIP50", "3tissue_PIP10", "3tissue_CS"), 
                consequence %ni% c("synonymous","missense","LoF")) %>% 
  dplyr::group_by(variant) %>%
  dplyr::mutate(active_any = ifelse(any(active==T),T, F)) %>%
  dplyr::mutate(pip = max(pip, na.rm = T)) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  dplyr::mutate(Skew_logPadj = ifelse(is.na(Skew_logPadj), 0, Skew_logPadj)) %>%
  arrange(-active,-Skew_logPadj) %>% 
  ungroup() %>%
  dplyr::mutate(causal = T) %>%
  distinct(variant, cell_type, active, log2FC, logPadj_BF, log2Skew, Skew_logPadj, causal)

# Get negative set variants
# Not enough variants to match if PIP < 0.01
negatives_df <- mpra %>%
  group_by(variant) %>%
  dplyr::mutate(active_any = ifelse(any(active==T),T, F)) %>%
  dplyr::mutate(pip = max(pip, na.rm = T)) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  dplyr::mutate(Skew_logPadj = ifelse(is.na(Skew_logPadj), 0, Skew_logPadj)) %>%
  arrange(-active,-Skew_logPadj) %>% 
  ungroup() %>%
  dplyr::filter(pip < 0.02, 
                type == '3tissue_CS', 
                consequence %ni% c("synonymous","missense","LoF"), 
                variant %ni% positives_df$variant) %>% 
  ungroup() %>%
  dplyr::mutate(causal = F) %>%
  distinct(variant, cell_type, active, log2FC, logPadj_BF, log2Skew, Skew_logPadj, causal)

# bind together
mpra_mc <- bind_rows(positives_df, negatives_df)

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
  dplyr::filter(cell_type %in% c("K562", "HepG2", "HCT116", "SK-N-DZ", "SK-N-MC")) %>%
  dplyr::mutate(cell_type = ifelse(cell_type == "HepG2", "HEPG2", cell_type),
                cell_type = ifelse(cell_type %in% c("SK-N-DZ", "SK-N-MC"), "SKNSH", cell_type))



# join with cell-type specific e2gs and filter down to only cell type specific e2g
mpra.overall <- mpra_mc %>% 
  left_join(e2g_filt_cts_df %>%
              distinct(variant, cell_type) %>%
              dplyr::mutate(is.e2g = 1),
            by = c('variant','cell_type')) %>%
  replace(is.na(.),0) 


prset_df <- mpra.overall %>%
  dplyr::mutate(is.e2g.active = ifelse(active + is.e2g == 2, T, F)) %>%
  group_by(variant) %>%
  arrange(-is.e2g.active,-Skew_logPadj) %>%
  filter(row_number() == 1) %>%
  ungroup()


print(prset_df %>% distinct(variant, causal) %>% count(causal))

# gridsearch
it_p <- rep(c(0.001,0.01,0.05,0.1,0.15,0.2,0.25,0.5,0.75,1),10)
it_a <- c(rep(c(0),10),rep(c(0.25),10),rep(c(0.5),10),rep(c(0.75),10),rep(c(1),10),rep(c(1.25),10),rep(c(1.5),10),rep(c(1.75),10),rep(c(2),10),rep(c(3),10))

prec_rec_emvar <- function(pthresh,athresh, df) {
  print(paste(pthresh, athresh, sep = ', '))
  # Get 2x2 table
  out_df <- df %>%
    # use below for TARV analysis
    dplyr::mutate(emVar_test = if_else(is.e2g & active & Skew_logPadj >= -log10(pthresh) & !is.na(Skew_logPadj) &  abs(log2Skew) >= athresh, T, F)) %>%
    dplyr::select(emVar_test, causal) %>%
    dplyr::mutate(annot = case_when(emVar_test == T ~ 1,
                                    emVar_test == F ~ 0)) %>%
    dplyr::filter(!is.na(emVar_test)) %>%
    dplyr::count(causal, annot)
  
  # Compute precision, recall, and 95% CIs
  out_df <- out_df %>%
    dplyr::summarize(TP = sum(n[causal == T & annot == T]),
                     FP = sum(n[causal == F & annot == T]),
                     FN = sum(n[causal == T & annot == F]),
                     TN = sum(n[causal == F & annot == F]),
                     total_pos = sum(n[annot == T]),
                     total = sum(n)) 

  # if there are more negatives:
  if (dim(df %>% dplyr::filter(causal == T))[1] < dim(df %>% dplyr::filter(causal == F))[1]){
    out_df <- out_df  %>%
      dplyr::mutate(FPR = FP / (FP + TN),
                    N = TP + FN,
                    FP_new = FPR * N,
                    TN_new = N - FP_new) %>%
      dplyr::mutate(FP = FP_new,
                    TN = TN_new) 
    # else if there are more positives
  } else if (dim(df %>% dplyr::filter(causal == T))[1] > dim(df %>% dplyr::filter(causal == F))[1]){
    out_df <- out_df  %>%
      dplyr::mutate(TPR = TP / (TP + FN),
                    N = TN + FP,
                    TP_new = TPR * N,
                    FN_new = N - TP_new) %>%
      dplyr::mutate(TP = TP_new,
                    FN = FN_new) 
  }

  out_df <- out_df %>%
    dplyr::mutate(precision = TP / (TP + FP),
                  recall = TP / (TP + FN),
                  prec_upper = binom.confint(TP, (TP + FP), method = "wilson")$upper,
                  prec_lower = binom.confint(TP, (TP + FP), method = "wilson")$lower,
                  rec_upper = binom.confint(TP, (TP + FN), method = "wilson")$upper,
                  rec_lower = binom.confint(TP, (TP + FN), method = "wilson")$lower) %>%
    dplyr::mutate(FDR = pthresh,
                  Skew = athresh)
  return(out_df)
}

data2 <- mapply(prec_rec_emvar,it_p,it_a, MoreArgs=list(df = prset_df))
tmp2 <- t(data2) %>% as.data.frame()
tmp2[is.na(tmp2)] <- 0

tmp2
# prepare for plotting
tmp2 <- transform(tmp2, FDR = as.character(FDR))%>% 
  transform(Skew = as.numeric(Skew)) %>%
  transform(precision = as.numeric(precision)) %>%
  transform(recall = as.numeric(recall))  %>%
  transform(TP = as.numeric(TP)) %>%
  transform(prec_upper = as.numeric(prec_upper)) %>%
  transform(prec_lower = as.numeric(prec_lower))

tmp2$FDR <- factor(tmp2$FDR, levels = c("0.001","0.01","0.05","0.1","0.15","0.2","0.25","0.5","0.75","1"))
tmp2$Skew <- factor(tmp2$Skew, levels = c("0","0.25","0.5","0.75","1","1.25","1.5","1.75","2","3"))

# plot
ph <- ggplot(tmp2,aes(FDR,Skew, fill= precision)) +
  geom_tile() + pretty_plot()+scale_fill_viridis_c(limits = c(0.5, 0.85)) +
  geom_text(aes(label = round(precision,2)), color = "white", size = 4)+
  ggtitle("Precision, 14,999 PIP > 0.9, 18,349,  cell-matched e2g + emVar,full sets")+
  coord_equal()

rh <- ggplot(tmp2,aes(FDR,Skew, fill= recall)) +
  geom_tile() + pretty_plot()+scale_fill_viridis_c(option = "magma", limits = c(0, 0.5)) +
  geom_text(aes(label = round(recall,2)), color = "white", size = 4) +
  ggtitle("Recall")+
  coord_equal()

plt_combined <- ggarrange(ph, rh,
                          labels = c("a","b"),
                          ncol = 2, nrow = 1)
plt_combined

cowplot::save_plot(
  'figures/s7ab.R',
  plt_combined,
  base_height = 10,
  base_width =10,
  device = cairo_pdf)



