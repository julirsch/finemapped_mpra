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
source("/code/utils.R")

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

#altius columns corresponding to tissues of interest
altius <- DHS_Meuleman.meta %>%
  dplyr::rename("tissue" = `Biosample name`,
                "celltype" = `Altius Aggregation ID`,
                "system" = System) %>%
  dplyr::select(tissue, celltype) %>%
  dplyr::mutate(tissue = ifelse(tissue == "HepG2", "HEPG2", tissue)) %>% 
  dplyr::filter(tissue %in% unique(mpra$cell_type))  %>%
  .$celltype


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


# join with cell-type specific CREs and filter down to only cell type specific CRE
mpra.overall <- mpra_mc %>% 
  left_join(out %>%
              dplyr::filter(cell_type == tissue, value > 0) %>%
              distinct(variant, cell_type, value),
            by = c('variant','cell_type')) %>%
  replace(is.na(.),0) %>%
  dplyr::rename(is.CRE = value) %>%
  dplyr::mutate(is.CRE = ifelse(is.CRE > 0, 1,0))


prset_df <- mpra.overall %>%
  dplyr::mutate(is.CRE.active = ifelse(active + is.CRE == 2, T, F)) %>%
  group_by(variant) %>%
  arrange(-is.CRE.active,-Skew_logPadj) %>%
  filter(row_number() == 1) %>%
  ungroup()


print(prset_df %>% distinct(variant, causal) %>% count( causal))

# gridsearch
it_p <- rep(c(0.001,0.01,0.05,0.1,0.15,0.2,0.25,0.5,0.75,1),10)
it_a <- c(rep(c(0),10),rep(c(0.25),10),rep(c(0.5),10),rep(c(0.75),10),rep(c(1),10),rep(c(1.25),10),rep(c(1.5),10),rep(c(1.75),10),rep(c(2),10),rep(c(3),10))

prec_rec_emvar <- function(pthresh,athresh, df) {
  print(paste(pthresh, athresh, sep = ', '))
  # Get 2x2 table
  out_df <- df %>%
    # use below for TARV analysis
    dplyr::mutate(emVar_test = if_else(is.CRE & active & Skew_logPadj >= -log10(pthresh) & !is.na(Skew_logPadj) &  abs(log2Skew) >= athresh, T, F)) %>%
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

# Run gridsearch
data2 <- mapply(prec_rec_emvar,it_p,it_a, MoreArgs=list(df = prset_df))
tmp2 <- t(data2) %>% as.data.frame()
tmp2[is.na(tmp2)] <- 0

tmp2

# Prepare for plotting
tmp2 <- transform(tmp2, FDR = as.character(FDR))%>% 
  transform(Skew = as.numeric(Skew)) %>%
  transform(precision = as.numeric(precision)) %>%
  transform(recall = as.numeric(recall))  %>%
  transform(TP = as.numeric(TP)) %>%
  #transform(total = as.numeric(total)) %>% 
  transform(prec_upper = as.numeric(prec_upper)) %>%
  transform(prec_lower = as.numeric(prec_lower))

tmp2$FDR <- factor(tmp2$FDR, levels = c("0.001","0.01","0.05","0.1","0.15","0.2","0.25","0.5","0.75","1"))
tmp2$Skew <- factor(tmp2$Skew, levels = c("0","0.25","0.5","0.75","1","1.25","1.5","1.75","2","3"))

# Plot
ph <- ggplot(tmp2,aes(FDR,Skew, fill= precision)) +
  geom_tile() + pretty_plot()+scale_fill_viridis_c(limits = c(0.5, 0.85)) +
  geom_text(aes(label = round(precision,2)), color = "white", size = 4)+
  ggtitle("Precision, 14,999 PIP > 0.9, 18,349,  cell-matched DHS + emVar,full sets")+
  coord_equal()

rh <- ggplot(tmp2,aes(FDR,Skew, fill= recall)) +
  geom_tile() + pretty_plot()+scale_fill_viridis_c(option = "magma", limits = c(0, 0.5)) +
  geom_text(aes(label = round(recall,2)), color = "white", size = 4) +
  ggtitle("Recall")+
  coord_equal()

plt_combined <- ggarrange(ph, rh,
                          labels = c("c","d"),
                          ncol = 2, nrow = 1)
plt_combined

cowplot::save_plot(
  paste0('figures/s7cd.pdf'),
  plt_combined,
  base_height = 10,
  base_width =10,
  device = cairo_pdf)

