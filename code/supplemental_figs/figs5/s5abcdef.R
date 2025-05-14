# Load libraries
library(tidyverse)
library(BuenColors)
library(rtracklayer)
library(reshape2)
library(ComplexHeatmap)
library(GenomicRanges)
library(ggpubr)
library(binom)
library(dplyr)
library(vroom)
source("code/utils.R")
options(stringsAsFactors = FALSE)

# Read in processed MPRA results and filter to GTEx
mpra <- vroom("/mnt/sdb/gwas_eqtl_mpra/data/preprocess/core_mpra.txt.gz") %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0)) %>%
  dplyr::filter(cohort=='GTEx')

# Let's look at cell type specific correlations
DHS_Meuleman.meta <- read_delim("/mnt/sdb/gwas_eqtl_mpra/data/annotations/meuleman/DHS_Index_and_Vocabulary_metadata.tsv", delim = "\t", col_names = T)
DHS_Meuleman.meta <- DHS_Meuleman.meta[-734,]

# Read in DHS quantitative
in1.gr <- readRDS("/mnt/sdb/gwas_eqtl_mpra/data/annotations/meuleman/meuleman_en_gr.rds")
in2.gr <- in1.gr[,0]

# Convert variants to granges
variant.gr <- mpra %>%
  dplyr::distinct(chromosome, position, variant) %>%
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


indf <- mpra %>%
  distinct(variant,cell_type) %>%
  left_join(dhs_vars,
            by = "variant") %>%
  dplyr::select(variant, cell_type, any_of(altius)) %>%
  distinct() 

# Get matching tissues for DHS
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
variant_gr <- mpra %>%
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

# Find overlaps with E2G
idx <- findOverlaps(variant_gr, e2g_gr)
e2g_filt_df <- bind_cols(variant_gr[idx@from] %>% 
                           as_tibble() %>%
                           dplyr::select(variant),
                         e2g_gr[idx@to] %>%
                           as_tibble()) %>%
  distinct()

# Keep E2G for MPRA tissues
e2g_filt_cts_df <- e2g_filt_df %>%
  dplyr::filter(cell_type %in% c("K562", "HepG2", "HCT116", "SK-N-DZ", "SK-N-MC")) %>%
  dplyr::mutate(cell_type = ifelse(cell_type == "HepG2", "HEPG2", cell_type),
                cell_type = ifelse(cell_type %in% c("SK-N-DZ", "SK-N-MC"), "SKNSH", cell_type))



# annotations dfs:
out
e2g_filt_cts_df

#################  choose 1 of the following ################

# s5a,b
# if just CREs (any):
# Get positive set variants
positives_df <- mpra %>%
  dplyr::filter(consequence %ni% c("synonymous","missense","LoF"),
                CRE > 0) %>% 
  dplyr::mutate(causal = T) %>%
  distinct(variant, cell_type, log2FC, logPadj_BF, CRE, causal, library)
# get negative variants
negatives_df <- mpra %>%
  dplyr::filter(consequence %ni% c("synonymous","missense","LoF"),
                CRE == 0) %>% 
  dplyr::mutate(causal = F) %>%
  distinct(variant, cell_type, log2FC, logPadj_BF, CRE, causal, library)

# s5c,d
# cell-matching e2g
# Get positive set variants
positives_df <- mpra %>%
  distinct(variant, cell_type, log2FC, logPadj_BF, CRE, consequence) %>%
  left_join(e2g_filt_cts_df %>%
              distinct(variant, cell_type) %>%
              dplyr::mutate(is.e2g = 1),
            by = c('variant','cell_type')) %>%
  replace(is.na(.),0) %>%
  dplyr::filter(consequence %ni% c("synonymous","missense","LoF"),
                is.e2g > 0) %>% 
  dplyr::mutate(causal = T) 
negatives_df <- mpra %>%
  distinct(variant, cell_type, log2FC, logPadj_BF, CRE, consequence) %>%
  left_join(e2g_filt_cts_df %>%
              distinct(variant, cell_type) %>%
              dplyr::mutate(is.e2g = 1),
            by = c('variant','cell_type')) %>%
  replace(is.na(.),0) %>%
  dplyr::filter(consequence %ni% c("synonymous","missense","LoF")) %>%
  dplyr::group_by(variant) %>%
  dplyr::filter(all(is.e2g == 0)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(causal = F) 


#s5e,f
# cell-matching DHS
# Get positive set variants
positives_df <- mpra %>%
  distinct(variant, cell_type, log2FC, logPadj_BF, consequence) %>%
  left_join(out %>%
              dplyr::filter(cell_type == tissue, value > 0) %>%
              distinct(variant, cell_type, value),
            by = c('variant','cell_type')) %>%
  replace(is.na(.),0) %>%
  dplyr::rename(is.CRE = value) %>%
  dplyr::mutate(is.CRE = ifelse(is.CRE > 0, 1,0)) %>%
  dplyr::filter(consequence %ni% c("synonymous","missense","LoF"),
                is.CRE > 0) %>% 
  dplyr::mutate(causal = T) 

negatives_df <- mpra %>%
  distinct(variant, cell_type, log2FC, logPadj_BF, consequence) %>%
  left_join(out %>%
              dplyr::filter(cell_type == tissue, value > 0) %>%
              distinct(variant, cell_type, value),
            by = c('variant','cell_type')) %>%
  replace(is.na(.),0) %>%
  dplyr::rename(is.CRE = value) %>%
  dplyr::mutate(is.CRE = ifelse(is.CRE > 0, 1,0)) %>%
  dplyr::filter(consequence %ni% c("synonymous","missense","LoF")) %>%
  dplyr::group_by(variant) %>%
  dplyr::filter(all(is.CRE == 0)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(causal = F) 

#################################

# Bind together positives and negatives
mpra_mc <- bind_rows(positives_df, negatives_df)

# Prepare precision/recall input dataframe
prset_df <- mpra_mc %>%
  dplyr::group_by(variant) %>%
  arrange(-logPadj_BF) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()


print(prset_df %>% distinct(variant, causal) %>% count(causal))


# precision and recall function for gridsearch
prec_rec_active <- function(pthresh,a, df) {
  print(paste(pthresh,a, sep = ', '))
  
  out_df <- df %>%
    dplyr::mutate(active_test = ifelse(logPadj_BF >= -log10(pthresh) & abs(log2FC) >= a, T, F)) %>%
    dplyr::select(active_test, causal) %>%
    dplyr::mutate(annot = case_when(active_test == T ~ 1,
                                    active_test == F ~ 0)) %>%
    dplyr::filter(!is.na(active_test)) %>%
    dplyr::count(causal, annot)
  
  # Compute precision, recall, and 95% CIs
  out_df <- out_df %>%
    dplyr::summarize(TP = sum(n[causal == T & annot == T]),
                     FP = sum(n[causal == F & annot == T]),
                     FN = sum(n[causal == T & annot == F]),
                     TN = sum(n[causal == F & annot == F]),
                     total_pos = sum(n[annot == T]),
                     total = sum(n)) 
  out_df <- out_df %>%
    dplyr::mutate(precision = TP / (TP + FP),
                  recall = TP / (TP + FN),
                  prec_upper = binom.confint(TP, (TP + FP), method = "wilson")$upper,
                  prec_lower = binom.confint(TP, (TP + FP), method = "wilson")$lower,
                  rec_upper = binom.confint(TP, (TP + FN), method = "wilson")$upper,
                  rec_lower = binom.confint(TP, (TP + FN), method = "wilson")$lower) %>%
    dplyr::mutate(FDR = pthresh,
                  activity = a)
  return(out_df)
}

# gridsearch across p-values and activity thresholds
it_p <- rep(c(0.00001,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1),10)
it_a <- c(rep(c(0),10),rep(c(0.1),10),rep(c(0.2),10),rep(c(0.3),10),rep(c(0.4),10),rep(c(0.5),10),rep(c(0.75),10),rep(c(1),10),rep(c(2),10),rep(c(3),10))
data <- mapply(prec_rec_active,it_p,it_a,  MoreArgs=list(df = prset_df))
tmp <- t(data) %>% as.data.frame()
tmp[is.na(tmp)] <- 0

# prepare data for plotting
tmp <- transform(tmp, FDR = as.character(FDR))%>%
  transform(TN = as.numeric(TN)) %>%
  transform(FN = as.numeric(FN)) %>%
  transform(FP = as.numeric(FP)) %>%
  transform(activity = as.numeric(activity)) %>%
  transform(precision = as.numeric(precision)) %>%
  transform(recall = as.numeric(recall))  %>%
  transform(TP = as.numeric(TP)) %>%
  transform(total = as.numeric(total)) %>% 
  transform(total_pos = as.numeric(total_pos)) %>% 
  transform(prec_upper = as.numeric(prec_upper)) %>% 
  transform(rec_upper = as.numeric(rec_upper)) %>% 
  transform(prec_lower = as.numeric(prec_lower)) %>% 
  transform(rec_lower = as.numeric(rec_lower)) 
tmp$FDR <- factor(tmp$FDR, levels = c("1e-05","1e-04","5e-04","0.001","0.005","0.01","0.05","0.1","0.5","1"))
tmp$activity <- factor(tmp$activity, levels = c("0","0.1", "0.2", "0.3",  "0.4","0.5","0.75","1","2","3"))


# Plot
# choose appropriate title 
ph <- ggplot(tmp,aes(FDR,activity, fill= precision)) +
  geom_tile(color='white', size = .1) + pretty_plot()+scale_fill_viridis_c(limits = c(0.2, 0.7),  oob = scales::squish) +
  geom_text(aes(label = round(precision,2)), color = "white", size = 4)+
  xlab('p-value')+
  #ggtitle('Precision, 41,572 positive in any CRE, 103,384 negative in no CRE')+
  #ggtitle('Precision, 14,995 positive in cell-matched e2g, 130,001 negative in no cell-matched e2g')+
  ggtitle('Precision, 8,877 positive in cell-matched DHS, 136,079 negative in no cell-matched DHS')+
  coord_equal()
rh <- ggplot(tmp,aes(FDR,activity, fill= recall)) +
  geom_tile(color = 'white',size = 0.1) + pretty_plot()+scale_fill_viridis_c(option = "magma", limits = c(0.1, 0.9), oob = scales::squish) +
  geom_text(aes(label = round(recall,2)), color = "white", size = 4)+
  xlab('p-value')+
  ggtitle('Recall')+
  coord_equal()

plt_combined <- ggarrange(ph, rh,
                          labels = c("a","b"),
                          ncol = 2, nrow = 1)
plt_combined

#outt <- "figures/s5b.pdf"
#outt <-"figures/s5cd.pdf"
outt <- "figures/s5ef.pdf"

print(outt)
cowplot::save_plot(
  outt,
  plt_combined,
  base_height = 10,
  base_width = 10,
  device = cairo_pdf)

