# Load libraries
library(tidyverse)
library(vroom)
library(GenomicRanges)
library(rtracklayer)
library(plyranges)
library(matrixStats)
library(annotables)
library(seqR)
library(Matrix)
library(mgcv)
library(BuenColors)
library(PNWColors)
library(cowplot)
library(patchwork)
library(ggExtra)
library(ggrastr)
source("code/utils.R")
source("code/theme.R")
options(stringsAsFactors = FALSE)

# Read in processed MPRA results
mpra_df <- vroom("data/preprocess/core_mpra.txt.gz")
mpra_meta_df <- vroom("data/preprocess/mpra_meta.txt.gz")
mpra_df <- mpra_df %>%
  left_join(mpra_meta_df,
            by = c("variant", "cohort")) %>%
  ungroup() 

# Get CRE variants
keep_vars <- mpra_df %>%
  filter(CRE > 0) %>%
  distinct(variant) %>%
  pull(variant)

# Get activity for the set of CRE variants in 3 cell-types
annots_cts_df <- mpra_df %>%
  filter(cell_type %in% c("K562", "HEPG2", "SKNSH")) %>%
  filter(cs_id > 0, CRE > 0, cohort %in% c("BBJ", "UKBB", "GTEx"), type %in% c("PIP10", "CS", "49tissue_PIP50", "3tissue_PIP10", "3tissue_CS")) %>%
  group_by(variant, cell_type, CRE) %>%
  summarize(A_log2FC = max(A_log2FC),
            B_log2FC = max(B_log2FC)) %>%
  ungroup() %>% 
  distinct(variant, cell_type, "allele1" = A_log2FC, "allele2" = B_log2FC) %>%
  pivot_longer(!c("variant", "cell_type")) %>%
  dplyr::rename("log2FC" = value, "allele" = name) %>%
  pivot_wider(names_from = "cell_type", values_from = "log2FC") %>%
  na.omit() 

# Read in oligo sequences
oligos_df <- vroom::vroom("/mnt/sdb/gwas_eqtl_mpra/tables/stable2.txt.gz", delim = "\t", col_names = T) %>%
  tidytable::select(-duplicate) %>%
  tidytable::pivot_longer(!variant) %>%
  tidytable::mutate(name = gsub("_oligo", "", name)) %>%
  tidytable::filter(!grepl("N", value))
oligos_df$length <- nchar(oligos_df$value)

# Add oligo sequences
annots_cts_df <- annots_cts_df %>%
  inner_join(oligos_df %>%
               select(variant, "allele" = name))

# Quantile normalize activity across cell-types
annots_cts_df[,3:dim(annots_cts_df)[2]] <- preprocessCore::normalize.quantiles(annots_cts_df[,3:dim(annots_cts_df)[2]] %>% as.matrix()) %>% tidytable::as_tidytable()

# Shred (count) oligos into k-mers
kmer_length <- 2
oligos_kmer_mat <- count_multimers(oligos_df$value, k_vector = c(1:kmer_length), with_kmer_names = T)
colnames(oligos_kmer_mat) <- gsub("\\.|_|0", "", colnames(oligos_kmer_mat))
kmer_alph <- colnames(oligos_kmer_mat)
kmer_alphr <- chartr("ATGC", "TACG", stringi::stri_reverse(kmer_alph))
colnames(oligos_kmer_mat) <- lapply(1:length(kmer_alph), function(x) {
  c(kmer_alph[x], kmer_alphr[x]) %>%
    unique() %>%
    sort() %>%
    paste0(., collapse = "_")
}) %>%
  unlist()

# Normalize k-mers by oligo length
agg_smatrix <- function(x) {
  nms <- colnames(x)
  uniquenms <- unique(nms)
  sparseMatrix(i = x@i + 1, 
               j = match(nms, uniquenms)[x@j + 1],
               x = x@x,
               dimnames = list(rownames(x), uniquenms),
               repr = "T")
}

oligos_kmer_mat <- agg_smatrix(as(oligos_kmer_mat, "TsparseMatrix")) %>%
  as(., "CsparseMatrix")
kmer_indx <- nchar(gsub("_.*", "", colnames(oligos_kmer_mat)))
oligos_kmer_mat <- lapply(1:kmer_length, function(x) {
  oligos_kmer_mat[, which(kmer_indx == x)] / (oligos_df$length - x + 1)
}) %>%
  do.call("cbind", .)

# Regress k-mer effects out of activity (GC and dinucleotides) with a GAM
in_df <- annots_cts_df %>%
  pivot_longer(!c(variant, allele)) %>%
  left_join(bind_cols(oligos_df %>% 
                        select(variant, name), 
                      oligos_kmer_mat %>% 
                        as.matrix() %>% 
                        as_tibble()),
            by = c("variant", "allele" = "name"))
form <- as.formula(paste0("value ~ ", paste0("s(", names(in_df)[5:dim(in_df)[2]], ")", collapse = " + ")))
set.seed(1234)
fit <- bam(form, data = in_df, method = "fREML")
summary(fit)
resid_fit <- resid(fit)

# Get adjusted and unadjusted activity values 
in_adj_df <- in_df %>%
  mutate(value = resid_fit) %>%
  select(variant, allele, name, value) %>%
  mutate(name = paste0(name, "_adj")) %>% 
  pivot_wider(names_from = "name", values_from = value)
annots_cts_df <- annots_cts_df %>%
  left_join(in_adj_df,
            by = c("variant", "allele"))

# Read in counts of TF motifs in each oligo
motifcount_df <- vroom::vroom("data/motifs/motifcount.txt.gz", altrep = T)

# Create motif matrix for regression
motifcount_wide_mat <- annots_cts_df %>%
  mutate(rowidx = 1:dim(annots_cts_df)[1]) %>%
  distinct(variant, allele, rowidx) %>%
  left_join(motifcount_df %>%
              filter(variant %in% unique(annots_cts_df$variant)) %>%
              select(-motif_db)) %>%
  tidytext::cast_sparse(rowidx, motif, nadj)

# Write out files for regression
annots_cts_df %>%
  vroom::vroom_write("data/motifs/mpra_CRE_120k_3cts.txt.gz", delim = "\t", col_names = T)
motifcount_wide_mat %>%
  as.matrix() %>%
  as_tibble() %>%
  vroom::vroom_write("data/motifs/motifcount_CRE_120k_3cts.txt.gz", delim = "\t", col_names = T)
#annots_cts_df <- vroom::vroom("data/motifs/mpra_CRE_120k_3cts.txt.gz")
#motifcount_wide_mat <- vroom::vroom("data/motifs/motifcount_CRE_120k_3cts.txt.gz")

# Create data frames for regression
y_dt <- annots_cts_df[,3:8] %>%
  collapse::qM()
x_dt <- motifcount_wide_mat %>%
  as.matrix()
gc()

# Relatively speedy parallelized OLS
flm_func2 <- function(x) {
  out <- lapply(1:dim(y_dt)[2], function(i) {
    tmp <- collapse::flm(y_dt[,i], x, method = "eigen", eigen.method = 2, return.raw = T, add.icpt = TRUE)[c(1, 2, 4)]
    pval <- -2 * pt(abs(tmp$coefficients[2] / tmp$se[2]), tmp$df.residual, lower.tail = F, log.p = T)
    tmp_prop_elem <- sum(x > 0) / length(x)
    list("beta" = as.numeric(tmp$coefficients[2]), "se" = as.numeric(tmp$se[2]), "log10_pval" = as.numeric(pval), "prop_elem" = tmp_prop_elem)
  }) %>% 
    collapse::rowbind(return = "list") %>%
    collapse::t_list() %>%
    set_names(colnames(y_dt))
  # Deal with memory used up by finished / phantom nodes every once in a while
  invisible(if(runif(1) > 0.95) {gc(v = FALSE)})
  return(out)
}

# Run OLS
n_cores <- 32
system.time(motifbeta_df <- collapse::dapply(x_dt, flm_func2, MARGIN = 2, mc.cores = n_cores, parallel = T, return = "data.frame") %>% 
              as.list() %>% 
              collapse::unlist2d(recursive = T, DT = T, idcols = c("motif", "cell_type")))

# Write out results of activity ~ motif count regression
motifbeta_df %>%
  vroom::vroom_write("data/motifs/motifbeta_CRE_120k_3cts.txt.gz")
#motifbeta_df <- vroom::vroom("data/motifs/motifbeta_CRE_120k_3cts.txt.gz")

# Load motif metadata
load("data/motifs/motifs_processed_20230225.RData")
motifs_pfml <- c(hocomoco_pfml, jaspar_pfml, cisbp_pfml, vierstra_pfml)

# Motif regression results long to wide
motifbeta_wide_df <- motifbeta_df %>%
  pivot_wider(names_from = c("cell_type"), values_from = c("beta", "se", "log10_pval")) %>%
  dplyr::filter(prop_elem > 0.001) %>%
  left_join(motifs_meta, 
            by = "motif") %>%
  dplyr::filter(!is.na(TF))

# Correlation of motif effects across cell-types
motifbeta_wide_df %>%
  dplyr::select(beta_HEPG2_adj, beta_K562_adj, beta_SKNSH_adj) %>%
  cor()

# Plot effects of motifs on activity (Fig. 1e)
p1 <- motifbeta_wide_df %>%
  dplyr::select(motif, TF, prop_elem, beta_SKNSH, beta_HEPG2_adj, beta_K562_adj) %>%
  ggplot(aes(x = beta_HEPG2_adj, y = beta_K562_adj)) + 
  geom_point(color = 'grey', size = 1, alpha = 0.7) +
  geom_point(data = motifbeta_wide_df %>% filter(grepl('ETS|ELS|ELK|KLF|SP|ELK', TF)), color = pnw_palette("Bay", 8, type = "continuous")[5], size = 1) +
  geom_point(data = motifbeta_wide_df %>% filter(grepl('GATA|TAL|STAT5|TRPS1', TF)), color = pnw_palette("Bay", 8, type = "continuous")[8], size = 1) +
  geom_point(data = motifbeta_wide_df %>% filter(grepl('GFI1', TF)), color = pnw_palette("Bay", 8, type = "continuous")[1], size = 1) +
  geom_point(data = motifbeta_wide_df %>% filter(grepl('HNF4|HNF1|CEBP|TP5|TP6|TP7|PPAR', TF)), color = pnw_palette("Bay", 8, type = "continuous")[6], size = 1) +
  geom_point(data = motifbeta_wide_df %>% filter(grepl('P53', TF)), color = pnw_palette("Bay", 8, type = "continuous")[3], size = 1) +
  geom_point(data = motifbeta_wide_df %>% filter(grepl('REST|SNAI', TF)), color = pnw_palette("Bay", 8, type = "continuous")[4], size = 1) +
  pretty_plot() +
  geom_abline() +
  ylab("Activity contribution (K562)") +
  xlab("Activity contribution (HepG2)")
p1

# Save plot (Fig. 1e)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(2))
cowplot::save_plot(
  paste0("figures/fig1/fig1e_motifbeta.pdf"),
  plt_combined,
  base_height = 2.2,
  base_width = 2.2,
  device = cairo_pdf
)





