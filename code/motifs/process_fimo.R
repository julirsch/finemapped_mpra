# Load libraries
library(tidyverse)
library(MotifDb)
library(universalmotif)
library(BiocParallel)
library(GenomicRanges)

# Parameters
args = commandArgs(trailingOnly = TRUE)
mdb <- args[1]
chrom <- args[2]
#mdb <- "vierstra"
#chrom <- 21
inmotifs <-  paste0("mpra_chr", chrom, "_", mdb, ".out.gz")
ncores <- 200 # allow 3GB per core for 10000 sequence blocks
outdir <- "/mnt/sdb/gtex_mpra/data/fimo/"

# Read in MPRA
mpra <- vroom::vroom("/mnt/sdb/gwas_eqtl_mpra/data/preprocess/core_mpra.txt.gz")

# Load motifs
load("/mnt/sdb/gtex_mpra/code/motifs_processed_20221025.RData")

# Read in oligos
oligos <- vroom::vroom("/mnt/sdb/gwas_eqtl_mpra/tables/stable2.txt.gz", delim = "\t", col_names = T) %>%
  dplyr::select(-duplicate) %>%
  pivot_longer(!variant) %>%
  dplyr::mutate(name = gsub("_oligo", "", name))

# Read in FIMO motif scores and positions
motif_df <- vroom::vroom(paste0("/mnt/sdb/gtex_mpra/data/fimo/", inmotifs), col_names = c("motif", "vtmp", "start", "stop", "strand", "score", "pval", "qval", "seq"), skip = 1) %>%
  dplyr::select(motif, vtmp, start, stop, strand, score, pval, seq) %>%
  dplyr::mutate(vtmp = gsub("_oligo", "", vtmp))
gc()

# Filter to only SNPs
snps <- mpra %>%
  distinct(variant, allele1, allele2) %>%
  dplyr::filter(nchar(allele1) == 1, nchar(allele2) == 1) %>%
  .$variant
rm(mpra)

# Get instances of difference by SNP
motif_a1_df <- motif_df %>%
  dplyr::mutate(variant = gsub("_.*", "", vtmp),
                allele = gsub(".*_", "", vtmp)) %>%
  dplyr::select(motif, vtmp, variant, allele, start, stop, strand, score, pval, seq) %>%
  dplyr::filter(allele == "allele1",
                score > 0,
                variant %in% snps)
motif_a2_df <- motif_df %>%
  dplyr::mutate(variant = gsub("_.*", "", vtmp),
                allele = gsub(".*_", "", vtmp)) %>%
  dplyr::select(motif, vtmp, variant, allele, start, stop, strand, score, pval, seq) %>%
  dplyr::filter(allele == "allele2",
                score > 0,
                variant %in% snps)
motif_variants_df <- motif_a1_df %>%
  full_join(motif_a2_df,
            by = c("motif", "variant", "start", "stop", "strand")) %>%
  dplyr::mutate(pval = pmin(pval.x, pval.y, na.rm = T)) %>%
  dplyr::select(motif, variant, start, stop, strand, pval, "seq_allele1" = seq.x, "score_allele1" = score.x, "seq_allele2" = seq.y, "score_allele2" = score.y) %>%
  distinct() %>%
  dplyr::filter(abs(score_allele1 - score_allele2) > 0.01 | is.na(score_allele1) | is.na(score_allele2))
rm(motif_a1_df); rm(motif_a2_df); gc()
motif_variants_seqs_df <- motif_variants_df %>%
  distinct(motif, "seq" = seq_allele1) %>%
  bind_rows(motif_variants_df %>%
              distinct(motif, "seq" = seq_allele2)) %>%
  na.omit() %>%
  distinct()

# Get missing sequences
motif_variants_missing_df <- bind_rows(motif_variants_df %>%
                                         dplyr::filter(is.na(seq_allele1)) %>%
                                         distinct(motif, variant, start, stop, strand, seq_allele2, seq_allele1) %>%
                                         left_join(oligos %>%
                                                     dplyr::filter(name == "allele1") %>%
                                                     dplyr::select(variant, value),
                                                   by = "variant") %>%
                                         dplyr::mutate(tmp = substr(value, start, stop)) %>%
                                         dplyr::select(-value) %>%
                                         dplyr::mutate(seq_allele1 = case_when(strand == "+" ~ tmp,
                                                                               strand == "-" ~ stringi::stri_reverse(chartr("ATGC", "TACG", tmp)))) %>%
                                         dplyr::select(-tmp),
                                       motif_variants_df %>%
                                         dplyr::filter(is.na(seq_allele2)) %>%
                                         distinct(motif, variant, start, stop, strand, seq_allele1, seq_allele2) %>%
                                         left_join(oligos %>%
                                                     dplyr::filter(name == "allele2") %>%
                                                     dplyr::select(variant, value),
                                                   by = "variant") %>%
                                         dplyr::mutate(tmp = substr(value, start, stop)) %>%
                                         dplyr::select(-value) %>%
                                         dplyr::mutate(seq_allele2 = case_when(strand == "+" ~ tmp,
                                                                               strand == "-" ~ stringi::stri_reverse(chartr("ATGC", "TACG", tmp)))) %>%
                                         dplyr::select(-tmp))

# Get unique sequences to rescore
# Optional: apply more stringent filters here to reduce burden on matches
motif_df <- motif_df %>%
  dplyr::filter(pval <= 0.0001, score > 0)
unique_ms_df <- motif_df %>%
  distinct(motif, seq)

# Combine sequences for matches and for disruption to rescore
unique_ms_df <- bind_rows(unique_ms_df,
                          motif_variants_seqs_df,
                          motif_variants_missing_df %>%
                            dplyr::select(motif, seq_allele1, seq_allele2) %>%
                            pivot_longer(!motif) %>%
                            dplyr::select(motif, "seq" = value) %>%
                            distinct()
                          )

# Set background for FIMO
bkgd <- c(0.281774, 0.222020, 0.228876, 0.267330)
adjbkgd <- c(0.274552, 0.225448, 0.225448, 0.274552)
names(bkgd) <- c("A", "C", "G", "T")
names(adjbkgd) <- c("A", "C", "G", "T")

# Adjust PWMs to ~FIMO scale for easy score computation
motif_mdb <- eval(parse(text = paste0(mdb, "_mdb")))
pwmlist_scale <- lapply(1:length(motif_mdb), function(y) {
  motif <- motif_mdb[[y]]
  adjmotif <- (motif * 20 + 0.1 * bkgd) / (20 + 0.1)
  value <- log2(adjmotif / adjbkgd)
  value
})
names(pwmlist_scale) <- names(motif_mdb)

# Get min and max for scaling scores to %iles
pwmlist_minmax <- lapply(1:length(pwmlist_scale), function(y) {
  motif <- pwmlist_scale[[y]]
  mmin <- sum(apply(motif, 2, min))
  mrand <- sum(colSums(motif) / 4)
  mmax <- sum(apply(motif, 2, max))
  data.frame("motif" = names(pwmlist_scale)[y], "random" = mrand, "min" = mmin, "max" = mmax)
}) %>%
  bind_rows() %>%
  as_tibble()

# Set up parallelization
indx <- seq(1:dim(unique_ms_df)[1])
blocks <- floor(dim(unique_ms_df)[1] / 10000)
indxl <- split(indx, ceiling(seq_along(indx)/(length(indx) / blocks)))
params <- BiocParallel::bpparam()
params$workers <- ncores
  
# Function to score motifs on each split
scoreRun <- function(z) {
  out <- lapply(z, function(y) {
    inmotif <- unique_ms_df[y,]$motif
    inseq <- strsplit(unique_ms_df[y,]$seq, "")[[1]]
    adjmotif <- pwmlist_scale[[inmotif]]
    tmpout <- lapply(1:length(inseq), function(x) {
      nuc <- inseq[x]
      adjmotif[nuc, x]
    }) %>%
      unlist() %>%
      as_tibble() %>%
      dplyr::summarize(score = sum(value))
  }) %>% 
    bind_rows() 
  unique_ms_df[z,] %>%
    bind_cols(out)
}

# Run function on each block and combine
results <- bplapply(indxl, scoreRun, BPPARAM = params) %>%
  bind_rows()
BiocParallel::bpstop()
gc()

# Fix score and compute percentiles
motif_df <- motif_df %>%
  dplyr::select(-score) %>%
  left_join(results %>%
              left_join(pwmlist_minmax,
                        by = "motif"),
            by = c("motif", "seq")) %>%
  dplyr::mutate(pct = (score - min) / (max - min),
                pcta = (score - random) / (max - random)) %>%
  dplyr::select(motif, vtmp, start, stop, strand, score, pct, pcta, pval, seq) %>%
  distinct()

# Add score for variants
motif_variants_df <- motif_variants_df %>%
  left_join(motif_variants_missing_df,
            by = c("motif", "variant", "start", "stop", "strand")) %>%
  dplyr::mutate(seq_allele1 = case_when(is.na(seq_allele1.x) ~ seq_allele1.y,
                                         T ~ seq_allele1.x),
                seq_allele2 = case_when(is.na(seq_allele2.x) ~ seq_allele2.y,
                                         T ~ seq_allele2.x)) %>%
  dplyr::select(motif, variant, start, stop, strand, pval, seq_allele1, seq_allele2) %>%
  left_join(results,
            by = c("motif", "seq_allele1" = "seq")) %>%
  dplyr::rename("score_allele1" = score) %>%
  left_join(results,
            by = c("motif", "seq_allele2" = "seq")) %>%
  dplyr::rename("score_allele2" = score) %>%
  left_join(pwmlist_minmax,
            by = "motif") %>%
  dplyr::mutate(pct_allele1 = (score_allele1 - min) / (max - min),
                pct_allele2 = (score_allele2 - min) / (max - min),
                pcta_allele1 = (score_allele1 - random) / (max - random),
                pcta_allele2 = (score_allele2 - random) / (max - random)) %>%
  dplyr::select(motif, variant, start, stop, strand, score_allele1, score_allele2, pct_allele1, pct_allele2, pcta_allele1, pcta_allele2, pval) %>%
  distinct() %>%
  dplyr::filter(abs(score_allele1 - score_allele2) > 0.01)

# Write out motif matches
motif_df <- motif_df %>%
  dplyr::mutate(variant = gsub("_.*", "", vtmp),
                allele = gsub(".*_", "", vtmp)) %>%
  dplyr::select(motif, variant, allele, start, stop, strand, score, pct, pcta, pval) %>%
  dplyr::filter(pval <= 0.0001, score > 0) %>%
  distinct()
vroom::vroom_write(motif_df, paste0(outdir, "mpra_chr", chrom, "_", mdb, "_match.txt.gz"), delim = "\t", col_names = T)

# Write out motif disruptions
motif_variants_df <- motif_variants_df %>%
  dplyr::filter(score_allele1 > 0 | score_allele2 > 0, pct_allele1 > 0 | pct_allele2 > 0) %>% 
  dplyr::filter(abs(score_allele1 - score_allele2) > 0.01)
vroom::vroom_write(motif_variants_df, paste0(outdir, "mpra_chr", chrom, "_", mdb, "_disrupt.txt.gz"), delim = "\t", col_names = T)

# Wriute out number of (adjusted) motifs in each element
motif_tmp <- motif_df %>%
  dplyr::filter(pval <= 0.0001, score > 0, pct >= 0.85) %>%
  distinct(variant, allele, motif, start, stop, strand) %>%
  dplyr::mutate(size = stop - start)
motif_span <- motif_tmp %>%
  dplyr::mutate("vam" = paste0(variant, "_", allele, ";", motif)) %>%
  dplyr::select(-strand) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "vam", 
                           start.field =  "start", 
                           end.field = "stop",
                           keep.extra.columns = F) %>%
  GenomicRanges::reduce() %>%
  as_tibble() %>%
  dplyr::select(seqnames, width) %>%
  group_by(seqnames) %>%
  dplyr::summarize(sum = sum(width)) %>%
  dplyr::mutate(va = gsub(";.*", "", seqnames),
                motif = gsub(".*;", "", seqnames)) %>%
  dplyr::mutate(variant = gsub("_.*", "", va),
                allele = gsub(".*_", "", va)) %>%
  dplyr::select(-seqnames, -va)
motif_count <- motif_span %>%
  left_join(motif_tmp %>%
              distinct(motif, size),
            by = "motif") %>%
  dplyr::mutate(nadj = sum / size) %>%
  distinct(variant, allele, motif, nadj)
vroom::vroom_write(motif_count, paste0(outdir, "mpra_chr", chrom, "_", mdb, "_count.txt.gz"), delim = "\t", col_names = T)

# Clean up
rm(list = ls())
gc()

