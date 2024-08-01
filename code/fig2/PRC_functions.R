# Load libraries
library(tidyverse)
library(binom)
source("code/utils.R")

# Function to create positive and negative variant sets
make_prset <- function(df, type) {
  
  # For eQTLs
  if(type == "eQTLs") {
    
    # Get positive set variants
    positives_df <- df %>%
      dplyr::filter(pip > 0.9, 
                    abs(z) > qnorm(1 - 5 * 10^-8), 
                    type %in% c("49tissue_PIP50", "3tissue_PIP10", "3tissue_CS"), 
                    consequence %ni% c("synonymous","missense","LoF")) %>% 
      dplyr::group_by(variant) %>%
      dplyr::mutate(pip = max(pip, na.rm = T)) %>%
      dplyr::filter(pip == max(pip),
                    !is.na(pip)) %>%
      filter(row_number() == 1) %>%
      ungroup() %>%
      dplyr::mutate(causal = T)
    
    # Get negative set variants
    # Not enough variants to match if PIP < 0.01
    negatives_df <- df %>%
      group_by(variant) %>%
      dplyr::mutate(pip = max(pip, na.rm = T)) %>%
      dplyr::filter(pip == max(pip),
                    !is.na(pip)) %>%
      filter(row_number() == 1) %>%
      ungroup() %>%
      dplyr::filter(pip < 0.02, 
                    type == '3tissue_CS', 
                    consequence %ni% c("synonymous","missense","LoF"), 
                    variant %ni% positives_df$variant) %>% 
      ungroup() %>%
      dplyr::mutate(causal = F)
    
    # Match positives and negatives 1:1
    set.seed(123)
    sampsize <- dim(positives_df)[1]
    prset_df <- positives_df %>%
      bind_rows(negatives_df %>%
                  sample_n(sampsize))
    return(prset_df)
  } 
  
  # For GWAS
  else if(type == "complex_traits") {
    
    # BBJ dataset is too small for PRC
    df <- df %>% 
      dplyr::filter(cohort != "BBJ")
    
    # Get positive set variants
    positives_df <- df %>%
      dplyr::filter(pip > 0.9, 
                    pchisq(chisq_marginal, 1, log.p = TRUE, lower.tail = F) / -log(10) > -log10(5 * 10^-8), 
                    type %in% c("CS", "PIP10"),
                    consequence %ni% c("synonymous", "missense", "LoF")) %>% 
      dplyr::group_by(variant) %>%
      dplyr::mutate(pip = max(pip, na.rm = T)) %>%
      dplyr::filter(pip == max(pip),
                    !is.na(pip)) %>%
      filter(row_number() == 1) %>%
      ungroup() %>%
      dplyr::mutate(causal = T)
    
    # Get negative set variants
    negatives_df <- df %>%
      group_by(variant) %>%
      dplyr::mutate(pip = max(pip, na.rm = T)) %>%
      dplyr::filter(pip == max(pip),
                    !is.na(pip)) %>%
      filter(row_number() == 1) %>%
      ungroup() %>%
      dplyr::filter(pip < 0.01, 
                    type == "CS", 
                    consequence %ni% c("synonymous", "missense", "LoF"), 
                    variant %ni% positives_df$variant) %>% 
      ungroup() %>%
      dplyr::mutate(causal = F)
    
    # Match positives and negatives 1:1
    set.seed(123)
    sampsize <- dim(positives_df)[1]
    prset_df <- positives_df %>%
      bind_rows(negatives_df %>%
                  sample_n(sampsize))
    return(prset_df)
  }
}

# Function to add joint annotation to df
add_annots <- function(annot_name, in_exprs, df) {
  
  # Add joint annotation
  df %>%
    dplyr::mutate(tmp = case_when(!! in_exprs[[1]] & !! in_exprs[[3]] ~ 1,
                                         !! in_exprs[[1]] & !! in_exprs[[4]] ~ 0,
                                         !! in_exprs[[2]] & !! in_exprs[[3]] ~ 0,
                                         !! in_exprs[[2]] & !! in_exprs[[4]] ~ 0,
                                         T ~ NA)) %>%
    dplyr::rename(!!annot_name := tmp)

}

# Function to compute precision and recall
prec_rec <- function(in_annot, df) {

  # Get 2x2 table
  out_df <- df %>%
    dplyr::select(!! sym(in_annot), causal) %>%
    dplyr::mutate(annot = case_when(!! sym(in_annot) > 0 ~ 1,
                                    !! sym(in_annot) == 0 ~ 0)) %>%
    dplyr::filter(!is.na(!! sym(in_annot))) %>%
    dplyr::count(causal, annot)
  
  # Compute precision, recall, and 95% CIs
  out_df %>%
    dplyr::summarize(TP = sum(n[causal == T & annot == T]),
                     FP = sum(n[causal == F & annot == T]),
                     FN = sum(n[causal == T & annot == F]),
                     total_pos = sum(n[annot == T]),
                     total = sum(n)) %>%
    dplyr::mutate(precision = TP / (TP + FP),
                  recall = TP / (TP + FN),
                  prec_upper = binom.confint(TP, (TP + FP), method = "wilson")$upper,
                  prec_lower = binom.confint(TP, (TP + FP), method = "wilson")$lower,
                  rec_upper = binom.confint(TP, (TP + FN), method = "wilson")$upper,
                  rec_lower = binom.confint(TP, (TP + FN), method = "wilson")$lower) %>%
    dplyr::mutate(annot = in_annot, .before = TP)
}

