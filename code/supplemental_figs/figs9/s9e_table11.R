
# Load in libraries
library(vroom)
library(ggplot2)
library(tidyverse)
library(BuenColors)
options(stringsAsFactors = FALSE)


# read in data - already-picked causal variants
causal_tables <- vroom("tables/stable9.txt")
mpra_mc <- causal_tables %>% dplyr::filter(var_type %in%c("eQTL control", "eQTL test")) %>% dplyr::select(-var_type)
pip9 <- causal_tables %>% dplyr::filter(var_type %in%c("Complex traits control", "Complex traits test")) %>% dplyr::select(-var_type)

#precision and recall for emVar and emVar + CRE

## background rate * 
## precision is true positive / total positive

#TP/(TP + FN)

## background = 1/cs size 

# First, looking at emVar status
in_annot = 'emVar_any'

# calculate precision and recall
out_df <- mpra_mc %>%
  dplyr::select(!! sym(in_annot), causal) %>%
  dplyr::mutate(annot = case_when(!! sym(in_annot) > 0 ~ 1,
                                  !! sym(in_annot) == 0 ~ 0)) %>%
  dplyr::filter(!is.na(!! sym(in_annot))) %>%
  dplyr::count(causal, annot)

out_df <- out_df %>%
  dplyr::summarize(TP = sum(n[causal == T & annot == T]),
                   FP = sum(n[causal == F & annot == T]),
                   FN = sum(n[causal == T & annot == F]),
                   TN = sum(n[causal == F & annot == F]),
                   total_pos = sum(n[annot == T]),
                   total = sum(n)) 

# Scale to 1:1
if (dim(mpra_mc %>% dplyr::filter(causal == T))[1] < dim(mpra_mc %>% dplyr::filter(causal == F))[1]){
  out_df <- out_df  %>%
    dplyr::mutate(FPR = FP / (FP + TN),
                  N = TP + FN,
                  FP_new = FPR * N,
                  TN_new = N - FP_new) %>%
    dplyr::mutate(FP = FP_new,
                  TN = TN_new) 
  # else if there are more positives
} else if (dim(mpra_mc %>% dplyr::filter(causal == T))[1] > dim(mpra_mc %>% dplyr::filter(causal == F))[1]){
  out_df <- out_df  %>%
    dplyr::mutate(TPR = TP / (TP + FN),
                  N = TN + FP,
                  TP_new = TPR * N,
                  FN_new = N - TP_new) %>%
    dplyr::mutate(TP = TP_new,
                  FN = FN_new) 
}

# define parameters for precision and recall
TP <- out_df$TP
FP <- out_df$FP
FN <- out_df$FN

# predict
bg = c()
emvar = c()
for (i in seq(2,20)){
  bg= append(bg, 1/i)
  emvar = append(emvar,TP/(TP + (i-1)*FP))
}

# repeat for intersection of CRE and emVar status
in_annot = 'CRE_yes_emVar_any_yes'
out_df <- mpra_mc %>%
  dplyr::select(!! sym(in_annot), causal) %>%
  dplyr::mutate(annot = case_when(!! sym(in_annot) > 0 ~ 1,
                                  !! sym(in_annot) == 0 ~ 0)) %>%
  dplyr::filter(!is.na(!! sym(in_annot))) %>%
  dplyr::count(causal, annot)

out_df <- out_df %>%
  dplyr::summarize(TP = sum(n[causal == T & annot == T]),
                   FP = sum(n[causal == F & annot == T]),
                   FN = sum(n[causal == T & annot == F]),
                   TN = sum(n[causal == F & annot == F]),
                   total_pos = sum(n[annot == T]),
                   total = sum(n)) 

if (dim(mpra_mc %>% dplyr::filter(causal == T))[1] < dim(mpra_mc %>% dplyr::filter(causal == F))[1]){
  out_df <- out_df  %>%
    dplyr::mutate(FPR = FP / (FP + TN),
                  N = TP + FN,
                  FP_new = FPR * N,
                  TN_new = N - FP_new) %>%
    dplyr::mutate(FP = FP_new,
                  TN = TN_new) 
  # else if there are more positives
} else if (dim(mpra_mc %>% dplyr::filter(causal == T))[1] > dim(mpra_mc %>% dplyr::filter(causal == F))[1]){
  out_df <- out_df  %>%
    dplyr::mutate(TPR = TP / (TP + FN),
                  N = TN + FP,
                  TP_new = TPR * N,
                  FN_new = N - TP_new) %>%
    dplyr::mutate(TP = TP_new,
                  FN = FN_new) 
}


TP <- out_df$TP
FP <- out_df$FP
FN <- out_df$FN
cre_emvar = c()
for (i in seq(2,20)){
  cre_emvar = append(cre_emvar,TP/(TP + (i-1)*FP))
}

# prepare for plotting
out2 <- data.frame(cbind("CRE_emVar" = cre_emvar,
                         "emVar" = emvar,
                         "bkg" = bg,
                         "size" = seq(2,20)))

# Plot
sp <- ggplot(out2)+
  geom_line(aes(y=CRE_emVar, x = size), color =  BuenColors::jdb_palette("corona")[3])+
  geom_line(aes(y=emVar, x = size), color =  BuenColors::jdb_palette("corona")[1])+
  geom_line(aes(y=bkg, x = size), color = 'grey')+
  pretty_plot()
sp

# Save off supplementary table 11
out2 %>% vroom::vroom_write(file = "tables/stable11.txt.gz", delim = "\t", col_names = T)

# save plot
library(patchwork)
plt_combined <- sp + plot_layout(nrow = 1, heights = c(6))
cowplot::save_plot('figures/s9e.pdf',
                   plt_combined,
                   base_height = 6,
                   base_width = 6,
                   device = cairo_pdf
)
