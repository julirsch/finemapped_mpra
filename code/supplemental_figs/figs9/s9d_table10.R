# Load in libraries
source("code/utils.R")
source('code/github_ready.prc.functions.R')
library(vroom)
library(ggplot2)
library(tidyverse)
options(stringsAsFactors = FALSE)

# Read in preprocessed eQTL MPRA data
gtex_mpra <- vroom('/mnt/sdb/gwas_eqtl_mpra/data/preprocess/gtex.txt.gz')

# make cell type combinations
one_g <- combn(c("K562","HEPG2","SKNSH","HCT116"),1,)
two_g <- combn(c("K562","HEPG2","SKNSH","HCT116"),2, FUN = paste, collapse = ',')
three_g <- combn(c("K562","HEPG2","SKNSH","HCT116"),3, FUN = paste, collapse = ',')
four_g <- combn(c("K562","HEPG2","SKNSH","HCT116"),4, FUN = paste, collapse = ',')

# prep output matrix
out <- data.frame(matrix(ncol = 20, nrow = 0))
colnames(out) <- c("cohort","annot","TP","FP","FN","TN","total_pos","total","FPR","N","FP_new","TN_new","precision","recall","prec_upper","prec_lower","rec_uper","rec_lower","ct","tt")

# iterate through combinations of 
for (cmb in c(one_g, two_g, three_g, four_g)){
  for (c in cmb){
    print(c)
    len = str_count(c, ',') + 1
    if (len == 1){
      cells = c
    }else {
      cells = str_split(c,pattern = ',', n = len)[[1]]
    }
    
    mpra <- gtex_mpra 
    
    # Filter to relevant cell types
    mpra <- mpra %>% dplyr::filter(cell_type %in% cells) %>%
      dplyr::select(-active_any, -emVar_any,-emVar_all) %>%
      dplyr::group_by(variant) %>%
      dplyr::mutate(active_any = ifelse(any(active), T, F),
                    emVar_any = ifelse(any(emVar), T, F),
                    emVar_all = ifelse(all(emVar),T,F)) %>%
      ungroup()
    
    # Add annotations
    mpra <- mpra %>% ungroup() %>% 
      add_annots("CRE_yes_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 1, CRE == 0), .) %>%
      add_annots("CRE_no_emVar_any_yes", rlang::exprs(emVar_any == 1, emVar_any == 0, CRE == 0, CRE == 1), .)
    
    # Make PrC set
    mpra_mc <- make_prset(mpra, "eQTLs")
    
    # Calculate precision and recall
    annotlist <- c("emVar_any", "CRE", "CRE_yes_emVar_any_yes", "CRE_no_emVar_any_yes")
    eqtl_out <- lapply(annotlist, function(x) {prec_rec(x, mpra_mc)}) %>%
      bind_rows() %>%
      dplyr::mutate(cohort = "eQTL", .before = annot) %>%
      dplyr::mutate(ct = len, tt = c)
    out <- rbind(out,eqtl_out)
  }
}

to_plot_out_gtex <- out

# Prepare for plotting
to_plot_out_gtex$ct <- as.character(to_plot_out_gtex$ct)

# Meta-analyze
meta_emvar_gtex <- to_plot_out_gtex %>%
  dplyr::filter(annot == 'emVar_any') %>%
  dplyr::mutate(pse2 = precision * (1-precision) / (TP +FP),
                rse2 = recall * (1-recall) / (TP +FN))%>%
  dplyr::select( -prec_upper, -prec_lower,-rec_upper, -rec_lower) %>%
  group_by(ct) %>% 
  dplyr::mutate(precision = mean(precision),
                recall = mean(recall),
                c = n()) %>%
  dplyr::mutate(pse2 = sqrt(sum(pse2))/c,
                rse2 = sqrt(sum(rse2))/c) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::mutate(tt = paste0('meta', ct),
                ct = 'meta') %>%
  dplyr::mutate(prec_upper = precision +pse2*1.96,
                prec_lower = precision - pse2*1.96,
                rec_upper = recall+ rse2*1.96,
                rec_lower = recall - rse2*1.96) %>%
  dplyr::select(-pse2, -rse2, -c) %>%
  ungroup()

meta_emvar_CRE_gtex <- to_plot_out_gtex %>%
  dplyr::filter(annot == 'CRE_yes_emVar_any_yes') %>%
  dplyr::mutate(pse2 = precision * (1-precision) / (TP +FP),
                rse2 = recall * (1-recall) / (TP +FN))%>%
  dplyr::select( -prec_upper, -prec_lower,-rec_upper, -rec_lower) %>%
  group_by(ct) %>% 
  dplyr::mutate(precision = mean(precision),
                recall = mean(recall),
                c = n()) %>%
  dplyr::mutate(pse2 = sqrt(sum(pse2))/c,
                rse2 = sqrt(sum(rse2))/c) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::mutate(tt = paste0('meta', ct),
                ct = 'meta') %>%
  dplyr::mutate(prec_upper = precision +pse2*1.96,
                prec_lower = precision - pse2*1.96,
                rec_upper = recall+ rse2*1.96,
                rec_lower = recall - rse2*1.96) %>%
  dplyr::select(-pse2, -rse2, -c) %>%
  ungroup()

# Combine
to_plot_out_gtex <- rbind(to_plot_out_gtex, meta_emvar_gtex, meta_emvar_CRE_gtex)

# fit a linear model
meta <- meta_emvar_CRE_gtex %>% dplyr::mutate(tt = as.numeric(str_split_fixed(meta_emvar_CRE_gtex$tt, 'meta', n=2)[,2]))
lp <- lm('precision ~ tt', data = meta)
lr <- lm('recall ~ tt', data = meta)
for (ct in seq(10)){
  print(paste0(ct,' celltypes'))
  print(paste0('Precision is: ', lp$coefficients[1]+lp$coefficients[2]*ct))
  print(paste0('Recall is: ', lr$coefficients[1]+lr$coefficients[2]*ct))
}

# get predictions
pred_prec <- lapply(seq(10), FUN = function(ct) lp$coefficients[1]+lp$coefficients[2]*ct)
pred_rec <- lapply(seq(10), FUN = function(ct) lr$coefficients[1]+lr$coefficients[2]*ct)
pred <- as.data.frame(cbind(pred_prec,pred_rec)) %>% dplyr::mutate(annot = 'CRE_yes_emVar_any_yes', cell_type = 'predicted')
pred$pred_prec <- as.numeric(pred$pred_prec)
pred$pred_rec <- as.numeric(pred$pred_rec)

# Save off parts of supplementary table 10
to_plot_out_gtex %>% dplyr::mutate(panel = 'FigS9 gtex cell type combo') %>% vroom::vroom_write(file = "tables/s9d.txt.gz", delim = "\t", col_names = T)
pred %>% dplyr::mutate(panel = 'FigS9 gtex cell type combo predictions') %>% vroom::vroom_write(file = "tables/s9d.pred.txt.gz", delim = "\t", col_names = T)

# Plot
p <- ggplot(data =  to_plot_out_gtex %>% dplyr::filter(annot == 'emVar_any'), aes(x = recall, y = precision, color=ct, shape = annot)) +
  geom_point(data = to_plot_out_gtex %>% dplyr::filter(annot == 'emVar_any',ct != 'meta'), size=3) +
  geom_errorbar(aes(ymin = prec_lower,ymax = prec_upper),width=0)+
  geom_errorbarh(aes(xmin = rec_lower,xmax = rec_upper),height=0)+
  geom_point(data = to_plot_out_gtex %>% dplyr::filter(ct == 'meta'), size = 1.5) +
  geom_line(data = to_plot_out_gtex %>% dplyr::filter(ct == 'meta')) +
  geom_point(data = to_plot_out_gtex %>% dplyr::filter(ct == 'meta',annot == 'CRE_yes_emVar_any_yes'), size = 3) +
  geom_errorbar(data = to_plot_out_gtex %>% dplyr::filter(ct == 'meta'), aes(ymin = prec_lower,ymax = prec_upper),width=0)+
  geom_errorbarh(data = to_plot_out_gtex %>% dplyr::filter(ct == 'meta'), aes(xmin = rec_lower,xmax = rec_upper),height=0)+
  BuenColors::pretty_plot(fontsize = 20)+
  theme(legend.position = 'right',
        aspect.ratio=1,
        panel.border = element_blank(),
        axis.line = element_line())+
  scale_color_manual(values = c(BuenColors::jdb_palette("brewer_spectra")[c(1,2,3,6)], 'black')) +
  scale_shape_manual(values=c(17,19))+
  geom_hline(yintercept = 0.5, color="grey50", linetype='dashed')+
  ylim(0.49,1)+xlim(0,0.28)+
  ggtitle('GTEx')+
  theme(legend.position = "none") 

p

# Save
library(patchwork)
plt_combined <- p + plot_layout(nrow = 1, heights = c(6))
cowplot::save_plot('figures/s9d.pdf',
                   #paste0("../tables/", pdfname, "_tracks.pdf"),
                   plt_combined,
                   base_height = 6,
                   base_width = 6,
                   device = cairo_pdf
)

