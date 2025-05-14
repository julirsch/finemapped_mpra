
source("code/utils.R")
library(PropCIs)
library(vroom)
library(tidyverse)
library(BuenColors)
library(ggplot2)
library(dplyr)
library(binom)
library(tidyr)
library(TFBSTools)
library(GenomicRanges)
options(stringsAsFactors = FALSE)

# read in mpra data
mpra_df <- vroom('data/preprocess/core_mpra.txt.gz') 

# Small collapsed MPRA
mpra_k562 <- mpra_df %>%
  dplyr::filter(cs_id > 0, cohort != "control") %>%
  distinct(variant, cell_type, active, active_any) %>%
  group_by(variant, cell_type) %>%
  dplyr::summarize(active = any(active),
                   active_any = any(active)) %>%
  ungroup() %>%
  dplyr::filter(cell_type == 'K562')

# Examining the utility of ChIP outside of tested cell type
# make variant granges
variant_gr <- mpra_k562 %>%
  distinct(variant) %>% 
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) 
variant_gr <- as.data.frame(variant_gr)
rownames(variant_gr) <- variant_gr$variant
variant_gr <- variant_gr %>% 
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "position", 
                           end.field = "position",
                           keep.extra.columns = F)

# read in ChIP-seq data
chip_gr <-readRDS("/mnt/sdb/ukbb-finemapping/data/chipatlas/ChIPAtlas.ChIP.rds")


# get only K562 ChIP-seq data
k562_gr <- chip_gr[which(elementMetadata(chip_gr)[,'celltype'] == 'K-562')]
k562_gr_list <- split(k562_gr, as.factor(k562_gr$antigen))
k562_df <- data.frame(lapply(1:length(k562_gr_list), 
                             function(x) {countOverlaps(variant_gr, k562_gr_list[[x]])}))
names(k562_df) <- paste("TF" ,names(k562_gr_list),sep=".")
k562_df[,1:190] <- ifelse(k562_df[,1:190] > 0, 1, 0)
k562_df$sumrow <- rowSums(k562_df)
k562_df <- k562_df %>% dplyr::mutate(sumrow = ifelse(sumrow > 0, 1, 0)) # binarize

# make a matrix without K562 , only subsetting to antigens in k562
no_k562_gr <- chip_gr[which(elementMetadata(chip_gr)[,'celltype'] != 'K-562' & 
                              elementMetadata(chip_gr)[,'antigen'] %in% names(k562_gr_list))]
no_k562_gr_list <- split(no_k562_gr, as.factor(no_k562_gr$antigen))
no_k562_df <- data.frame(lapply(1:length(no_k562_gr_list), 
                                function(x) {countOverlaps(variant_gr, no_k562_gr_list[[x]])}))
names(no_k562_df) <- paste("TF" ,names(no_k562_gr_list),sep=".")
dim(no_k562_df)
no_k562_df[,1:171] <- ifelse(no_k562_df[,1:171] > 0, 1, 0)
no_k562_df$sumrow <- rowSums(no_k562_df)
no_k562_df <- no_k562_df %>% dplyr::mutate(sumrow = ifelse(sumrow > 0, 1, 0)) # binarize

# in k562:
mpra_k562_in_chip_k562 <- mpra_k562 %>% 
  left_join(tibble::rownames_to_column(k562_df, 'variant') %>% 
              dplyr::rename(k562_chip_any = sumrow), 
            by='variant')  %>% 
  dplyr::select(-k562_chip_any) %>%
  pivot_longer(cols = starts_with('TF'),
               names_to = 'TF',
               names_prefix = 'TF.',
               values_to = 'in.chip') 


# get difference in proportion
k562_diff <- mpra_k562_in_chip_k562 %>% 
  dplyr::group_by(TF) %>%
  dplyr::mutate(num_in_chip = sum(in.chip > 0)) %>%
  dplyr::filter(num_in_chip > 20) %>%
  dplyr::select(-num_in_chip) %>%
  count(active,in_chip) %>% 
  #dplyr::filter(n_distinct(in.chip, active) == 4) %>% 
  ungroup() %>%
  dplyr::group_by(TF, active) %>%
  dplyr::mutate(num.in.group = sum(n)) %>%
  ungroup() %>%
  dplyr::filter(in.chip == 1) %>%
  dplyr::rename('num.in.chip' = 'n') %>% 
  pivot_wider(id_cols = TF, names_from = c(active), values_from = c(num.in.chip,num.in.group), names_sep = '.') %>%
  na.omit() %>%
  dplyr::group_by(TF) %>%
  dplyr::summarize(diff = PropCIs::wald2ci(num.in.chip.TRUE, num.in.group.TRUE, num.in.chip.FALSE, num.in.group.FALSE, conf.level = 0.95, adjust = "AC")$estimate,
                   diff_lower = PropCIs::wald2ci(num.in.chip.TRUE, num.in.group.TRUE, num.in.chip.FALSE, num.in.group.FALSE, conf.level = 0.95, adjust = "AC")$conf.int[1],
                   diff_upper = PropCIs::wald2ci(num.in.chip.TRUE, num.in.group.TRUE, num.in.chip.FALSE, num.in.group.FALSE, conf.level = 0.95, adjust = "AC")$conf.int[2],
                   num.in.chip.FALSE = num.in.chip.FALSE,
                   num.in.group.FALSE = num.in.group.FALSE,
                   num.in.chip.TRUE = num.in.chip.TRUE,
                   num.in.group.TRUE = num.in.group.TRUE)

# get proportions  
k562_prop <- mpra_k562_in_chip_k562 %>% 
  dplyr::group_by(TF) %>%
  dplyr::mutate(num_in_chip = sum(in.chip > 0)) %>%
  dplyr::filter(num_in_chip > 20) %>%
  dplyr::select(-num_in_chip) %>%
  count(active,in.chip) %>% 
  dplyr::filter(n_distinct(in.chip, active) == 4) %>%
  dplyr::summarise(p.val = prop.test(x=matrix(c(n[in.chip == 1 & active == T],
                                                n[in.chip == 1 & active == F],
                                                n[in.chip == 0 & active == T],
                                                n[in.chip == 0 & active == F]),
                                              nrow=2))$p.value,
                   prop.active = prop.test(x=matrix(c(n[in.chip == 1 & active == T],
                                                      n[in.chip == 1 & active == F],
                                                      n[in.chip == 0 & active == T],
                                                      n[in.chip == 0 & active == F]),
                                                    nrow=2))$estimate[[1]],
                   prop.inactive = prop.test(x=matrix(c(n[in.chip == 1 & active == T],
                                                        n[in.chip == 1 & active == F],
                                                        n[in.chip == 0 & active == T],
                                                        n[in.chip == 0 & active == F]),
                                                      nrow=2))$estimate[[2]])


# now, filter out any variant in k562 ChIP-seq peaks
mpra_no_k562 <- mpra_k562 %>%
  dplyr::filter(variant %ni% (tibble::rownames_to_column(k562_df, 'variant') %>%
                                dplyr::filter(sumrow != 0) %>% 
                                distinct(variant) %>% 
                                .$variant)) %>%
  left_join(tibble::rownames_to_column(no_k562_df, 'variant') %>% 
              dplyr::rename(nok562.chip.any = sumrow), 
            by='variant')  %>% 
  dplyr::select(-nok562.chip.any) %>%
  pivot_longer(cols = starts_with('TF'),
               names_to = 'TF',
               names_prefix = 'TF.',
               values_to = 'in.chip') 

# get difference in proportion
nok562_diff <- mpra_no_k562 %>% 
  dplyr::group_by(TF) %>%
  dplyr::mutate(num_in_chip = sum(in.chip > 0)) %>%
  dplyr::filter(num_in_chip > 20) %>%
  dplyr::select(-num_in_chip) %>%
  count(active,in.chip) %>% 
  #dplyr::filter(n_distinct(in.chip, active) == 4) %>% 
  ungroup() %>%
  dplyr::group_by(TF, active) %>%
  dplyr::mutate(num.in.group = sum(n)) %>%
  ungroup() %>%
  dplyr::filter(in.chip == 1) %>%
  dplyr::rename('num.in.chip' = 'n') %>% 
  pivot_wider(id_cols = TF, names_from = c(active), values_from = c(num.in.chip,num.in.group), names_sep = '.') %>%
  na.omit() %>%
  dplyr::group_by(TF) %>%
  dplyr::summarize(diff = PropCIs::wald2ci(num.in.chip.TRUE, num.in.group.TRUE, num.in.chip.FALSE, num.in.group.FALSE, conf.level = 0.95, adjust = "AC")$estimate,
                   diff_lower = PropCIs::wald2ci(num.in.chip.TRUE, num.in.group.TRUE, num.in.chip.FALSE, num.in.group.FALSE, conf.level = 0.95, adjust = "AC")$conf.int[1],
                   diff_upper = PropCIs::wald2ci(num.in.chip.TRUE, num.in.group.TRUE, num.in.chip.FALSE, num.in.group.FALSE, conf.level = 0.95, adjust = "AC")$conf.int[2],
                   num.in.chip.FALSE = num.in.chip.FALSE,
                   num.in.group.FALSE = num.in.group.FALSE,
                   num.in.chip.TRUE = num.in.chip.TRUE,
                   num.in.group.TRUE = num.in.group.TRUE)

nok562_prop <- mpra_no_k562 %>% 
  dplyr::mutate(num_in_chip = sum(in.chip > 0)) %>%
  dplyr::filter(num_in_chip > 20) %>%
  dplyr::select(-num_in_chip) %>%
  dplyr::group_by(TF) %>%
  count(active,in.chip) %>% 
  dplyr::filter(n_distinct(in.chip, active) == 4) %>%
  dplyr::summarise(p.val = prop.test(x=matrix(c(n[in.chip == 1 & active == T],
                                                n[in.chip == 1 & active == F],
                                                n[in.chip == 0 & active == T],
                                                n[in.chip == 0 & active == F]),
                                              nrow=2))$p.value,
                   prop.active = prop.test(x=matrix(c(n[in.chip == 1 & active == T],
                                                      n[in.chip == 1 & active == F],
                                                      n[in.chip == 0 & active == T],
                                                      n[in.chip == 0 & active == F]),
                                                    nrow=2))$estimate[[1]],
                   prop.inactive = prop.test(x=matrix(c(n[in.chip == 1 & active == T],
                                                        n[in.chip == 1 & active == F],
                                                        n[in.chip == 0 & active == T],
                                                        n[in.chip == 0 & active == F]),
                                                      nrow=2))$estimate[[2]])



# Prepare for plotting
joint <- rbind(right_join(nok562_prop, nok562_diff, by = 'TF') %>% dplyr::mutate(K562_chip_overlap = 'No'),
               right_join(k562_prop,k562_diff, by = 'TF') %>% dplyr::mutate(K562_chip_overlap = 'Yes')) %>%
  dplyr::mutate(num.in.chip = num.in.chip.FALSE + num.in.chip.TRUE)  %>%
  na.omit()

joint_pivot <- joint %>% 
  dplyr::select(-c(num.in.chip.FALSE, num.in.group.FALSE,num.in.chip.TRUE,num.in.group.TRUE)) %>%
  arrange(TF) %>% 
  pivot_wider(id_cols = 'TF', names_from = 'K562_chip_overlap',values_from = c(p.val, prop.active, prop.inactive, diff, diff_lower, diff_upper, num.in.chip)) %>% 
  arrange(-diff_Yes)
to_plot <- joint %>%
  group_by(K562_chip_overlap) %>%
  arrange(-diff) %>%
  dplyr::mutate(id = row_number()) %>%
  dplyr::mutate(important = ifelse(TF %in% c('GATA1',"TAL1","CEBPA","EP300","SPI1","HIF1A"),1,0))

# plot
p <-ggplot(to_plot,aes(x=id, group = K562_chip_overlap, color = K562_chip_overlap,y=diff, label = TF))+
  geom_point()+
  geom_text(data=to_plot[to_plot$important == 1,], angle = 90, hjust = 1
  ) +
  geom_errorbar(aes(ymin = diff_lower, ymax = diff_upper), width = 0)+
  xlab("TF")+
  ylab('Diff in proportion, active vs inactive')+
  scale_color_manual(values = c('grey','maroon'))+
  theme_minimal()+ #removes axes
  geom_hline(aes(yintercept = 0))+ #add y-axis where you want
  theme( axis.line.y = element_line(color="black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())


p

#save
cowplot::save_plot(paste0('figures/s3b.chip.pdf'),
                   p,
                   base_height = 6,
                   base_width = 12,
                   device = cairo_pdf
)



## save off important dfs
# saveRDS(k562_df, '/mnt/sdb/gwas_eqtl_mpra/reviews/data/k562.df.RDS')
# saveRDS(mpra_k562, '/mnt/sdb/gwas_eqtl_mpra/reviews/data/mpra.k562.RDS')
# saveRDS(mpra_no_k562, '/mnt/sdb/gwas_eqtl_mpra/reviews/data/mpra.no.k562.RDS')
# saveRDS(no_k562_df, '/mnt/sdb/gwas_eqtl_mpra/reviews/data/no.k562.df.RDS')



# ORs in general - k562s
test1 <- mpra_k562 %>% 
  left_join(tibble::rownames_to_column(k562_df, 'variant') %>% 
              dplyr::rename(k562.chip.any = sumrow), 
            by='variant')  %>% 
  dplyr::select(-k562.chip.any) %>%
  pivot_longer(cols = starts_with('TF'),
               names_to = 'TF',
               names_prefix = 'TF.',
               values_to = 'in.chip') %>%
  dplyr::group_by(TF) %>%
  dplyr::mutate(num_in_chip = sum(in.chip > 0)) %>%
  dplyr::filter(num_in_chip > 20) %>%
  dplyr::select(-num_in_chip) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(variant,active) %>%
  dplyr::summarize(in.chip = max(in.chip)) %>%
  dplyr::ungroup() %>%
  count(active,in.chip) %>%
  arrange(-in.chip, active) 

test1


prop.test(matrix(test1$n, nrow =2))


# OR in general - excluding variants in k562 chip-seq peaks, 
# now looking at chip-seq peaks for the same antigens but in non-k562 cell types
test2 <- mpra_k562 %>%
  dplyr::filter(variant %ni% (tibble::rownames_to_column(k562_df, 'variant') %>%
                                dplyr::filter(sumrow != 0) %>% 
                                distinct(variant) %>% 
                                .$variant)) %>%
  left_join(tibble::rownames_to_column(no_k562_df, 'variant') %>% 
              dplyr::rename(nok562.chip.any = sumrow), 
            by='variant')  %>% 
  dplyr::select(-nok562.chip.any) %>%
  pivot_longer(cols = starts_with('TF'),
               names_to = 'TF',
               names_prefix = 'TF.',
               values_to = 'in.chip') %>%
  dplyr::group_by(TF) %>%
  dplyr::mutate(num_in_chip = sum(in.chip > 0)) %>%
  dplyr::filter(num_in_chip > 20) %>%
  dplyr::select(-num_in_chip) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(variant,active) %>%
  dplyr::summarize(in.chip = max(in.chip)) %>%
  dplyr::ungroup() %>%
  count(active,in.chip) %>%
  arrange(-in.chip, active) 

prop.test(matrix(test2$n, nrow =2))


prop.test(matrix(mpra_k562 %>%
                   dplyr::filter(cell_type == 'K562', cohort == 'GTEx') %>%
                   dplyr::filter(variant %ni% (tibble::rownames_to_column(k562_df, 'variant') %>%
                                                 dplyr::filter(sumrow != 0) %>% 
                                                 distinct(variant) %>% 
                                                 .$variant)) %>%
                   left_join(tibble::rownames_to_column(no_k562_df, 'variant') %>% 
                               dplyr::select(variant,sumrow) %>%
                               dplyr::rename(nok562.chip.any = sumrow), 
                             by='variant')  %>%
                   count(active,nok562.chip.any) %>%
                   arrange(-nok562chip.any, active) %>%
                   .$n, nrow =2))





