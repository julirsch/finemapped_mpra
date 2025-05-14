# Load in libraries
source("code/utils.R")
library(tidyverse)
library(BuenColors)
library(reticulate)
library(binom)
library(ggplot2)
library(cowplot)
library(patchwork)
library(vroom)
library(GenomicRanges)
library(ComplexHeatmap)
library(ggpubr)
options(stringsAsFactors = FALSE)

# read in processed MPRA data
mpra <- vroom('/mnt/sdb/gwas_eqtl_mpra/data/preprocess/core_mpra.txt.gz') %>%
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0)) 

# filter down to inactive variants
mpra.df <- mpra %>% 
  dplyr::filter(cohort != 'control') %>%
  dplyr::filter(active_any == F,
                pip > 0.5)

variant.gr <- mpra.df %>%
  distinct(variant, chromosome, position) %>%
  dplyr::mutate(chromosome = str_split_fixed(variant, ":", 4)[, 1],
                position = as.numeric(str_split_fixed(variant, ":", 4)[, 2])) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "position", 
                           end.field = "position",
                           keep.extra.columns = T)



# get annotations
atac <- readRDS("/mnt/sdb/ukbb-finemapping/data/annotations/meuleman/meuleman_cpm_gr.rds")
df <- as.data.frame(countOverlaps(atac,variant.gr)) %>%
  dplyr::rename("overlap" = "countOverlaps(atac, variant.gr)")
df <- tibble::rownames_to_column(df, "identifier")


# read in metadata from https://www.meuleman.org/DHS_Index_and_Vocabulary_metadata.html
metadata <- vroom::vroom("/mnt/sdb/gwas_eqtl_mpra/reviews/data/meuleman_dhs_metadata.txt") %>%
  dplyr::rename('biosample' = 'Biosample name',
                'altius' = "Altius Aggregation ID",
                'sample.type'= "Biosample type",
                'state'="Biological state")


# join in atac peak overlaps with peak information
test <- as.data.frame(atac) %>%
  tibble::rownames_to_column(var = 'identifier')%>% 
  dplyr::filter(identifier %in% (df %>% dplyr::filter(overlap >0) %>% .$identifier)) %>%
  left_join(df, by = 'identifier') %>%
  dplyr::select(-c('seqnames','start','start','end','width','strand'))

choose <- test %>% 
  pivot_longer(cols = colnames(test %>% dplyr::select(-identifier,-overlap )),
               names_to ='altius',
               values_to = 'cpm'
  ) %>%

  dplyr::filter(cpm > 0) %>%
  group_by(altius) %>%
  summarise(n = sum(overlap)) %>%
  ungroup() %>%
  left_join(metadata %>% distinct(Organ, altius, sample.type, biosample),
            by = 'altius') %>%
  arrange(-n) %>% print(n=500)

# add a column for all
l <- c('all',sum(df$overlap),'all',NA,'all')
names(l) <- c('altius','n','Organ','sample.type','biosample')  
choose <- rbind(l,choose)
choose$n <- as.numeric(choose$n)

# convert to proportions
choose <- choose %>% 
  dplyr::mutate(tot = dim(mpra.df %>% distinct(variant))[1]) %>%
  dplyr::mutate(prop = n/tot)

# define colors
clrs <- c('#5AA245','#BA2E2C','#FDD380',"#78a996","#4c6db8","#8C61B6","#CE72BB","#ae49bb",'firebrick',"#6bcde1","#FD6467","#222957","#80584E","#4C72AE","#fadc13",'black')
names(clrs) <- c('Hematopoietic','Hepatic','Genitourinary',"Digestive","Nervous","Epithelial","Integumentary","Connective","Cardiovascular","Embryonic","Respiratory", "Renal","Musculoskeletal","Endocrine","Fetal Life Support",NA)

p1 <- choose %>% dplyr::mutate(rn = row_number()) %>%
  left_join(metadata %>% distinct(altius, Organ, System), by = c('altius','Organ') )%>%
  dplyr::filter(rn <52) %>%
  ggplot(.,aes(x=rn, y = prop, fill = System))+
  geom_col()+
  pretty_plot() +
  scale_fill_manual(values = clrs) +
  #theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  #geom_errorbar(position = 'dodge', aes(ymin = lower, ymax = upper))+
  labs(fill = 'Organ system')+
  xlab('DHS sample')+
  ylab('Proportion inactive variants')

p1

# save
cowplot::save_plot(paste0('figures/s16b.pdf'),
                   #paste0("../tables/", pdfname, "_tracks.pdf"),
                   p1,
                   base_height = 4,
                   base_width = 6,
                   device = cairo_pdf
)

# plot out legend annotation

clrs2 <- c("cyan","yellow","magenta") 
names(clrs2) <- c("Cancer",'Lines',"Primary")

p2<- choose %>% dplyr::mutate(rn = row_number()) %>%
  dplyr::mutate(val= 1) %>%
  dplyr::filter(rn <  52) %>%
  ggplot(., aes(fill = sample.type, x = rn, y = val)) + 
  geom_col() + 
  scale_fill_manual(values = clrs2) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  pretty_plot()
p2

# save off legend
cowplot::save_plot(paste0('figures/s16b.legend.pdf'),
                   #paste0("../tables/", pdfname, "_tracks.pdf"),
                   p2,
                   base_height = 4,
                   base_width = 6,
                   device = cairo_pdf
)

# save off supplementary table
choose %>% vroom::vroom_write('tables/table27.txt')

# test for enricihment in top 50
choose %>% head(n=51) %>% left_join(metadata %>% distinct(altius, Organ, System), by = c('altius','Organ') )%>% count(System)
metadata %>% count(System) %>% print(n=500)
dhyper(16,91, 734-91,50)
