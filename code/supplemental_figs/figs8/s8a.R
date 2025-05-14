
# Load in libraries
source("code/utils.R")
source("code/preprocess/read_in_functions.R")
library(vroom)
library(tidyverse)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)

# Read in data without best library selected
cohort = 'gtex'
paired = 'paired'

path <- "/mnt/sdb/gtex_mpra/final_mpra_data/processed/"
print(paste0(path,cohort,"_mpra_",paired,"_final20230117.rds"))
mpra <- readRDS(paste0(path,'gtex',"_mpra_",'paired',"_final20230117.rds")) %>% 
  dplyr::filter(type != "other_test") %>%
  dplyr::rename(log2Skew = Log2Skew)  %>% 
  dplyr::mutate(cs_uid = case_when(cohort == 'control' ~ NA_character_,
                                   cs_id == -1 ~ NA_character_,
                                   TRUE ~ paste0(cohort, ";", tissue, ";", gene, ";", cs_id))) 

cohort = 'traits'
paired = 'paired'

path <- "/mnt/sdb/gtex_mpra/final_mpra_data/processed/"
print(paste0(path,cohort,"_mpra_",paired,"_final20230117.rds"))
traits_mpra_df <- readRDS(paste0(path,cohort,"_mpra_",paired,"_final20230117.rds")) %>% 
  dplyr::filter(type != "other_test") %>% 
  dplyr::filter(cell_type != "GM12878") %>%
  dplyr::rename(log2Skew = Log2Skew) %>% 
  dplyr::mutate(cs_uid = case_when(cohort == 'control' ~ NA_character_,
                                   cs_id == -1 ~ NA_character_,
                                   TRUE ~ paste0(cohort, ";", trait, ";", region, ";", cs_id)))

mpra_df <- bind_rows(traits_mpra_df, mpra) %>%
  ungroup()
rm(traits_mpra_df); rm(mpra)

mpra_df <- def_emVars(mpra_df)
mpra_df <- mpra_df %>% ungroup()

# get background rates of emVars per library
# using low-pip CS (non-control) variants across GTEx, UKBB, and BBJ
mpra_spec <- mpra_df %>%
  dplyr::filter(cohort %in% c('GTEx',"UKBB","BBJ"), !grepl('ctrl',type), pip < 0.01) %>%
  distinct(variant, cell_type, library, emVar) %>% 
  dplyr::count(cell_type, library, emVar)%>% 
  group_by(cell_type, library) %>%
  dplyr::mutate(tot = sum(n)) %>%
  dplyr::filter(emVar == TRUE) %>%
  dplyr::select(-emVar) %>%
  dplyr::mutate(prop = n / tot,
                lower = binom.confint(n, tot, method='wilson')$lower,
                upper = binom.confint(n, tot, method='wilson')$upper) %>%
  arrange(-prop) 

# mean plasmid ref and alt counts
mpra_rna <- mpra_df %>%
  dplyr::filter(cohort %in% c('GTEx',"UKBB","BBJ"), !grepl('ctrl',type), pip < 0.01) %>%
  distinct(variant, mean_RNA_ref, mean_RNA_alt,cell_type, library, emVar) %>% 
  dplyr::mutate(avg_rna = (mean_RNA_ref + mean_RNA_alt)/2) %>%
  dplyr::group_by(cell_type,library) %>%summarize(avg = mean(avg_rna))

# add standard error
#mpra_spec <- mpra_spec %>% dplyr::mutate(se = sqrt((prop * (1-prop))/tot))
#mpra_spec<- mpra_spec %>% dplyr::mutate(se = max(abs(upper-prop)/1.96, abs(lower-prop)/1.96))

mpra_spec <- left_join(mpra_spec,mpra_rna, by = c("cell_type","library"))

# generate heatplot
p<- ggplot(mpra_spec, aes(cell_type, forcats::fct_rev(library), fill = log2(avg), size =prop)) +
  geom_point(shape = 21, stroke = 0) +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(1, 15), breaks = c(0.02, 0.04, 0.06, 0.08 )) +
  scale_fill_gradient(high = BuenColors::jdb_palette('brewer_spectra')[1], low = "#F9F5FF") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25),
                             label.position = "bottom",
                             title.position = "right",
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  labs(size = "Proportion emVars ", fill = "mean log2(RNA count)", x = NULL, y = NULL)
p

plt_combined <- p

plt_combined

cowplot::save_plot(
  'figures/figs2/libraries_heatmap.pdf',
  plt_combined,
  base_height = 5.8,
  base_width = 6,
  device = cairo_pdf)

