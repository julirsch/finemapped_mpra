# Load libraries
library(tidyverse)
library(vroom)
library(GenomicRanges)
library(rtracklayer)
library(plyranges)
library(matrixStats)
library(annotables)
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

# Filter to core trait and control variants
mpra_df <- mpra_df %>% 
  dplyr::filter(cohort == "control" | (pip > 0.1 | cs_id > 0))

# Collapse variants across libraries, preferring emVars with more RNA counts
overview_df <- mpra_df %>% 
  dplyr::filter(mean_RNA_alt > 20, mean_RNA_ref > 20) %>%
  group_by(variant, cell_type) %>% 
  arrange(-emVar, -pmin(mean_RNA_alt, mean_RNA_ref)) %>% 
  dplyr::filter(row_number() == 1) %>% 
  ungroup() %>%  
  dplyr::mutate(name = paste(variant, cell_type, cohort, sep  = ";")) %>%
  distinct(variant, cell_type, library, log2FC, log2Skew, emVar, name, active_any, emVar_any, logPadj_BF, Skew_logPadj) 

# Collapse variants across cell-types, preferring emVars
overview_df <- overview_df %>%
  group_by(variant) %>%
  arrange(-emVar, -logPadj_BF, -Skew_logPadj) %>% 
  dplyr::filter(row_number() == 1) %>%
  ungroup()

# Plot activity and allelic effect for representative variants (Fig. 1f)
p1 <- overview_df %>%
  ggplot(aes(x = log2FC, y = log2Skew, color = emVar)) + 
  rasterize(geom_point(size = 0.2, alpha = 0.2), dpi = 1000) +
  rasterize(geom_point(data = overview_df %>% dplyr::filter(emVar == FALSE), aes(x = log2FC, y = log2Skew), size = 0.2, alpha = 0.2), dpi = 1000) +
  rasterize(geom_point(data = overview_df %>% dplyr::filter(emVar == TRUE), aes(x = log2FC, y = log2Skew), size = 0.2, alpha = 0.2), dpi = 1000) +
  pretty_plot() + 
  theme(legend.position = "bottom") +  
  scale_color_manual(values = c(pnw_palette("Sunset2")[2], pnw_palette("Sunset2")[4])) +
  guides(colour = "none") +
  theme(legend.title = element_blank()) +
  ylab("Allelic effect") +
  xlab("Activity")

# Add density plots
p1 <- ggMarginal(p1, groupColour = TRUE, groupFill = TRUE)

# Save plot (Fig. 1f)
plt_combined <- patchwork::wrap_elements(p1) + plot_layout(nrow = 1, heights = c(4))
cowplot::save_plot(
  paste0("/mnt/sdb/gwas_eqtl_mpra/figures/fig1/fig1f_scatter.pdf"),
  plt_combined,
  base_height = 4,
  base_width = 4.5,
  device = cairo_pdf
)
