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

# Read in DHS count matrix (from Meuleman)
dhs_gr <- readRDS("data/annotations/meuleman/meuleman_cpm_gr.rds")
dhs_small_gr <- dhs_gr[,0]

# Get DHS max across cell_types
values(dhs_small_gr) <- DataFrame(dhs_max = rowMaxs(dhs_gr@elementMetadata %>% as.matrix()))

# Convert MPRA variants to GRanges
variant_gr <- mpra_df %>%
  dplyr::distinct(chromosome, position, variant) %>%
  na.omit() %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "position", 
                           end.field = "position",
                           keep.extra.columns = T)

# Overlap MPRA variants with DHS
idx <- findOverlaps(variant_gr, dhs_small_gr)
dhsvars_df <- bind_cols(variant_gr[idx@from] %>%
                        as_tibble() %>%
                        dplyr::select(variant),
                      dhs_small_gr[idx@to] %>%
                        as_tibble(.name_repait = "minimal") %>%
                        dplyr::select(-seqnames, -start, -end, -width, -strand))

# Filter to non-coding, fine-mapped variants
out_df <- mpra_df %>%
  dplyr::filter(cs_id > 0, cohort != "control") %>%
  dplyr::filter(CRE_Meuleman == 1, consequence2 %ni% c("LoF", "missense", "synonymous")) %>%
  left_join(dhsvars_df,
            by = "variant") %>%
  dplyr::select(variant, log2FC_meta, dhs_max, promoter) %>%
  distinct() %>%
  na.omit()

# Get correlations between DHS and activity (Fig. 1d)
out_df %>% 
  dplyr::group_by(promoter) %>%
  dplyr::summarize(cor = Hmisc::rcorr(log2FC_meta, dhs_max, type = "pearson")$r[2],
                   pval = cor.test(log2FC_meta, dhs_max, type = "pearson")$p.value)

# Truncate outlier values
out_df <- out_df %>%
  dplyr::mutate(dhs_max = ifelse(dhs_max >= quantile(dhs_max, 0.99), quantile(dhs_max, 0.99), dhs_max),
                log2FC_meta = ifelse(log2FC_meta >= quantile(log2FC_meta, 0.99), quantile(log2FC_meta, 0.99), log2FC_meta))

# Plot correlations (Fig. 1d)
p1 <- out_df %>%
  ggplot(aes(x = dhs_max, y = log2FC_meta, group = promoter, color = as.factor(promoter))) +
  xlab("Chromatin accessibility (max)") + 
  ylab("Activity (Log2FC)") + 
  geom_smooth(method = "gam") +
  scale_color_manual(values = jdb_palette("FantasticFox")[3:4]) +
  pretty_plot() +
  plot_theme +
  theme(legend.position = "none")

# Save plot (Fig. 1d)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(2))
cowplot::save_plot(
  paste0("figures/fig1/fig1c_dhscorr.pdf"),
  plt_combined,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf
)

