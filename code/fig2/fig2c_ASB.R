# Load libraries
library(tidyverse)
library(vroom)
library(GenomicRanges)
library(rtracklayer)
library(plyranges)
library(ggrastr)
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

# Read in DHS allelic imbalance 
cnames <- c("chromosome", "start", "end", "ID", "ref", "alt",	
            "AAF", "RAF", "mean_BAD", "nSNPs", "max_cover",	"mean_cover", 
            "footprints_n",	"hotspots_n",	"group_id",	"pval_ref_combined",	"pval_alt_combined",	
            "es_combined", "logit_es_combined",	"min_pval",	"min_fdr")
dhs_asb_gr <- vroom("data/asb/all.aggregation.dnase.bed.gz", delim = "\t", col_names = cnames, skip = 1) %>%
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chromosome", 
                           start.field =  "start", 
                           end.field = "end", 
                           keep.extra.columns = T) 

# Liftover DHS allelic imbalnce to hg19
ch <- import.chain("resources/liftover/hg38ToHg19.over.chain")
dhs_asb_gr <- liftOver(dhs_asb_gr, ch) %>%
  unlist() %>%
  plyranges::mutate(var = paste0(seqnames, ":", end),
                variant = paste0(seqnames, ":", end, ":", ref, ":", alt))
dhs_asb_gr <- dhs_asb_gr[!(duplicated(dhs_asb_gr$variant) | duplicated(dhs_asb_gr$variant, fromLast = TRUE)),]
names(dhs_asb_gr) <- dhs_asb_gr$variant
dhs_asb_df <- dhs_asb_gr %>%
  as_tibble()

# Merge with MPRA data
out_df <- mpra_df %>%
  left_join(dhs_asb_df %>%
              #dplyr::filter(nSNPs > 10, max_cover > 100) %>%
              dplyr::filter(nSNPs > 20) %>%
              dplyr::mutate(caAI_fdr = min_fdr, caAI_es = -1 * logit_es_combined, caAI_thresh = ifelse(caAI_fdr < 0.25, T, F)) %>%
              dplyr::select(variant, caAI_fdr, caAI_es, caAI_thresh),
            by = "variant")

# Filter to non-coding, fine-mapped variants
out_df <- out_df %>% 
  dplyr::filter(cs_id > 0, cohort %in% c("UKBB", "GTEx", "BBJ"),
                consequence2 %ni% c("LoF", "missense", "synonymous")) %>%
  dplyr::select(variant, emVar_meta, emVar_any, log2Skew_meta, caAI_es, caAI_fdr, caAI_thresh) %>%
  distinct() %>% 
  na.omit()

# Write out variants used in analysis
out_df %>%
  distinct(variant) %>%
  data.table::fwrite("data/asb/all.aggregation.dnase.variants.txt", sep = "\t", col.names = T)

# Correlations between allelic imbalance and MPRA allelic effect
out_df %>%
  dplyr::group_by(caAI_thresh, emVar_any) %>%
  dplyr::summarize(n = length(variant),
                   cor = Hmisc::rcorr(log2Skew_meta, caAI_es, type = "pearson")$r[2],
                   pval = cor.test(log2Skew_meta, caAI_es, type = "pearson")$p.value)

# Plot correlation (Fig. 2c)
p1 <- out_df %>%
  dplyr::filter(caAI_thresh == T) %>%
  ggplot(aes(x = caAI_es, y = log2Skew_meta, color = emVar_any)) +
  xlab("Chromatin accessibility") + 
  ylab("Allelic effect") + 
  rasterize(geom_point(size = 0.5, alpha = 0.3), dpi = 1000) +
  geom_smooth(data = out_df %>% dplyr::filter(caAI_thresh == T, emVar_any == F), aes(x = caAI_es, y = log2Skew_meta), method = "lm", color = PNWColors::pnw_palette('Sunset2')[2], fill = PNWColors::pnw_palette('Sunset2')[2], se = T) +
  geom_smooth(data = out_df %>% dplyr::filter(caAI_thresh == T, emVar_any == T), aes(x = caAI_es, y = log2Skew_meta), method = "lm", color = PNWColors::pnw_palette('Sunset2')[4], fill = PNWColors::pnw_palette('Sunset2')[4], se = T) +
  scale_color_manual(values = c(PNWColors::pnw_palette('Sunset2')[2], PNWColors::pnw_palette('Sunset2')[4])) +
  pretty_plot() +
  coord_cartesian(xlim = c(-2.1, 2.1), ylim = c(-2.5, 2.5)) +
  plot_theme +
  theme(legend.position = "none")
p1

# Save plot (Fig. 2c)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(2))
cowplot::save_plot(
  paste0("figures/fig2/fig2c_caAI.pdf"),
  plt_combined,
  base_height = 2,
  base_width = 2,
  device = cairo_pdf
)


