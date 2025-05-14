library(tidyverse)
library(collapse)

# Set parameters
outdir <- "figures/"
code_path<-""

source(paste0(code_path,"code/Preprocess/utils.R"))

# Load motifs
#load(paste0(code_path,"data/preprocess/motifs_processed_20230225.RData"))
#motifs_pfml <- c(hocomoco_pfml, jaspar_pfml, cisbp_pfml, vierstra_pfml)

# Read in MPRA data
satmut_emVar_df <- vroom::vroom(paste0(code_path,"data/preprocess/satmut_emvars_table.txt"))

# Setup theme
my_theme <-
  BuenColors::pretty_plot(fontsize = 8) +
  BuenColors::L_border() +
  theme(
    plot.background = element_blank(),
    plot.margin = margin(0, 0.1, 0, 0.1, unit = "cm"),
    plot.title = element_text(hjust = 4e-3, margin = margin(b = -12)),
    # legend.position = "none",
    legend.position = c(0.3, 1),
    legend.justification = c(1, 1),
    legend.title = element_text(margin = margin(0, 0, 0, 0)),
    legend.background = element_blank(),
    legend.key.size = unit(0.2, "cm")
  )
# Basic GB numbers for emVars
satmut_emVar_df %>%
  collapse::replace_NA(0) %>%
  dplyr::summarize(mean(in_gb_or_gbmotif))
satmut_emVar_df %>%
  collapse::replace_NA(0) %>%
  dplyr::group_by(occupied_disruptmotif != 0) %>% 
  dplyr::summarize(mean(in_gb_or_gbmotif))

# Plot proportion of mechanisms that explain SatMut variants
in_plot <- satmut_emVar_df %>%
  collapse::replace_NA(0) %>%
  dplyr::filter(variant != "chr6:160560897:CTGGTAAGT:C") %>%
  dplyr::mutate(type = ifelse(occupied_disruptmotif == 0, "unknown", "canonical")) %>%
  #dplyr::mutate(in_occupied_disruptweakmotif = ifelse(occupied_disruptweakmotif == 0, 0, 1)) %>%
  dplyr::select(variant, type, in_fp, enf_chip, in_gb_or_gbmotif) %>%
  pivot_longer(cols = -c("variant", "type")) %>%
  dplyr::group_by(type, name) %>%
  dplyr::summarize(n = sum(value),
                   tot = length(value)) %>% 
  ungroup() %>% 
  #dplyr::group_by(name) %>%
  #dplyr::mutate(tot = sum(tot)) %>%
  dplyr::mutate(type = ordered(type, levels = c("unknown", "canonical"))) %>%
  arrange(desc(type)) %>%
  group_by(name) %>%
  dplyr::mutate(prop = n / tot,
                lower = binom.confint(n, tot, method='wilson')$lower,
                upper = binom.confint(n, tot, method='wilson')$upper,
                pvalue = binom.confint(n, tot, method='wilson')$p.value)
p1 <- in_plot %>%
  ggplot(., aes(x = fct_relevel(name, c("in_fp", "enf_chip", "in_gb_or_gbmotif")), 
                y = prop, fill = type)) + 
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                size = 0.5,
                width = 0,
                col = "black",
                position = position_dodge(0.9, preserve = "single")) +
  pretty_plot() +
  scale_fill_manual(values = BuenColors::jdb_palette("brewer_spectra")[c(2,6)]) +
  coord_cartesian(ylim = c(0.3, 1)) +
  my_theme +
  ylab("Proportion of emVars") +
  xlab("")
prop.test(c(39, 52), c(50, 52))

yes
