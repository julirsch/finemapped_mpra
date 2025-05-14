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

# Get counts of variants for categories
counts_df <- mpra_df %>%  
  dplyr::mutate(name = case_when(cohort == "GTEx" ~ "eQTLs",
                                 cohort %in% c("UKBB", "BBJ") ~ "Complex Traits",
                                 type %in% c("loc_CS", "loc_PIP10", "annot_CS", "annot_PIP10",  "null_CS", "null_PIP10") ~ "Controls: Complex Traits",
                                 type %in% c("3tissue_locctrl", "3tissue_annotctrl", "49tissue_locctrl", "49tissue_annotctrl") ~ "Controls: eQTLs")) %>% 
  distinct(variant, name) %>%
  dplyr::mutate(count = 1) %>%
  pivot_wider(names_from = name, values_from = count) %>%
  replace(is.na(.), 0) %>% 
  dplyr::mutate(`Controls: Both` = ifelse(`Controls: eQTLs` == 1 & `Controls: Complex Traits` == 1, 1, 0)) %>%
  dplyr::mutate(`Controls: eQTLs` = ifelse(`Controls: Both` == 1, 0, `Controls: eQTLs`)) %>% 
  dplyr::mutate(`Controls: Complex Traits` = ifelse(`Controls: Both` == 1, 0, `Controls: Complex Traits`)) %>% 
  dplyr::select(-variant) %>% 
  dplyr::summarize(across(everything(), sum)) %>%
  t() %>%
  as_tibble(rownames = "name") %>%
  dplyr::rename("count" = "V1") %>%
  dplyr::mutate(name = factor(name, levels = c("Controls: Both",
                                               "Controls: Complex Traits",
                                               "Controls: eQTLs",
                                               "Complex Traits", 
                                               "eQTLs")))

# Define colors for plot
plot_colors <- c("grey50", "#B9B5C3", "#AAC5CF", "#5E4FA2", "#3F96B7")
names(plot_colors) <- c("Controls: Both", "Controls: Complex Traits", "Controls: eQTLs", "Complex Traits", "eQTLs")

# Plot proportion of variants by category (Fig. 1b)
p1 <- counts_df %>%
  ggplot(aes(fill = name, x = "1", y = count)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_color_manual(values = plot_colors, aesthetics = c("color", "fill")) +
  coord_flip() + 
  BuenColors::pretty_plot(fontsize = 20) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position = "none")
p1

# Save plot (Fig. 1b - single variant bar plot)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot(
  paste0("figures/fig1/fig1b_singlevar.pdf"),
  plt_combined,
  base_height = 3,
  base_width = 3,
  device = cairo_pdf
)

# Read in trait domain info
domain_df <- vroom::vroom("data/ukbb/UKBB_96traits_release1.traits") %>%
  left_join(domain_colors %>% 
              as_tibble(rownames = "domain"), 
            by = "domain") %>% 
  dplyr::select(trait, domain, "color" = value)

# Select high-PIP complex trait variants and add domain info
mpra_ukbb_df <- mpra_df %>% 
  dplyr::filter(cohort %in% c("UKBB", "BBJ"), pip > 0.5) %>%
  distinct(variant, trait, cohort) %>% 
  left_join(domain_df,
            by = "trait") %>%
  dplyr::mutate(domain = ifelse(trait == "Glaucoma", "Other", domain))

# Get domain(s) for each variant
counts_df <- mpra_ukbb_df %>% 
  distinct(variant, domain) %>%
  group_by(variant) %>%
  dplyr::mutate(domain = case_when(length(variant) > 1 ~ "Multiple",
                                   T ~ domain)) %>%
  ungroup() %>%
  distinct(variant, domain) %>%
  dplyr::count(domain) %>%
  dplyr::mutate(domain = factor(domain, levels = c("Multiple", "Hematopoietic", "Other", "Metabolic", "Hepatic", "Renal", "Lipids",
                                                   "Skeletal", "Cardiovascular", "Immunological", "Behavioral", "Neurological")))

# Plot domains (Fig. 1b  - domain pie chart)
p1 <- counts_df %>%
  ggplot(aes(fill = domain, x = "", y = n)) +
  geom_bar(position = "stack", stat = "identity") +
  coord_polar("y", start = 0) +
  scale_color_manual(values = domain_colors, aesthetics = c("color", "fill")) +
  theme_void()
p1

# Save plot (Fig. 1b  - domain pie chart)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot(
  paste0("figures/fig1/fig1b_domainpie.pdf"),
  plt_combined,
  base_height = 3,
  base_width = 3,
  device = cairo_pdf
)

# Define cohort colors
cohort_colors <- c("#5e4fa2" ,"#a298c4", "#C5BFDC")
names(cohort_colors) <- c("Both", "UKBB", "BBJ")

# Get cohort(s) for each variant
counts_df <- mpra_ukbb_df %>%
  distinct(variant, cohort) %>%
  group_by(variant) %>%
  dplyr::mutate(cohort = case_when(length(variant) > 1 ~ "Both",
                                   T ~ cohort)) %>%
  ungroup() %>%
  distinct(variant, cohort) %>%
  dplyr::count(cohort) 

# Plot cohort pie chart (Fig. 1b  - cohort pie chart)
p1 <- counts_df %>%
  ggplot(aes(fill = cohort, x = "", y = n)) +
  geom_bar(position = "stack", stat = "identity") +
  coord_polar("y", start = 0) +
  scale_color_manual(values = cohort_colors, aesthetics = c("color", "fill")) +
  theme_void()
p1

# Save plot (Fig. 1b  - cohort pie chart)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot(
  paste0("figures/fig1/fig1b_cohortpie.pdf"),
  plt_combined,
  base_height = 3,
  base_width = 3,
  device = cairo_pdf
)

# Select high-PIP eQTL variants
gtex_df <- mpra_df %>%
  dplyr::filter(cohort == "GTEx", pip > 0.5) %>%
  distinct(variant, pip, system)

# Get system(s) for each variant
counts_df <- gtex_df %>% 
  distinct(variant, system) %>%
  group_by(variant) %>%
  dplyr::mutate(system = case_when(length(variant) > 1 ~ "Multiple",
                                   T ~ system)) %>%
  ungroup() %>%
  distinct(variant, system) %>%
  dplyr::count(system) %>%
  dplyr::mutate(system = factor(system, levels = c("Multiple", "Nervous", "Digestive", "Integumentary", "Circulatory", "Reproductive",
                                                   "Endocrine", "Muscular", "Respiratory", "Immune", "Exocrine", "Renal")))

# Plot systems (Fig. 1b  - system pie chart)
p1 <- ggplot(counts_df, aes(fill = system, x = "", y = n)) +
  geom_bar(position = "stack", stat = "identity") +
  coord_polar("y", start = 0)+
  scale_color_manual(values = system_colors, aesthetics = c("color", "fill")) +
  theme_void()
p1

# Save plot (Fig. 1b  - system pie chart)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot(
  paste0("figures/fig1/fig1b_systempie.pdf"),
  plt_combined,
  base_height = 3,
  base_width = 3,
  device = cairo_pdf
)

# Read in MPRA haplotype data
haplos_df <- vroom("/mnt/sdb/gwas_eqtl_mpra/data/haplos/haplos_master_table.txt.gz")

# Define cohort colors
cohort_colors <- c("#88CFA4", "#5E4FA2", "#3F96B7")
names(cohort_colors) <- c("Both", "Complex traits", "eQTLs")

# Get cohort(s) for each haplotype
counts_df <- haplos_new_df %>%
  distinct(v1, v2, library, cell_type) %>%
  dplyr::mutate(cohort = ifelse(library %in% c("OL41", "OL42"), "eQTLs", "Complex traits")) %>%
  distinct(v1, v2, cohort) %>%
  group_by(v1, v2) %>%
  dplyr::mutate(cohort = case_when(length(v1) > 1 ~ "Both",
                                   T ~ cohort)) %>%
  ungroup() %>%
  distinct(v1, v2, cohort) %>%
  dplyr::count(cohort) %>%
  dplyr::mutate(cohort = factor(cohort, levels = c("Both", "Complex traits", "eQTLs")))

# Plot proportion of haplotypes by cohort (Fig. 1b - haplotype bar plot)
p1 <- counts_df %>%
  ggplot(aes(fill = cohort, x = "", y = n)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_color_manual(values = cohort_colors, aesthetics = c("color", "fill")) +
  coord_flip() + 
  BuenColors::pretty_plot(fontsize = 20) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position = "none")
p1

# Save plot (Fig. 1b - haplotype bar plot)
plt_combined <- p1 + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot(
  "figures/fig1/fig1b_haplos.pdf",
  plt_combined,
  base_height = 5,
  base_width = 10,
  device = cairo_pdf
)





# get controls PIP
bar_plot_ctrls <- mpra_df %>%  
  dplyr::filter(cohort == "control") %>% 
  dplyr::mutate(cohort = case_when(type %in% c("loc_CS", "loc_PIP10",
                                               "annot_CS", "annot_PIP10",
                                               "null_CS","null_PIP10") ~ "Controls: Complex Traits",
                                   type %in% c("3tissue_locctrl","3tissue_annotctrl",
                                               "49tissue_locctrl","49tissue_annotctrl") ~ "Controls: eQTLs")) %>% 
  distinct(variant, cohort) %>%
  dplyr::mutate(count = 1) %>%
  pivot_wider(names_from = cohort,
              values_from = count) %>%
  replace(is.na(.),0) %>% 
  dplyr::mutate(`Controls: Both` = ifelse(`Controls: eQTLs` == 1 &
                                            `Controls: Complex Traits` == 1,
                                          1, 0)) %>%
  dplyr::mutate(`Controls: eQTLs` = ifelse(`Controls: Both` == 1, 0, `Controls: eQTLs`)) %>% 
  dplyr::mutate(`Controls: Complex Traits` = ifelse(`Controls: Both` == 1, 0, `Controls: Complex Traits`)) %>% 
  dplyr::select(-variant) %>% 
  colSums()
bar_plot_ctrls <-   melt(bar_plot_ctrls) %>% 
  tibble::rownames_to_column() %>%
  dplyr::rename(cat = rowname,
                n = value)


# Test variants
# get PIP > 0.1, in CRE
pip_bin_breaks <- c(0, 0.1, 1.0)
bar_plot_pips <-  mpra_df %>% 
  dplyr::filter(cohort != "control") %>% 
  dplyr::mutate(cohort = case_when(cohort == "GTEx" ~ "eQTLs",
                                   cohort %in% c("UKBB", "BBJ") ~ "Complex Traits")) %>% 
  distinct(variant, pip, cohort, CRE) %>%
  group_by(variant, cohort) %>%
  dplyr::filter(pip == max(pip)) %>%
  ungroup() %>%
  dplyr::mutate(pip_bin = cut(pip, pip_bin_breaks, include.lowest = T, ordered_result = T)) %>%
  dplyr::mutate(cat = paste(cohort, pip_bin, CRE, sep = " ")) %>%
  distinct(variant, cat) %>%
  count(cat)

# combine for plotting
bar_plot_counts <- rbind(bar_plot_pips, bar_plot_ctrls)

# prepare to plot
bar_plot_counts$cat<- factor(bar_plot_counts$cat, levels = c("Controls: Both",
                                                             "Controls: Complex Traits",
                                                             "Controls: eQTLs",
                                                             "Complex Traits [0,0.1] 0",
                                                             "Complex Traits [0,0.1] 1",
                                                             "Complex Traits (0.1,1] 0",
                                                             "Complex Traits (0.1,1] 1",
                                                             "eQTLs [0,0.1] 0",
                                                             "eQTLs [0,0.1] 1",
                                                             "eQTLs (0.1,1] 0",
                                                             "eQTLs (0.1,1] 1"))

# define colors
colours <- c("grey50", # controls that overlap between eQTL controls and complex trait controls
             "#B9B5C3", # complex trait controls
             "#AAC5CF", # eQTL controls
             "#DEDBEC","#9E95C7","#7E72B4","#5E4FA2", # complex traits 
             "#BEDDE9","#9DCCDD","#7CBBD2","#3F96B7") # eQTLs

# plot 
p <- bar_plot_counts %>%
  dplyr::mutate(dummy = "dummy") %>%
  ggplot(., aes(fill = cat, x = dummy, y = n)) + 
  geom_bar(position="stack",stat="identity") + 
  scale_color_manual(values = colours,aesthetics = c("color","fill")) +
  coord_flip() + 
  BuenColors::pretty_plot(fontsize = 20)+
  xlab(NULL)+
  ylab(NULL)+
  labs(fill = "Cohort")+
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position="none")


p

# save
plt_combined <- p + plot_layout(nrow = 1, heights = c(3))
cowplot::save_plot("figures/s1a.haplo.pip.pdf",
                   plt_combined,
                   base_height = 5,
                   base_width = 10,
                   device = cairo_pdf
)




