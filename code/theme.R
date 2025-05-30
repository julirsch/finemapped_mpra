# Setup standard theme
plot_theme <-
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

# Define trait domain colors
domain_colors <- c("#4C72AE", "#E4812F", "#F3DE67", "#6BBCCC", "#5AA245", "#BA2E2C", "#8C61B6", "#80584E", "#CE72BB", "#BACD3C", "#7F7F7F", "black")
names(domain_colors) <- c("Metabolic", "Lipids", "Cardiovascular", "Immunological", "Hematopoietic", "Hepatic", "Renal", "Skeletal", "Neurological", "Behavioral", "Other", "Multiple")

# Definte tissue system colors
system_colors <- c("#86CB66", "#4C72AE", "#F3DE67", "#CE72BB", "#b565b1", "#6BBCCC", "#78a996", "#80584E", "#FD6467", "#80584E", "#FDD380", "black")
names(system_colors) <- c("Integumentary", "Endocrine", "Circulatory", "Nervous", "Exocrine", "Immune", "Digestive", "Renal", "Respiratory", "Muscular", "Reproductive", "Multiple")

# Define cell type colors
ct_colors <- c(
  pnw_palette('Bay',8,type='continuous')[5],
  pnw_palette('Bay',8,type='continuous')[3],
  pnw_palette('Starfish')[5],
  pnw_palette('Bay',8,type='continuous')[2],
  pnw_palette('Bay',8,type='continuous')[8]
)
names(ct_colors) <- c("A549", "HCT116", "HEPG2", "K562", "SKNSH")
