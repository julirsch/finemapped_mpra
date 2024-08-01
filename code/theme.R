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
 
# Cell type colors
ct_colors <- c(
  pnw_palette('Bay',8,type='continuous')[5],
  pnw_palette('Bay',8,type='continuous')[3],
  pnw_palette('Starfish')[5],
  pnw_palette('Bay',8,type='continuous')[2],
  pnw_palette('Bay',8,type='continuous')[8]
)
names(ct_colors) <- c("A549", "HCT116", "HEPG2", "K562", "SKNSH")
