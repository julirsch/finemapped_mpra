library(reshape2)
library(ggplot2)
library(GGally)
library(dplyr)
library(tidyr)
library(stringr)
library(ggrastr)

#Code for various QC plots. 
#Code for Supplemental Figure 2c is written as OL32_rep_norm_count_log2_r123.png
#Code for supplemental figure 2d is written as Positive_controls_correlation_All.pdf
#Code for supplemental figure 2e is written as supplemental_1c_rna.pdf

#Folder containing input data files gtex_mpra_20230117.txt.gz, haplos_merged_20240202.txt, traits_mpra_20230117.txt.gz
filedir="/Users/tewher/Google_Drive_JAX/Manuscripts/2023_UKBB-BBJ-GTEx_MPRA/2023_UKBB_GTEx_MPRA_Manuscript/data/"

#Output folder
dir="/Users/tewher/Google_Drive_JAX/Manuscripts/2023_UKBB-BBJ-GTEx_MPRA/2023_UKBB_GTEx_MPRA_Manuscript/Analysis/Ryan/QC_figures/"
setwd(dir)

cor_func <- function(data, mapping, method, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, method=method, use='complete.obs')
  
  ggally_text(
    label = as.character(round(corr, 3)), 
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black',
    ...
  ) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
}

raster_points <- function(data, mapping, ...){
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  p <- ggplot(data) + geom_point_rast(size=0.01, alpha=0.1,raster.dpi=1000) + theme_minimal()
  p
}

########################################################
## Supplemental Figure 1 (MA Plot) and positive control correlations

activity_locs <- read.delim("activity.files", stringsAsFactors = F, header = F)
act_lists <- list()

for(i in 1:nrow(activity_locs)){
  lib <- activity_locs[i,2]
  cell <- activity_locs[i,3]
  message(paste0(lib," ",cell))
  act_lists[[lib]][[cell]] <- read.delim(paste0(filedir,activity_locs[i,1]), stringsAsFactors = F)
}

act_count_comp <- data.frame()
for(lib in unique(activity_locs[,2])){
  message(lib)
  for(cell in names(act_lists[[lib]])){
    if(cell=="GM12878") next
    message(cell)
    #temp <- act_lists[[lib]][[cell]][which(act_lists[[lib]][[cell]]$DNA_mean>30),c("ID","ctrl_mean","log2FoldChange","pvalue")]
    #temp <- temp[which(temp$ID %in% rownames(exp_means[[lib]])[exp_means[[lib]][,tolower(cell)]>0]),]
    temp <- act_lists[[lib]][[cell]][,c("ID","SNP","DNA_mean","ctrl_mean","exp_mean","log2FoldChange","pvalue")]
    temp$adj <- p.adjust(temp$pvalue, method="bonferroni")
    temp$library <- lib
    temp$celltype <- cell
    temp$sig <- "Not Significant"
    temp$sig[which(temp$adj < 0.01 & abs(temp$log2FoldChange) > 1)] <- "Active"
    temp$sig[temp$ID %in% act_lists[[lib]][[cell]]$ID[act_lists[[lib]][[cell]]$ctrl_mean < 30]] <- "Failed QC - Excluded"
    temp$sig[temp$ID %in% act_lists[[lib]][[cell]]$ID[act_lists[[lib]][[cell]]$exp_mean == 0]] <- "Failed QC - Excluded"
    
    failed_qc_snps <- temp %>%
      filter(sig == "Failed QC - Excluded") %>%
      pull(SNP)  # Extract unique SNPs
    
    temp <- temp %>%
      mutate(sig = ifelse(SNP %in% failed_qc_snps, "Failed QC - Excluded", sig))
    
    act_count_comp <- rbind(act_count_comp, temp)
  }
}

act_count_comp <- act_count_comp[act_count_comp$sig != "Failed QC - Excluded", ]

tmp_plotA<-ggplot(act_count_comp,aes(x=DNA_mean,y=log2FoldChange,color=sig)) +
  theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) +
  #scale_colour_manual(values=c("Not Significant"=rgb(0,0,0,200,maxColorValue=255),"Active"=rgb(55,126,184,200,maxColorValue=255),"Failed QC - Excluded"=rgb(230,149,57,200,maxColorValue=255))) +
  scale_colour_manual(values=c("Not Significant"=rgb(0,0,0,200,maxColorValue=255),"Active"=rgb(55,126,184,200,maxColorValue=255))) +
  geom_point_rast(alpha = .3,size=1, raster.dpi=1000) +
  scale_x_log10() +
  #coord_cartesian(xlim = c(10, 1000),ylim = c(-1.5,7.5)) +
  xlab("Normalized Tag Count - Plasmids") + ylab("MPRA Expression Fold Change\nlog_2(RNA/Plasmid)") +
  theme(legend.position = c(.75, .15),
        legend.key = element_blank(),
        legend.background = element_rect(color=rgb(0,0,0,150,maxColorValue=255), fill = "white", size = .5, linetype = "solid")) +
  guides(colour = guide_legend(override.aes = list(size=3,alpha=.7), title=NULL)) +
  geom_abline(intercept = 0, slope = 0,linetype = 1, size=.75, color=rgb(255,140,0,150,maxColorValue=255))

pdf(paste0(dir,"supplemental_1c_rna.pdf"), width=(100/25.4), height=(100/25.4))
print(tmp_plotA)
dev.off()

## Positive Controls 

pos_ctrls <- read.delim("positive_controls.out", stringsAsFactors = F, header = F)
pos_ctrls <- pos_ctrls$V1

act_count_comp_ctrl <- act_count_comp %>%
  filter(ID %in% pos_ctrls)

act_count_comp_ctrl <- act_count_comp_ctrl %>%
  mutate(lib_celltype = paste(celltype, library, sep = "_"))

# Function to add a 45-degree identity line
identity_line <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    coord_cartesian(xlim = c(0, 7), ylim = c(0, 7))  # Set fixed axis limits
}

# Function for rasterized points with subtle white gridlines in lower panels
raster_pointsB <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point_rast(size = 0.5, alpha = 0.3, raster.dpi = 1000) +
    coord_cartesian(xlim = c(0, 7), ylim = c(0, 7)) +  # Set fixed axis limits
    theme_minimal()  # Ensure a clean background
}

# Get unique cell types
celltype <- unique(act_count_comp_ctrl$celltype)

for (cell in celltype) {
  message(cell)
  data_filtered <- act_count_comp_ctrl %>% filter(celltype == cell)
  act_count_comp_ctrl_pivot <- data_filtered %>%
    select(ID, lib_celltype, log2FoldChange) %>%
    pivot_wider(names_from = lib_celltype, values_from = log2FoldChange)
  
  act_count_comp_ctrl_pivot$ID <- NULL
  
  # Filter data for the specific cell type
  
  # Generate ggpairs plot
  ggp_temp <- ggpairs(
    act_count_comp_ctrl_pivot, 
    diag = list(continuous = "blankDiag"),  # Removes diagonal density plot
    upper = list(continuous = wrap(cor_func, method = "pearson", size = 11)), 
    lower = list(combo = raster_pointsB, continuous = identity_line)
  ) + 
    theme_classic() +  # Ensure a clean theme
    theme(
      strip.text = element_text(size = 12),
      strip.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_line(size = 0.1, color = "grey80"),  # Subtle grid
      panel.grid.minor = element_line(size = 0.1, color = "grey80"),   # Lighter grid
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  # Determine the number of variables for aspect ratio
  num_vars <- ncol(act_count_comp_ctrl_pivot)  # Get number of columns
  square_size <- max(8, 2 * num_vars)  # Adjust scaling, ensuring a reasonable minimum size
  
  ggsave(filename = paste0(dir, "Positive_controls_correlation_", cell, ".png"), plot = ggp_temp, 
         width = square_size, height = square_size, dpi = 300)
  ggsave(filename = paste0(dir, "Positive_controls_correlation_", cell, ".pdf"), plot = ggp_temp, 
         width = square_size, height = square_size)
}


## Full Correlations
act_count_comp_ctrl_pivot <- act_count_comp_ctrl %>%
  select(ID, lib_celltype, log2FoldChange) %>%
  pivot_wider(names_from = lib_celltype, values_from = log2FoldChange)
act_count_comp_ctrl_pivot$ID <- NULL

cor_matrix <- cor(act_count_comp_ctrl_pivot, method = "pearson", use = "pairwise.complete.obs")

# Convert correlation matrix to a long format for ggplot
cor_data <- melt(cor_matrix)

# Convert Var1 and Var2 into factors with correct ordering
cor_data <- cor_data %>%
  mutate(col1_split = Var1, col2_split = Var2) %>%  # Duplicate original columns
  separate(col1_split, into = c("cell1", "library1"), sep = "_", extra = "merge") %>%
  separate(col2_split, into = c("cell2", "library2"), sep = "_", extra = "merge")

cor_data <- cor_data %>%
  arrange(cell1, cell2,library1,library2)

correct_order <- unique(c(cor_data$Var1, cor_data$Var2))

correct_order <- rev(correct_order)  # ⬅️ Flips the order to lower-left

cor_data <- cor_data %>%
  mutate(
    Var1 = factor(Var1, levels = correct_order, ordered = TRUE),  
    Var2 = factor(Var2, levels = correct_order, ordered = TRUE) 
  ) %>%
  arrange(Var1, Var2)

cust_breaks<-unique(cor_data$Var1)
cust_breaks_rev<-rev(unique(cor_data$Var1))

# Keep only the lower triangle (including excluding the diagonal)
cor_data <- cor_data %>%
  filter(as.numeric(Var1) >= as.numeric(Var2))  # Keep lower triangle only

# Create a correlation heatmap plot
# Generate clean correlation heatmap
ggp_temp<-ggplot(cor_data, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white", size = 0.25) +  # Thin gridlines for separation
  scale_fill_gradient(low = "blue", high = "orange", limits = c(0, 1), name = "Pearson") +
  theme_minimal(base_size = 14) +
  scale_x_discrete(limits=cust_breaks) + 
  scale_y_discrete(limits=rev(cust_breaks)) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotated axis labels for clarity
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),  # Remove axis titles for a clean look
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) 
  #labs(title = "Correlation of Positive Control Sequences")

# Count number of unique labels (to set aspect ratio)
num_vars <- length(unique(c(cor_data$Var1, cor_data$Var2)))
square_size <- 0.3 * num_vars  # Adjust multiplier for best fit

ggsave(filename = paste0(dir, "Positive_controls_correlation_All.png"), plot = ggp_temp, 
       width = square_size, height = square_size, dpi = 300)
ggsave(filename = paste0(dir, "Positive_controls_correlation_All.pdf"), plot = ggp_temp, 
       width = square_size, height = square_size)

########################################################
##Supplelental Figure 1 (Correaltion with colored squares)

count_files <- read.delim("norm_counts.files", header=F, stringsAsFactors = F)
correlation_list <- list()

pairwise_correlations <- function(data, lib) {
  # Compute correlation matrix
  corr_matrix <- cor(data, method="pearson", use='complete.obs')
  
  # Get column names
  cols <- colnames(data)
  
  # Generate all unique pairs
  results <- data.frame(Lib = character(), Column_1 = character(), Column_2 = character(), Correlation = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:(length(cols) - 1)) {
    for (j in (i + 1):length(cols)) {
      results <- rbind(results, data.frame(
        Lib = lib,
        Column_1 = cols[i],
        Column_2 = cols[j],
        Correlation = corr_matrix[i, j]
      ))
    }
  }
  return(results)
}


final_cor_results <- data.frame(Lib = character(), Column_1 = character(), Column_2 = character(), Correlation = numeric(), stringsAsFactors = FALSE)
final_corlog_results <- data.frame(Lib = character(), Column_1 = character(), Column_2 = character(), Correlation = numeric(), stringsAsFactors = FALSE)

for(i in 1:nrow(count_files)){
  count_plot <- read.delim(paste0(filedir,count_files[i,1]), stringsAsFactors = F, row.names=1)
  lib <- count_files[i,2]
  message(lib)
  count_plot <- count_plot[,!grepl("gm12878", colnames(count_plot), ignore.case = T)]
  cols_split <- colsplit(colnames(count_plot), "_", c("cell","rep"))
  cols_split$name <- paste(cols_split$cell, cols_split$rep, sep = "_")

  ggp_temp <- ggpairs(count_plot, diag=list("naDiag"), upper=list(continuous = wrap("cor", size=6)), lower=list(combo = raster_points))
  png(paste0(dir,lib,"_rep_norm_count.png"), height = 2048, width = 2048)
  print(ggp_temp + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)))#, strip.text.x=element_text(size=2), strip.text.y=element_text(size=2), axis.text=element_text(size=2),))
  dev.off()
 
  temp_result <- pairwise_correlations(count_plot, lib)
  final_cor_results <- rbind(final_cor_results, temp_result)  
  
  count_plot <- log2(count_plot+1)
  ggp_temp <- ggpairs(count_plot, diag=list("naDiag"), upper=list(continuous=wrap(cor_func, method="pearson", size=8)), lower=list(combo = raster_points), showStrips=F) + theme_minimal() + theme(strip.text = element_text(size = 12))
  png(paste0(dir,lib,"_rep_norm_count_log2.png"), height = 2048, width = 2048)
  print(ggp_temp + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), panel.grid.minor = element_blank(), 
                         panel.grid.major = element_blank()))#, strip.text.x=element_text(size=2), strip.text.y=element_text(size=2), axis.text=element_text(size=2),panel.spacing=grid::unit(0.1,"lines")))
  dev.off()
  
  temp_result <- pairwise_correlations(count_plot, lib)
  final_corlog_results <- rbind(final_corlog_results, temp_result)  
  
  count_plot <- count_plot[,which(colnames(count_plot) %in% cols_split$name[which(cols_split$rep %in% c("r1","r2","r3"))])]
  ggp_temp <- ggpairs(count_plot, diag=list("naDiag"), upper=list(continuous=wrap(cor_func, method="pearson", size=13)), lower=list(combo = raster_points), showStrips=F) + theme_minimal() + theme(strip.text = element_text(size = 12))
  png(paste0(dir,lib,"_rep_norm_count_log2_r123.png"), height = 2048, width = 2048)
  print(ggp_temp + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), panel.grid.minor = element_blank(), 
                         panel.grid.major = element_blank()))#, strip.text.x=element_text(size=2), strip.text.y=element_text(size=2), axis.text=element_text(size=2),panel.spacing=grid::unit(0.1,"lines")))
  dev.off()
  
  p1 <- ggcorr(count_plot, hjust=0.75,layout.exp=1.5,legend.position = "bottom", mid = "blue", high = "orange")
  p2 <- ggcorr(count_plot, legend.position = "none", size=0, mid = "blue", high = "orange")
  
  pdf(paste0(dir,lib,"_rep_norm_count_comps_colored_log2scores.pdf"), width=100/25.4, height = 100/25.4)
  print(p1)
  dev.off()
 
  pdf(paste0(dir,lib,"_rep_norm_count_comps_colored_log2scores.noLabels.pdf"), width=100/25.4, height = 100/25.4)
  print(p2)
  dev.off()
}

##Write Correlation Values

write.table(final_corlog_results, file = paste0(dir, "correlation_values_log.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(final_cor_results, file = paste0(dir, "correlation_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


##Box Plot of log2 transformed counts across all libraries

df_filtered <- final_corlog_results %>%
  mutate(
    Sample_1 = str_extract(Column_1, "^[^_]+"),  # Extracts everything before the first "_"
    Sample_2 = str_extract(Column_2, "^[^_]+")   # Extracts everything before the first "_"
  )

# Filter rows where Column_1 and Column_2 belong to the same sample but different replicates
df_filtered <- df_filtered %>%
  filter(Sample_1 == Sample_2)

correlation_summary <- df_filtered %>%
  group_by(Lib) %>%
  summarize(
    Min_Correlation = min(Correlation),
    Max_Correlation = max(Correlation),
    Mean_Correlation = mean(Correlation),
    Median_Correlation = median(Correlation),
    .groups = "drop"
  )

correlation_summary2 <- df_filtered %>%
  summarize(
    Min_Correlation = min(Correlation),
    Max_Correlation = max(Correlation),
    Mean_Correlation = mean(Correlation),
    Median_Correlation = median(Correlation),
    .groups = "drop"
  )

# Plot boxplots for each library

color_palette <- c(
  "plasmid" = "#000000",
  "a549" = "#EDC132",
  "hct116" = "#2E9093",
  "hepg2" = "#89689D",
  "k562" = "#086B8B",
  "sknsh" = "#DD4124"
)

# Reorder libraries so "OL41_42" is last
df_filtered$Lib <- factor(df_filtered$Lib, levels = c(setdiff(unique(df_filtered$Lib), "OL41_42"), "OL41_42"))
df_filtered$Sample_1 <- factor(df_filtered$Sample_1, levels = c("plasmid", "a549", "hct116", "hepg2", "k562", "sknsh"))

# Create the plot
p <- ggplot(df_filtered, aes(x = Lib, y = Correlation)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  # Black & white boxplot, no outliers
  geom_jitter(aes(color = Sample_1), width = 0.2, size = 2, alpha = 0.4) +  # Color points by Sample_1
  labs(title = "Pairwise Correlation of Replicates by Library",
       x = "Library",
       color = "Cell type",
       y = "Pearson correlation of log_2(counts)") +  # Change legend title to "Cell type"
  ylim(0,1) +  # Set y-axis limits from 0 to 1
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Keep only axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  scale_color_manual(values = color_palette)  # Apply custom colors

# Display the plot
pdf(paste0(dir,"All_libraries_replicate_cor.pdf"), width=8, height = 4.5)
print(p)
dev.off()

##Box Plot of untransformed counts across all libraries

df_lin_filtered <- final_cor_results %>%
  mutate(
    Sample_1 = str_extract(Column_1, "^[^_]+"),  # Extracts everything before the first "_"
    Sample_2 = str_extract(Column_2, "^[^_]+")   # Extracts everything before the first "_"
  )

# Filter rows where Column_1 and Column_2 belong to the same sample but different replicates
df_lin_filtered <- df_lin_filtered %>%
  filter(Sample_1 == Sample_2)

correlation_summary <- df_lin_filtered %>%
  group_by(Lib) %>%
  summarize(
    Min_Correlation = min(Correlation),
    Max_Correlation = max(Correlation),
    Mean_Correlation = mean(Correlation),
    Median_Correlation = median(Correlation),
    .groups = "drop"
  )

correlation_summary2 <- df_lin_filtered %>%
  summarize(
    Min_Correlation = min(Correlation),
    Max_Correlation = max(Correlation),
    Mean_Correlation = mean(Correlation),
    Median_Correlation = median(Correlation),
    .groups = "drop"
  )

# Reorder libraries so "OL41_42" is last
df_lin_filtered$Lib <- factor(df_lin_filtered$Lib, levels = c(setdiff(unique(df_lin_filtered$Lib), "OL41_42"), "OL41_42"))
df_lin_filtered$Sample_1 <- factor(df_lin_filtered$Sample_1, levels = c("plasmid", "a549", "hct116", "hepg2", "k562", "sknsh"))

# Create the plot
p <- ggplot(df_lin_filtered, aes(x = Lib, y = Correlation)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  # Black & white boxplot, no outliers
  geom_jitter(aes(color = Sample_1), width = 0.2, size = 2, alpha = 0.4) +  # Color points by Sample_1
  labs(title = "Pairwise Correlation of Replicates by Library (untransformed counts)",
       x = "Library",
       color = "Cell type",
       y = "Pearson correlation of counts")  +  # Change legend title to "Cell type"
  ylim(0,1) +  # Set y-axis limits from 0 to 1
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Keep only axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  scale_color_manual(values = color_palette)  # Apply custom colors

# Display the plot
pdf(paste0(dir,"All_libraries_replicate_corLIN.pdf"), width=8, height = 4.5)
print(p)
dev.off()


##Combined Plot
df_filtered$Transformation <- factor("log_2 (Sup Fig 1)", levels = c("Untransformed","log_2 (Sup Fig 1)"))
df_lin_filtered$Transformation <- factor("Untransformed", levels = c("Untransformed","log_2 (Sup Fig 1)"))

# Combine transformed and untransformed data
df_combined <- bind_rows(df_filtered, df_lin_filtered)

# Reorder libraries so "OL41_42" is last
df_combined$Lib <- factor(df_combined$Lib, levels = c(setdiff(unique(df_combined$Lib), "OL41_42"), "OL41_42"))

# Create the merged boxplot
p <- ggplot(df_combined, aes(x = Lib, y = Correlation, color = Transformation)) + 
  geom_boxplot(outlier.shape = NA, fill = NA, size = .7) +  # Transparent fill, outline color defined by aes()
  labs(title = "Pairwise Correlation of Replicates by Library",
       x = "Library",
       y = "Pearson Correlation",
       color = "Data Transformation") +  # Legend title
  ylim(0, 1) +  # Set y-axis limits from 0 to 1
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Keep only axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.position = "right"  # Move legend to the right
  ) +
  scale_color_manual(values = c("log_2 (Sup Fig 1)" = "red", "Untransformed" = "darkgrey"))

# Save the plot
pdf(paste0(dir, "All_libraries_replicate_cor_combined.pdf"), width = 8, height = 4.5)
print(p)
dev.off()