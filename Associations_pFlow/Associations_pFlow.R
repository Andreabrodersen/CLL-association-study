#Associations_pFlow

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(ComplexHeatmap)
library(pheatmap)
library(ComplexUpset)
library(units)
#files 
mutations_CLL <- read_csv("Mutations_177.csv")
Mut<-mutations_CLL
rm(mutations_CLL)

colnames_Mut<-colnames(Mut)
Mut<-data.frame(Mut)
rownames_Mut<-as.character(Mut$...1)
Mut$...1<-NULL
rownames(Mut)<-rownames_Mut
colnames_Mut<-colnames_Mut[-1]
colnames(Mut)<-colnames_Mut

Mut$pid<-NULL
Mut <- mutate(Mut, CK = rowSums(select(Mut, c(1:5, 7)), na.rm = TRUE))
Mut <- Mut %>%
  mutate(CK = ifelse(CK >= 3, 1, 0))
Mut$Survival_state<-NULL



#reading in flow ctrl: 
flow_ctrl<-read_csv("Pflow_177_nonoise_log_scaled.csv")
colnames_flow<-colnames(flow_ctrl)
flow_ctrl<-data.frame(flow_ctrl)
rownames_flow<-as.character(flow_ctrl$...1)
flow_ctrl$...1<-NULL
rownames(flow_ctrl)<-rownames_flow
colnames_flow<-colnames_flow[-1]
colnames(flow_ctrl)<-colnames_flow

common_ids <- intersect(rownames(Mut), rownames(flow_ctrl))



# WR tests------
prot_name <- colnames(flow_ctrl)
mutation_names <- colnames(Mut)

# Loop through each mutation and create a volcano plot
for (mutation_name in colnames(Mut)) {
  # Initialize empty lists to store data for each drug
  drug_data <- list()
  
  # Iterate through each drug to compare drug responses for patients with and without the mutation
  for (prot_name in colnames(flow_ctrl)) {
    patients_with_mutation <- as.matrix(flow_ctrl[which(Mut[, mutation_name] == 1), prot_name])
    patients_without_mutation <- as.matrix(flow_ctrl[which(Mut[, mutation_name] == 0), prot_name])
    
    # Perform a statistical test (e.g., Wilcoxon test) and calculate p-value
    wilcox_result <- wilcox.test(patients_with_mutation, patients_without_mutation)
    p_value <- wilcox_result$p.value
    
    # Store drug name and p-value in the list
    drug_data[[prot_name]] <- data.frame(
      Drug = prot_name,
      PValue = p_value,
      ChangeInProt = mean(patients_with_mutation) - mean(patients_without_mutation)
    )
  }
  
  # Combine the drug data into a single data frame
  drug_data_df <- do.call(rbind, drug_data)
  
  # Create the volcano plot using ggplot2
  volcano_plot <- ggplot(drug_data_df, aes(x = ChangeInProt, y = -log10(PValue), color = factor(PValue < 0.05))) +
    geom_point(size = 3) +
    scale_color_manual(values = c("TRUE" = "#93b9ed", "FALSE" = "black")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      x = "Change in Z-score(log10(MFI))",
      y = "-log10(p)",
      title = paste("Volcano Plot for", mutation_name)
    ) +
    theme_minimal()+
    coord_cartesian(xlim = c(-0.4, 0.4))  # setting the x-axis limits
  
  # Add labels for significant drugs
  volcano_plot <- volcano_plot +
    geom_text(
      data = subset(drug_data_df, PValue < 0.05),
      aes(label = Drug),
      vjust = +0.5,
      hjust = 0.75
    )
  
  # Save the plot to a file in directory
  plot_filename <- file.path("/Users/AmoBro/Desktop/PAPER 3 FINAL/Figures_remake_240722", paste0(mutation_name, "_pflow_volcano_plot20240722.pdf"))
  ggsave(filename = plot_filename, plot = volcano_plot, width = 8, height = 6)
}


#heatmaps----
# Create a heatmap using pheatmap


# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
col = list(TP53.mutation = c("0" = "#7199A1", "1" = "#A35B4F", "NA"= "gray"),
           X17p13.deletion..TP53. = c("0" = "#7199A1", "1" = "#A35B4F", "NA"= "gray"), 
           X11q.deletion = c("0" = "#7199A1", "1" = "#A35B4F", "NA"= "gray"),
           trisomy.12 = c("0" = "#7199A1", "1" = "#A35B4F", "NA"= "gray"),
           IGHV = c("0" = "#7199A1", "1" = "#A35B4F", "NA"= "gray"),
           X13q14.deletion = c("0" = "#7199A1", "1" = "#A35B4F", "NA"= "gray")
)


#can be added if igvh_tp53 is of interest
#IGVH_TP53 = c("0" = "#4f97a3", "1" = "#A35B4F", "NA"= "gray"),


# Create the heatmap annotation

new_height <- unit(2.5, "mm")

ha <- HeatmapAnnotation(
  X11q.deletion = Mut$`11q deletion`, 
  trisomy.12 = Mut$`trisomy 12`, 
  X13q14.deletion = Mut$`13q14 deletion`,
  X17p13.deletion..TP53. = Mut$`17p13 deletion (TP53)`,
  TP53.mutation = Mut$`TP53-mutation`, 
  IGHV = Mut$IGHV, 
  #IGVH_TP53 = Mut_sub$IGVH_TP53, 
  
  col = col,
  annotation_height = new_height 
)

Heatmap(as.matrix(t(flow_ctrl)), cluster_row_slices = FALSE,
        name = "flow", col = c('#1111ed','#ffd866','#ed1111'),
        row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize =6), 
        top_annotation = ha,
        column_title = "Patients",
        row_title = "Drugs", column_names_rot = 90)

