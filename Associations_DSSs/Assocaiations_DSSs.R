#Assocaiations DSSs

#library
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(ComplexUpset)
library(units)
library(openxlsx2)
library(dplyr)


#files----
DSSs<-read_csv("DSS_177_standardized.csv")
colnames_DSS<-colnames(DSSs)
DSSs<-data.frame(DSSs)
rownames_DSS<-as.character(DSSs$...1)
DSSs$...1<-NULL
rownames(DSSs)<-rownames_DSS
colnames_DSS<-colnames_DSS[-1]
colnames(DSSs)<-colnames_DSS
DSS_sub_combo <- DSSs[grep("&", colnames(DSSs))]
DSS_sub_single <- DSSs[-grep("&", colnames(DSSs))]


mutations_CLL <- read_csv("Mutations_177.csv")
Mut<-mutations_CLL
rm(mutations_CLL)

#Colnames
Combodrugs<-colnames(DSS_sub_combo)
Singledrugs<-colnames(DSS_sub_single)


#formatting mut:
colnames_Mut<-colnames(Mut)
Mut<-data.frame(Mut)
rownames_Mut<-as.character(Mut$...1)
Mut$...1<-NULL
rownames(Mut)<-rownames_Mut
colnames_Mut<-colnames_Mut[-1]
colnames(Mut)<-colnames_Mut
Mut_sub<-Mut
Mut_sub$pid<-NULL
Mut_sub <- mutate(Mut_sub, CK = rowSums(select(Mut_sub, c(1:5, 7)), na.rm = TRUE))
Mut_sub <- Mut_sub %>%
  mutate(CK = ifelse(CK >= 3, 1, 0))


#Analysis single----
p.mat.single.z <- matrix(NA, ncol = ncol(Mut_sub), nrow = ncol(DSS_sub_single))

for (i in 1:ncol(Mut_sub)) {
  for (j in 1:ncol(DSS_sub_single)) {
    drug_response <- DSS_sub_single[, j]
    gene_expression <- Mut_sub[, i]
    
    # Combine gene_expression and drug_response, excluding rows with NAs in either of them
    combined_data <- data.frame(gene_expression, drug_response)
    combined_data <- na.omit(combined_data)
    
    if (nrow(combined_data) < 2) {
      p.mat.single.z[j, i] <- NA  # Set the value to NA if there are not enough data points
    } else {
      # Debugging print statements
      print(paste("Gene:", i, "Drug:", j))
      print(table(combined_data$gene_expression))
      print(table(combined_data$drug_response))
      
      wilcox_result <- wilcox.test(combined_data$drug_response[combined_data$gene_expression == 0], 
                                   combined_data$drug_response[combined_data$gene_expression == 1])
      
      p.mat.single.z[j, i] <- wilcox_result$p.value
    }
  }
}

drug_names <- colnames(DSS_sub_single)
mutation_names <- colnames(Mut_sub)

# Set column names (drugs)
colnames(p.mat.single.z) <- mutation_names

# Set row names (mutations)
rownames(p.mat.single.z) <- drug_names



#directionality
# Assuming Singledrugs is the correct set of drug names
directionality_single_z <- matrix(NA, ncol = ncol(Mut_sub), nrow = length(Singledrugs))

# Loop through each gene in 'genes'
for (i in 1:ncol(Mut_sub)) {
  # Loop through each drug in 'drugs'
  for (j in 1:length(Singledrugs)) {
    p_value <- p.mat.single.z[j, i]
    
    if (!is.na(p_value) && p_value <= 0.05) {
      # Calculate the median of drug_response for two groups
      median_high_gene <- median(DSS_sub_single[gene_expression > median(gene_expression), Singledrugs[j]], na.rm = TRUE)
      median_low_gene <- median(DSS_sub_single[gene_expression <= median(gene_expression), Singledrugs[j]], na.rm = TRUE)
      
      # Determine directionality based on medians
      if (!is.na(median_high_gene) && !is.na(median_low_gene)) {
        if (median_high_gene > median_low_gene) {
          directionality_single_z[j, i] <- "Positive"
        } else if (median_high_gene < median_low_gene) {
          directionality_single_z[j, i] <- "Negative"
        } else {
          directionality_single_z[j, i] <- "No Association"
        }
      } else {
        directionality_single_z[j, i] <- "Missing Median Values"
      }
    } else {
      directionality_single_z[j, i] <- "Not Significant"
    }
  }
}

# Set column names (genes)
colnames(directionality_single_z) <- mutation_names

# Set row names (drugs)
rownames(directionality_single_z) <- Singledrugs


#Adjusting P-values-----
adjusted_p.mat_single_z <- p.mat.single.z
for (i in 1:ncol(p.mat.single.z)) {
  adjusted_p.mat_single_z[, i] <- p.adjust(p.mat.single.z[, i], method = "fdr")
}





#Analysis combo----
p.mat.combo.z <- matrix(NA, ncol = ncol(Mut_sub), nrow = ncol(DSS_sub_combo))

for (i in 1:ncol(Mut_sub)) {
  for (j in 1:ncol(DSS_sub_combo)) {
    drug_response <- DSS_sub_combo[, j]
    gene_expression <- Mut_sub[, i]
    
    # Combine gene_expression and drug_response, excluding rows with NAs in either of them
    combined_data <- data.frame(gene_expression, drug_response)
    combined_data <- na.omit(combined_data)
    
    if (nrow(combined_data) < 2) {
      p.mat.combo.z[j, i] <- NA  # Set the value to NA if there are not enough data points
    } else {
      # Debugging print statements
      print(paste("Gene:", i, "Drug:", j))
      print(table(combined_data$gene_expression))
      print(table(combined_data$drug_response))
      
      wilcox_result <- wilcox.test(combined_data$drug_response[combined_data$gene_expression == 0], 
                                   combined_data$drug_response[combined_data$gene_expression == 1])
      
      p.mat.combo.z[j, i] <- wilcox_result$p.value
    }
  }
}

drug_names <- colnames(DSS_sub_combo)
mutation_names <- colnames(Mut_sub)

# Set column names (drugs)
colnames(p.mat.combo.z) <- mutation_names

# Set row names (mutations)
rownames(p.mat.combo.z) <- drug_names

#Adjusting P-values-----
adjusted_p.mat_combo_z <- p.mat.combo.z
for (i in 1:ncol(p.mat.combo.z)) {
  adjusted_p.mat_combo_z[, i] <- p.adjust(p.mat.combo.z[, i], method = "fdr")
}


write_xlsx(p.mat.combo.z, "Significant_drugs_table_paper3.xlsx",col_names=TRUE,row_names=TRUE)
write_xlsx(p.mat.single.z, "Significant_drugs_paper3.xlsx",col_names=TRUE,row_names=TRUE)
write.csv(data.frame(p.mat.combo.z),"sign_combo_z")
write.csv(data.frame(p.mat.single.z),"sign_single_z")

#volcanos----
dir.create("Volcano_plots")
#singel_volcanos----
# Loop through each mutation and create a volcano plot
for (mutation_name in colnames(Mut_sub)) {
  # Initialize empty lists to store data for each drug
  drug_data <- list()
  
  # Iterate through each drug to compare drug responses for patients with and without the mutation
  for (drug in colnames(DSS_sub_single)) {
    patients_with_mutation <- DSS_sub_single[which(Mut_sub[, mutation_name] == 1), drug]
    patients_without_mutation <- DSS_sub_single[which(Mut_sub[, mutation_name] == 0), drug]
    
    # Perform a statistical test (e.g., Wilcoxon test) and calculate p-value
    wilcox_result <- wilcox.test(patients_with_mutation, patients_without_mutation)
    p_value <- wilcox_result$p.value
    
    # Store drug name and p-value in the list
    drug_data[[drug]] <- data.frame(
      Drug = drug,
      PValue = p_value,
      ChangeInDSS = mean(patients_with_mutation) - mean(patients_without_mutation)
    )
  }
  
  # Combine the drug data into a single data frame
  drug_data_df <- do.call(rbind, drug_data)
  
  # Create the volcano plot using ggplot2
  volcano_plot <- ggplot(drug_data_df, aes(x = ChangeInDSS, y = -log10(PValue), color = factor(PValue < 0.05))) +
    geom_point(size = 3, shape = ifelse(drug_data_df$PValue < 0.05, 17, 16))  +
    scale_color_manual(values = c("TRUE" = "#93b9ed", "FALSE" = "black")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      x = "Change in Z-score (DSS)",
      y = "-log10(p)",
      title = paste("Volcano Plot for", mutation_name)
    ) +
    theme_minimal()+
    scale_x_continuous(limits = c(-0.5, 0.75))  # Set X-axis limits
  
  # Add labels for significant drugs
  volcano_plot <- volcano_plot +
    geom_text(
      data = subset(drug_data_df, PValue < 0.05),
      aes(label = Drug),
      vjust = +0.5,
      hjust = 0.75
    )+
    geom_text_repel(data = subset(drug_data_df, Drug=="Ibrutinib&Venetoclax"),
                    aes(label = Drug),
                    col="darkgray",
                    min.segment.length = 5,
                    vjust = +0.5,
                    hjust = 0.75
    )
  
  # Save the plot to a file in the "volcano_plots" directory
  plot_filename <- file.path("/Users/AmoBro/Desktop/PAPER 3 FINAL/FILE_EXTRACTION/Volcano_plots", paste0(mutation_name, "Single_z_volcano_plot.pdf"))
  ggsave(filename = plot_filename, plot = volcano_plot, width = 8, height = 6)
}


#combo vocanos----

# Loop through each mutation and create a volcano plot
for (mutation_name in colnames(Mut_sub)) {
  # Initialize empty lists to store data for each drug
  drug_data <- list()
  
  # Iterate through each drug to compare drug responses for patients with and without the mutation
  for (drug in colnames(DSS_sub_combo)) {
    patients_with_mutation <- DSS_sub_combo[which(Mut_sub[, mutation_name] == 1), drug]
    patients_without_mutation <- DSS_sub_combo[which(Mut_sub[, mutation_name] == 0), drug]
    
    # Perform a statistical test (e.g., Wilcoxon test) and calculate p-value
    wilcox_result <- wilcox.test(patients_with_mutation, patients_without_mutation)
    p_value <- wilcox_result$p.value
    
    # Store drug name and p-value in the list
    drug_data[[drug]] <- data.frame(
      Drug = drug,
      PValue = p_value,
      ChangeInDSS = mean(patients_with_mutation) - mean(patients_without_mutation)
    )
  }
  
  # Combine the drug data into a single data frame
  drug_data_df <- do.call(rbind, drug_data)
  
  # Create the volcano plot using ggplot2
  volcano_plot <- ggplot(drug_data_df, aes(x = ChangeInDSS, y = -log10(PValue), color = factor(PValue < 0.05))) +
    geom_point(size = 3, shape = ifelse(drug_data_df$PValue < 0.05, 17, 16)) +
    scale_color_manual(values = c("TRUE" = "#93b9ed", "FALSE" = "black")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      x = "Change in Z-score (DSS)",
      y = "-log10(p)",
      title = paste("Volcano Plot for", mutation_name)
    ) +
    theme_minimal()+
    scale_x_continuous(limits = c(-0.5, 0.5))  # Set X-axis limits
  
  # Add labels for significant drugs
  volcano_plot <- volcano_plot +
    geom_text_repel(
      data = subset(drug_data_df, PValue < 0.05),
      aes(label = Drug),
      min.segment.length = 0,
      vjust = +0.5,
      hjust = 0.75,
      max.overlaps = 70
    )+
    geom_text_repel(data = subset(drug_data_df, Drug=="Ibrutinib&Venetoclax"),
                    aes(label = Drug),
                    col="darkgray",
                    min.segment.length = 5,
                    vjust = +0.5,
                    hjust = 0.75,
                    max.overlaps = 70
    )
  
  
  
  # Save the plot to a file in the "volcano_plots" directory
  plot_filename <- file.path("/Users/AmoBro/Desktop/PAPER 3 FINAL/FILE_EXTRACTION/Volcano_plots", paste0(mutation_name, "Combo_z_volcano_plot.pdf"))
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
  X11q.deletion = Mut_sub$`11q deletion`, 
  trisomy.12 = Mut_sub$`trisomy 12`, 
  X13q14.deletion = Mut_sub$`13q14 deletion`,
  X17p13.deletion..TP53. = Mut_sub$`17p13 deletion (TP53)`,
  TP53.mutation = Mut_sub$`TP53-mutation`, 
  IGHV = Mut_sub$IGHV, 
  #IGVH_TP53 = Mut_sub$IGVH_TP53, 
  
  col = col,
  annotation_height = new_height 
)


#single reordering---- 
set.seed(321)
kclus <- kmeans(t(DSS_sub_single), 3)
kclus$cluster

split_default <- paste0("Cluster\n", kclus$cluster)
default.hmap.single <- Heatmap(as.matrix(t(DSS_sub_single)), split=split_default, cluster_row_slices = FALSE,
                               name = "Combo", col = c('#1111ed','#ffd866','#ed1111'),
                               row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize =6), 
                               top_annotation = ha,
                               column_title = "Patients",
                               row_title = "Drugs", column_names_rot = 90)

split_reorder <- factor(paste0("Cluster\n", kclus$cluster), levels=c("Cluster\n2","Cluster\n1","Cluster\n3"))
reorder.hmap.single <- Heatmap(as.matrix(t(DSS_sub_single)), split=split_reorder, cluster_row_slices = FALSE,
                               name = "Combo", col = c('#1111ed','#ffd866','#ed1111'),
                               row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize =6), 
                               top_annotation = ha,
                               column_title = "Patients",
                               row_title = "Drugs", column_names_rot = 90)



#combo reordered----
set.seed(4321)
kclus <- kmeans(t(DSS_sub_combo), 3)
kclus$cluster

split_default <- paste0("Cluster\n", kclus$cluster)
default.hmap.combo <- Heatmap(as.matrix(t(DSS_sub_combo)), split=split_default, cluster_row_slices = FALSE,
                              name = "Combo", col = c('#1111ed','#ffd866','#ed1111'),
                              row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize =6), 
                              top_annotation = ha,
                              column_title = "Patients",
                              row_title = "Drugs", column_names_rot = 90)

split_reorder <- factor(paste0("Cluster\n", kclus$cluster), levels=c("Cluster\n2","Cluster\n3","Cluster\n1"))
reorder.hmap.combo <- Heatmap(as.matrix(t(DSS_sub_combo)), split=split_reorder, cluster_row_slices = FALSE,
                              name = "Combo", col = c('#1111ed','#ffd866','#ed1111'),
                              row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize =6), 
                              top_annotation = ha,
                              column_title = "Patients",
                              row_title = "Drugs", column_names_rot = 90)

