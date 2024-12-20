
#library-----
library(caret)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(randomForest)
library(shapper)
library(shapley)
library(DALEX)

#files----
#mutation file: 
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
Mut <- mutate(Mut, CK = rowSums(select(Mut, c(1:4)), na.rm = TRUE))
Mut <- Mut %>%
  mutate(CK = ifelse(CK >= 3, 1, 0))
Mut$Survival_state<-NULL
Mut$PFS_category_median<-NULL
Mut$`GC (Number of CNA)`<-NULL
Mut$CK<-NULL #only seven cases

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




#DSS file: 
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

Combodrugs<-colnames(DSS_sub_combo)
Singledrugs<-colnames(DSS_sub_single)


#make one big df sinbgle: 
rf_df_flow_single<-bind_cols(Mut,flow_ctrl,DSS_sub_single)
colnames_bigdf<-colnames(rf_df_flow_single)
rf_df_flow_single<-data.frame(rf_df_flow_single)
colnames(rf_df_flow_single)<-colnames_bigdf

colnames(rf_df_flow_single)=gsub("\\(|\\)|\\-|\\/|\\&|\\.","", colnames(rf_df_flow_single))
colnames(rf_df_flow_single)=gsub(" ","", colnames(rf_df_flow_single))
rf_df_flow_single_na_omit<-na.omit(rf_df_flow_single)



#make one big df coombo: 
rf_df_flow_combo<-bind_cols(Mut,flow_ctrl,DSS_sub_combo)
colnames_bigdf_combo<-colnames(rf_df_flow_combo)
rf_df_flow_combo<-data.frame(rf_df_flow_combo)
colnames(rf_df_flow_combo)<-colnames_bigdf_combo

colnames(rf_df_flow_combo)=gsub("\\(|\\)|\\-|\\/|\\&|\\.","", colnames(rf_df_flow_combo))
colnames(rf_df_flow_combo)=gsub(" ","", colnames(rf_df_flow_combo))
rf_df_flow_combo_na_omit<-na.omit(rf_df_flow_combo)


write.csv(rf_df_flow_combo_na_omit, "Rf_df_combo_nona_177.csv")
write.csv(rf_df_flow_single_na_omit, "Rf_df_single_nona_177.csv")


#randomforeset
#####################################


set.seed(111)

colnames(rf_df_flow_single_na_omit)=gsub("\\(|\\)|\\-|\\/|\\&|\\.|\\_","", colnames(rf_df_flow_single_na_omit))
column_names <- colnames(rf_df_flow_single_na_omit)
column_names_adjusted <- gsub("^([0-9])", "a\\1", column_names)

colnames(rf_df_flow_single_na_omit) <- column_names_adjusted

#predicting flow from single DSSs 
#looping over each column: ----
# Specify the range of columns to loop over
start_col <- 42
end_col <- ncol(rf_df_flow_single_na_omit)
df_ImpData_single<-c()

# Create a folder to save the plots
dir.create("volcanos_rf_single", showWarnings = FALSE)

# Loop over the columns
for (col_index in start_col:end_col) {
  
  # Get the column name for the current col_index
  col_name <- colnames(rf_df_flow_single_na_omit)[col_index]
  
  # Create random forest model
  rf.fit <- randomForest(rf_df_flow_single_na_omit[, col_index] ~ ., 
                         data = rf_df_flow_single_na_omit[, -c( 42:ncol(rf_df_flow_single_na_omit))], 
                         ntree = 500,
                         keep.forest = FALSE, 
                         importance = TRUE)
  
  # Extract variable importance data
  ImpData <- as.data.frame(importance(rf.fit))
  ImpData$Var.Names <- row.names(ImpData)
  ImpData$Var.Names <- factor(ImpData$Var.Names, levels = ImpData$Var.Names[order(ImpData$`%IncMSE`)])
  ImpData$Drug <- col_name
  # Create and save the plot
  plot <- ggplot(ImpData, aes(x = Var.Names, y = `%IncMSE`)) +
    geom_segment(aes(x = Var.Names, xend = Var.Names, y = 0, yend = `%IncMSE`), color = "skyblue") +
    geom_point(aes(size = IncNodePurity), color = "blue", alpha = 0.6) +
    theme_light() +
    coord_flip() +
    labs(
      x = "Variable names",
      y = "%IncMSE",
      title = paste("Plot for", col_name)
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  # Save the plot
  ggsave(paste0("volcanos_rf_single/plot_column_", col_name, ".png"), plot, width = 8, height = 6)
  
  df_ImpData_single <- rbind(df_ImpData_single,ImpData)
}


df_ImpData_sum_single <- df_ImpData_single %>% group_by(Var.Names) %>% 
  summarise(Mean=mean(`%IncMSE`),
            SD=sd(`%IncMSE`), 
            Max=max(`%IncMSE`),
            Which.max=Drug[which.max(`%IncMSE`)])
df_ImpData_sum_single$tval <- df_ImpData_sum_single$Mean/df_ImpData_sum_single$SD



df_ImpData_sum_single$Var.Names <- factor(df_ImpData_sum_single$Var.Names, levels=df_ImpData_sum_single$Var.Names[order(df_ImpData_sum_single$Max)])
ggplot(df_ImpData_sum_single, aes(x = Var.Names, y = Max)) +
  geom_point(aes(), color = "#93b9ed", alpha = 0.7, size=2) +
  geom_text(aes(y=Max+3,label=Which.max), size = 3) +
  theme_light() +
  coord_flip() +
  labs(
    x = "Phosphorylation target",
    y = "%IncMSE",
    title = "Predictive power of phosphorylation
    targets on single treatments"
  ) +
  ylim(min(df_ImpData_sum_single$Max)-3, max(df_ImpData_sum_single$Max)+6)+
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )




#predicting DSSs from flow COMBO------
set.seed(4321)

colnames(rf_df_flow_combo_na_omit)=gsub("\\(|\\)|\\-|\\/|\\&|\\.","", colnames(rf_df_flow_combo_na_omit))
column_names <- colnames(rf_df_flow_combo_na_omit)
column_names_adjusted <- gsub("^([0-9])", "a\\1", column_names)
colnames(rf_df_flow_combo_na_omit)<-column_names_adjusted
df_ImpData_combo<-c()


start_col <- 42
end_col <- ncol(rf_df_flow_combo)

# Create a folder to save the plots
dir.create("plots_combo", showWarnings = FALSE)

# Loop over the columns
df_ImpData_combo <- c()
for (col_index in start_col:end_col) {
  
  # Get the column name for the current col_index
  col_name <- colnames(rf_df_flow_combo_na_omit)[col_index]
  
  # Create random forest model
  rf.fit <- randomForest(rf_df_flow_combo_na_omit[, col_index] ~ ., 
                         data = rf_df_flow_combo_na_omit[, -c( 42:ncol(rf_df_flow_combo_na_omit))], 
                         ntree = 500,
                         keep.forest = FALSE, 
                         importance = TRUE)
  
  # Extract variable importance data
  ImpData <- as.data.frame(importance(rf.fit))
  ImpData$Var.Names <- row.names(ImpData)
  ImpData$Var.Names <- factor(ImpData$Var.Names, levels = ImpData$Var.Names[order(ImpData$`%IncMSE`)])
  ImpData$Drug <- col_name
  # Create and save the plot
  plot <- ggplot(ImpData, aes(x = Var.Names, y = `%IncMSE`)) +
    geom_segment(aes(x = Var.Names, xend = Var.Names, y = 0, yend = `%IncMSE`), color = "skyblue") +
    geom_point(aes(size = IncNodePurity), color = "blue", alpha = 0.6) +
    theme_light() +
    coord_flip() +
    labs(
      x = "Variable names",
      y = "%IncMSE",
      title = paste("Plot for", col_name)
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  # Save the plot
  ggsave(paste0("plots_combo/plot_column_", col_name, ".png"), plot, width = 8, height = 6)
  
  df_ImpData_combo <- rbind(df_ImpData_combo,ImpData)
}


df_ImpData_sum_combo <- df_ImpData_combo %>% group_by(Var.Names) %>% 
  summarise(Mean=mean(`%IncMSE`),
            SD=sd(`%IncMSE`), 
            Max=max(`%IncMSE`),
            Which.max=Drug[which.max(`%IncMSE`)])
df_ImpData_sum_combo$tval <- df_ImpData_sum_combo$Mean/df_ImpData_sum_combo$SD


df_ImpData_sum_combo$Var.Names <- factor(df_ImpData_sum_combo$Var.Names, levels=df_ImpData_sum_combo$Var.Names[order(df_ImpData_sum_combo$Max)])
ggplot(df_ImpData_sum_combo, aes(x = Var.Names, y = Max)) +
  geom_point(aes(), color = "#93b9ed", alpha = 0.7, size=2) +
  geom_text(aes(y=Max+6,label=Which.max), size = 3) +
  theme_light() +
  coord_flip() +
  labs(
    x = "Phosphorylation target",
    y = "%IncMSE",
    title = "Predictive power of variables on combination treatments"
  ) +
  ylim(min(df_ImpData_sum_single$Max)-4, max(df_ImpData_sum_single$Max)+6)+
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )




