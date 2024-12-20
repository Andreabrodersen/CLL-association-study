#RF_177_with_CV

# Load necessary libraries
library(randomForest)
library(ggplot2)
library(caret)
library(dplyr)
library(progress)
library(tidyverse)
library(ggpubr)


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


Mut <- mutate(Mut, CK = rowSums(select(Mut, c(1:4)), na.rm = TRUE))
Mut <- Mut %>%
  mutate(CK = ifelse(CK >= 3, 1, 0))
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

set.seed(4321)

# Clean column names
colnames(rf_df_flow_single_na_omit) <- gsub("\\(|\\)|\\-|\\/|\\&|\\.|\\_", "", colnames(rf_df_flow_single_na_omit))
column_names <- colnames(rf_df_flow_single_na_omit)
column_names_adjusted <- gsub("^([0-9])", "a\\1", column_names)
colnames(rf_df_flow_single_na_omit) <- column_names_adjusted



# Create train-control object for 10-fold cross-validation
train_control <- trainControl(method = "cv", number = 10)

start_col <- 42
end_col <- ncol(rf_df_flow_single_na_omit)
df_ImpData_single <- data.frame()



#10-fold testing: 

for (col_index in start_col:end_col) {
  
  # Get the column name for the current col_index
  col_name <- colnames(rf_df_flow_single_na_omit)[col_index]
  
  # Create random forest model
  rf.fit <- randomForest(rf_df_flow_single_na_omit[, col_index] ~ ., 
                         data = rf_df_flow_single_na_omit[, -c(42:ncol(rf_df_flow_single_na_omit))], 
                         ntree = 500,
                         keep.forest = FALSE, 
                         importance = TRUE)
  
  
  # Extract variable importance data
  ImpData_nocv <- as.data.frame(importance(rf.fit))
  ImpData_nocv$Var.Names <- row.names(ImpData_nocv)
  ImpData_nocv$Var.Names <- factor(ImpData_nocv$Var.Names, levels = ImpData_nocv$Var.Names[order(ImpData_nocv$`%IncMSE`)])
  ImpData_nocv$Drug <- col_name
  df_ImpData_single_nocv <- rbind(df_ImpData_single,ImpData_nocv)
}

# Predicting flow from single DSSs
start_col <- 42
end_col <- ncol(rf_df_flow_single_na_omit)
df_ImpData_single <- data.frame()

# Initialize progress bar


# List to store all cross-validation results
cv_summary_list <- list()
cv_results_list <- list()
oob_results_list <-list()

for (col_index in start_col:end_col) {
  col_name <- colnames(rf_df_flow_single_na_omit)[col_index]
  
  
  tunegrid <- expand.grid(.mtry = 13)
  # Create the dataset for the current column
  response <- rf_df_flow_single_na_omit[, col_index]
  predictors <- rf_df_flow_single_na_omit[, 1:41]
  
  data_subset <- data.frame(response, predictors)
  
  # Create random forest model with cross-validation
  rf.fit <- train(response ~ ., 
                  data = data_subset, 
                  method = "rf", 
                  trControl = train_control,
                  tuneGrid=tunegrid,
                  importance = TRUE)
  
  
  # Cross-validation results
  cv_results <- rf.fit$resample
  cv_summary <- rf.fit$results
  oob_results<- data.frame(RMSE=sqrt(rf.fit$finalModel$mse),Rsquared=rf.fit$finalModel$rsq)
  
  # Store cv_summary in list
  cv_summary_list[[col_name]] <- cv_summary
  cv_results_list[[col_name]]<- cv_results
  oob_results_list[[col_name]]<- oob_results
  

}


# Combine all cv_summary into a single data frame
cv_summary_combined <- bind_rows(cv_summary_list, .id = "Response")
cv_results_combined <- bind_rows(cv_results_list, .id = "Response")
oob_results_combined <- bind_rows(oob_results_list, .id = "Response")
oob_summary_combined<- oob_results_combined %>% group_by(Response)%>%summarize_if(is.numeric,mean)

cv_results_combined$Response <- factor(cv_results_combined$Response, 
                                       levels=cv_summary_combined$Response[order(cv_summary_combined$RMSE)])
oob_results_combined$Response <- factor(oob_results_combined$Response, 
                                        levels=cv_summary_combined$Response[order(cv_summary_combined$RMSE)])





ggarrange(
  ggplot(bind_rows(cv_summary_combined %>% mutate(data="10-fold testing"),oob_summary_combined %>% mutate(data="OOB")),
         aes(x = data, y = RMSE)) +
    geom_line(aes(group = Response)) +
    geom_point(color="#93b9ed") +
    labs(
      x = "",
      y = "Mean RMSE",
      #color = "Response"
    ) +
    theme_minimal(),
  
  
  ggplot(bind_rows(cv_summary_combined %>% mutate(data="10-fold testing")),
         aes(x = data, y = Rsquared)) +
    geom_point(color="#93b9ed") +
    labs(
      x = "",
      y = "Mean R-squared",
      #color = "Response"
    ) +
    theme_minimal(), common.legend = TRUE
)


# optional plots for inspecting RMSE in test vs OOB in detail---- 


# ggplot(bind_rows(cv_results_combined %>% mutate(data="10-fold testing"),oob_results_combined %>% mutate(data="OOB")),
#        aes(x = Response, y = RMSE)) +
#   geom_boxplot(aes(fill=data), outlier.size=0.1) +
#   facet_grid(~data)+
#   labs(
#     title = "RMSE Scores",
#     x = "",
#     y = "Mean RMSE",
#     #color = "Response"
#   ) +
#   theme_minimal()+
#   theme(axis.text.x = element_blank())


# ggplot(bind_rows(cv_results_combined %>% mutate(data="10-fold testing"),oob_results_combined %>% mutate(data="OOB")),
#        aes(x = Response, y = Rsquared)) +
#   geom_boxplot(aes(fill=data), outlier.size=0.1) +
#   facet_grid(~data)+
#   labs(
#     title = "Mean R-squared",
#     x = "",
#     y = "Mean R-squared",
#     #color = "Response"
#   ) +
#   theme_minimal()+
#   theme(axis.text.x = element_blank())

# ggplot(bind_rows(cv_summary_combined %>% mutate(data="10-fold testing"),oob_summary_combined %>% mutate(data="OOB")),
#        aes(x = data, y = Rsquared, color = Response)) +
#   geom_line(aes(group = Response)) +
#   geom_point() +
#   labs(
#     x = "",
#     y = "Mean R-squared",
#     color = "Response"
#   ) +
#   theme_minimal()


# Summary plot of RMSE vs. Number of Trees
# summary_plot <- ggplot(cv_summary_combined, aes(x = mtry, y = RMSE, color = Response)) +
#   geom_line() +
#   geom_point() +
#   labs(
#     title = "Cross-Validation RMSE Scores",
#     x = "Number of Trees",
#     y = "RMSE",
#     color = "Response"
#   ) +
#   theme_minimal()
# 
# # Save summary plot
# ggsave("cv_summary_plot.png", summary_plot, width = 10, height = 6)
# 
# # Print summary plot
# print(summary_plot)


# Extract variable importance data -------
ImpData <- as.data.frame(varImp(rf.fit)$importance)
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- factor(ImpData$Var.Names, levels = ImpData$Var.Names[order(ImpData$Overall)])
ImpData$Drug <- col_name

df_ImpData_sum_single <- df_ImpData_single_nocv %>% group_by(Var.Names) %>% 
  summarise(Mean = mean(Overall),
            SD = sd(Overall), 
            Max = max(Overall),
            Which.max = Drug[which.max(Overall)])
df_ImpData_sum_single$tval <- df_ImpData_sum_single$Mean / df_ImpData_sum_single$SD

df_ImpData_sum_single$Var.Names <- factor(df_ImpData_sum_single$Var.Names, levels = df_ImpData_sum_single$Var.Names[order(df_ImpData_sum_single$Max)])
ggplot(df_ImpData_sum_single, aes(x = Var.Names, y = Max)) +
  geom_point(aes(), color = "#93b9ed", alpha = 0.7, size = 2) +
  geom_text(aes(y = Max + 3, label = Which.max), size = 3) +
  theme_light() +
  coord_flip() +
  labs(
    x = "Phosphorylation target",
    y = "Overall Importance",
    title = "Predictive power of phosphorylation targets on single treatments"
  ) +
  ylim(min(df_ImpData_sum_single$Max) - 2, max(df_ImpData_sum_single$Max) + 6) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )






# Predicting DSSs from flow COMBO TEST-----
train_control <- trainControl(method = "cv", number = 10)

start_col <- 42
end_col <- ncol(rf_df_flow_combo_na_omit)
df_ImpData_single <- data.frame()



#10-fold testing: 

for (col_index in start_col:end_col) {
  
  # Get the column name for the current col_index
  col_name <- colnames(rf_df_flow_combo_na_omit)[col_index]
  
  # Create random forest model
  rf.fit <- randomForest(rf_df_flow_combo_na_omit[, col_index] ~ ., 
                         data = rf_df_flow_combo_na_omit[, -c(42:ncol(rf_df_flow_combo_na_omit))], 
                         ntree = 500,
                         keep.forest = FALSE, 
                         importance = TRUE)
  
  
  # Extract variable importance data
  ImpData_nocv <- as.data.frame(importance(rf.fit))
  ImpData_nocv$Var.Names <- row.names(ImpData_nocv)
  ImpData_nocv$Var.Names <- factor(ImpData_nocv$Var.Names, levels = ImpData_nocv$Var.Names[order(ImpData_nocv$`%IncMSE`)])
  ImpData_nocv$Drug <- col_name
  df_ImpData_combo_nocv <- rbind(df_ImpData_combo,ImpData_nocv)
}

# Predicting flow from single DSSs
start_col <- 42
end_col <- ncol(rf_df_flow_combo_na_omit)
df_ImpData_combo <- data.frame()

# Initialize progress bar


# List to store all cross-validation results
cv_summary_list <- list()
cv_results_list <- list()
oob_results_list <-list()

for (col_index in start_col:end_col) {
  col_name <- colnames(rf_df_flow_single_na_omit)[col_index]
  
  
  tunegrid <- expand.grid(.mtry = 13)
  # Create the dataset for the current column
  response <- rf_df_flow_single_na_omit[, col_index]
  predictors <- rf_df_flow_single_na_omit[, 1:41]
  
  data_subset <- data.frame(response, predictors)
  
  # Create random forest model with cross-validation
  rf.fit <- train(response ~ ., 
                  data = data_subset, 
                  method = "rf", 
                  trControl = train_control,
                  tuneGrid=tunegrid,
                  importance = TRUE)
  
  
  # Cross-validation results
  cv_results <- rf.fit$resample
  cv_summary <- rf.fit$results
  oob_results<- data.frame(RMSE=sqrt(rf.fit$finalModel$mse),Rsquared=rf.fit$finalModel$rsq)
  
  # Store cv_summary in list
  cv_summary_list[[col_name]] <- cv_summary
  cv_results_list[[col_name]]<- cv_results
  oob_results_list[[col_name]]<- oob_results
  
  
}


# Combine all cv_summary into a single data frame
cv_summary_combined <- bind_rows(cv_summary_list, .id = "Response")
cv_results_combined <- bind_rows(cv_results_list, .id = "Response")
oob_results_combined <- bind_rows(oob_results_list, .id = "Response")
oob_summary_combined<- oob_results_combined %>% group_by(Response)%>%summarize_if(is.numeric,mean)

cv_results_combined$Response <- factor(cv_results_combined$Response, 
                                       levels=cv_summary_combined$Response[order(cv_summary_combined$RMSE)])
oob_results_combined$Response <- factor(oob_results_combined$Response, 
                                        levels=cv_summary_combined$Response[order(cv_summary_combined$RMSE)])





ggarrange(
  ggplot(bind_rows(cv_summary_combined %>% mutate(data="10-fold testing"),oob_summary_combined %>% mutate(data="OOB")),
         aes(x = data, y = RMSE)) +
    geom_line(aes(group = Response)) +
    geom_point(color="#93b9ed") +
    labs(
      x = "",
      y = "Mean RMSE",
      #color = "Response"
    ) +
    theme_minimal(),
  
  
  ggplot(bind_rows(cv_summary_combined %>% mutate(data="10-fold testing")),
         aes(x = data, y = Rsquared)) +
    geom_point(color="#93b9ed") +
    labs(
      x = "",
      y = "Mean R-squared",
      #color = "Response"
    ) +
    theme_minimal(), common.legend = TRUE
)






#Original------
colnames(rf_df_flow_combo_na_omit) = gsub("\\(|\\)|\\-|\\/|\\&|\\.", "", colnames(rf_df_flow_combo_na_omit))
column_names <- colnames(rf_df_flow_combo_na_omit)
column_names_adjusted <- gsub("^([0-9])", "a\\1", column_names)
colnames(rf_df_flow_combo_na_omit) <- column_names_adjusted

start_col <- 42
end_col <- ncol(rf_df_flow_combo_na_omit)
df_ImpData_combo <- data.frame()

# Create a folder to save the plots
dir.create("plots_combo", showWarnings = FALSE)

# Loop over the columns
for (col_index in start_col:end_col) {
  col_name <- colnames(rf_df_flow_combo_na_omit)[col_index]
  
  # Create random forest model with cross-validation
  rf.fit <- train(rf_df_flow_combo_na_omit[, col_index] ~ ., 
                  data = rf_df_flow_combo_na_omit[, -c(start_col:end_col)], 
                  method = "rf", 
                  trControl = train_control,
                  importance = TRUE)
  
  # Extract variable importance data
  ImpData <- as.data.frame(varImp(rf.fit)$importance)
  ImpData$Var.Names <- row.names(ImpData)
  ImpData$Var.Names <- factor(ImpData$Var.Names, levels = ImpData$Var.Names[order(ImpData$Overall)])
  ImpData$Drug <- col_name
  
  # # Create and save the plot
  # plot <- ggplot(ImpData, aes(x = Var.Names, y = Overall)) +
  #   geom_segment(aes(x = Var.Names, xend = Var.Names, y = 0, yend = Overall), color = "skyblue") +
  #   geom_point(aes(size = Overall), color = "blue", alpha = 0.6) +
  #   theme_light() +
  #   coord_flip() +
  #   labs(
  #     x = "Variable names",
  #     y = "Overall Importance",
  #     title = paste("Plot for", col_name)
  #   ) +
  #   theme(
  #     legend.position = "bottom",
  #     panel.grid.major.y = element_blank(),
  #     panel.border = element_blank(),
  #     axis.ticks.y = element_blank()
  #   )
  # 
  # # Save the plot
  # ggsave(paste0("plots_combo/plot_column_", col_name, ".png"), plot, width = 8, height = 6)
  
  df_ImpData_combo <- rbind(df_ImpData_combo, ImpData)
}

df_ImpData_sum_combo <- df_ImpData_combo %>% group_by(Var.Names) %>% 
  summarise(Mean = mean(Overall),
            SD = sd(Overall), 
            Max = max(Overall),
            Which.max = Drug[which.max(Overall)])
df_ImpData_sum_combo$tval <- df_ImpData_sum_combo$Mean / df_ImpData_sum_combo$SD

df_ImpData_sum_combo$Var.Names <- factor(df_ImpData_sum_combo$Var.Names, levels = df_ImpData_sum_combo$Var.Names[order(df_ImpData_sum_combo$Max)])
ggplot(df_ImpData_sum_combo, aes(x = Var.Names, y = Max)) +
  geom_point(aes(), color = "#93b9ed", alpha = 0.7, size = 2) +
  geom_text(aes(y = Max + 6, label = Which.max), size = 3) +
  theme_light() +
  coord_flip() +
  labs(
    x = "Phosphorylation target",
    y = "Overall Importance",
    title = "Predictive power of variables on combination treatments"
  ) +
  ylim(min(df_ImpData_sum_single$Max) - 2, max(df_ImpData_sum_single$Max) + 6) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )
