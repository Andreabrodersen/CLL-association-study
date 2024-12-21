#Radom forest OOB vs 10-fold testing
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

# 
# write.csv(rf_df_flow_combo_na_omit, "Rf_df_combo_nona_177.csv")
# write.csv(rf_df_flow_single_na_omit, "Rf_df_single_nona_177.csv")

set.seed(4321)

# Clean column names annd ass a before alll columnns begining with a number
colnames(rf_df_flow_single_na_omit) <- gsub("\\(|\\)|\\-|\\/|\\&|\\.|\\_", "", colnames(rf_df_flow_single_na_omit))
column_names <- colnames(rf_df_flow_single_na_omit)
column_names_adjusted <- gsub("^([0-9])", "a\\1", column_names)
colnames(rf_df_flow_single_na_omit) <- column_names_adjusted

colnames(rf_df_flow_combo_na_omit) <- gsub("\\(|\\)|\\-|\\/|\\&|\\.|\\_", "", colnames(rf_df_flow_combo_na_omit))
column_names <- colnames(rf_df_flow_combo_na_omit)
column_names_adjusted <- gsub("^([0-9])", "a\\1", column_names)
colnames(rf_df_flow_combo_na_omit) <- column_names_adjusted






#SINGLE: 
# Create train-control object for 10-fold cross-validation
train_control <- trainControl(method = "cv", number = 10)
start_col <- 42
end_col <- ncol(rf_df_flow_single_na_omit)
df_ImpData_single <- data.frame()
for (col_index in start_col:end_col) {
  
  # Get the column name for the current col_index
  col_name <- colnames(rf_df_flow_single_na_omit)[col_index]
  
  # Create random forest model
  rf.fit_no_cv <- randomForest(rf_df_flow_single_na_omit[, col_index] ~ ., 
                               data = rf_df_flow_single_na_omit[, -c(42:ncol(rf_df_flow_single_na_omit))], 
                               ntree = 500,
                               keep.forest = FALSE, 
                               importance = TRUE)
  
  
  # Extract variable importance data
  ImpData_nocv <- as.data.frame(importance(rf.fit_no_cv))
  ImpData_nocv$Var.Names <- row.names(ImpData_nocv)
  ImpData_nocv$Var.Names <- factor(ImpData_nocv$Var.Names, levels = ImpData_nocv$Var.Names[order(ImpData_nocv$`%IncMSE`)])
  ImpData_nocv$Drug <- col_name
  df_ImpData_single_nocv <- rbind(df_ImpData_single,ImpData_nocv)
}

# Predicting flow from single DSSs
start_col <- 42
end_col <- ncol(rf_df_flow_single_na_omit)
df_ImpData_single <- data.frame()



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




#combo----

# Create train-control object for 10-fold cross-validation
train_control <- trainControl(method = "cv", number = 10)
start_col <- 42
end_col <- ncol(rf_df_flow_combo_na_omit)
df_ImpData_combo <- data.frame()
for (col_index in start_col:end_col) {
  
  # Get the column name for the current col_index
  col_name <- colnames(rf_df_flow_combo_na_omit)[col_index]
  
  # Create random forest model
  rf.fit_no_cv <- randomForest(rf_df_flow_combo_na_omit[, col_index] ~ ., 
                               data = rf_df_flow_combo_na_omit[, -c(42:ncol(rf_df_flow_combo_na_omit))], 
                               ntree = 500,
                               keep.forest = FALSE, 
                               importance = TRUE)
  
  
  # Extract variable importance data
  ImpData_nocv <- as.data.frame(importance(rf.fit_no_cv))
  ImpData_nocv$Var.Names <- row.names(ImpData_nocv)
  ImpData_nocv$Var.Names <- factor(ImpData_nocv$Var.Names, levels = ImpData_nocv$Var.Names[order(ImpData_nocv$`%IncMSE`)])
  ImpData_nocv$Drug <- col_name
  df_ImpData_combo <- rbind(df_ImpData_combo,ImpData_nocv)
}

# Predicting flow from single DSSs
start_col <- 42
end_col <- ncol(rf_df_flow_combo_na_omit)
df_ImpData_single <- data.frame()



# List to store all cross-validation results
cv_summary_list <- list()
cv_results_list <- list()
oob_results_list <-list()

for (col_index in start_col:end_col) {
  col_name <- colnames(rf_df_flow_combo_na_omit)[col_index]
  
  
  tunegrid <- expand.grid(.mtry = 13)
  # Create the dataset for the current column
  response <- rf_df_flow_combo_na_omit[, col_index]
  predictors <- rf_df_flow_combo_na_omit[, 1:41]
  
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
