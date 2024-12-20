#Random forest permutation tests


#library(caret)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(randomForest)
library(shapper)
library(DALEX)

set.seed(4321)

rf_df_flow_combo_na_omit<-read.csv("Rf_df_combo_nona_177.csv")
rf_df_flow_single_na_omit<-read.csv("Rf_df_single_nona_177.csv")

rownames(rf_df_flow_combo_na_omit)<-rf_df_flow_combo_na_omit$X
rf_df_flow_combo_na_omit$X<-NULL

rownames(rf_df_flow_single_na_omit)<-rf_df_flow_single_na_omit$X
rf_df_flow_single_na_omit$X<-NULL

#single----
start_col <- 42
end_col <- ncol(rf_df_flow_single_na_omit)
df_ImpData_single<-c()
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
  ImpData <- as.data.frame(importance(rf.fit))
  ImpData$Var.Names <- row.names(ImpData)
  ImpData$Var.Names <- factor(ImpData$Var.Names, levels = ImpData$Var.Names[order(ImpData$`%IncMSE`)])
  ImpData$Drug <- col_name
  df_ImpData_single <- rbind(df_ImpData_single,ImpData)
}



df_ImpData_sum_single <- df_ImpData_single %>% group_by(Var.Names) %>% 
  summarise(Mean=mean(`%IncMSE`),
            SD=sd(`%IncMSE`), 
            Max=max(`%IncMSE`),
            Which.max=Drug[which.max(`%IncMSE`)])
df_ImpData_sum_single$tval <- df_ImpData_sum_single$Mean/df_ImpData_sum_single$SD


#combo----
start_col <- 42
end_col <- ncol(rf_df_flow_combo_na_omit)
df_ImpData_combo<-c()
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
  ImpData <- as.data.frame(importance(rf.fit))
  ImpData$Var.Names <- row.names(ImpData)
  ImpData$Var.Names <- factor(ImpData$Var.Names, levels = ImpData$Var.Names[order(ImpData$`%IncMSE`)])
  ImpData$Drug <- col_name
  df_ImpData_combo <- rbind(df_ImpData_combo,ImpData)
}


df_ImpData_sum_combo <- df_ImpData_combo %>% group_by(Var.Names) %>% 
  summarise(Mean=mean(`%IncMSE`),
            SD=sd(`%IncMSE`), 
            Max=max(`%IncMSE`),
            Which.max=Drug[which.max(`%IncMSE`)])
df_ImpData_sum_combo$tval <- df_ImpData_sum_combo$Mean/df_ImpData_sum_combo$SD


#single H0----
start_col <- 42
end_col <- ncol(rf_df_flow_single_na_omit)
df_ImpData_single_h0<-c()
for (col_index in start_col:end_col) {
  set.seed(col_index)
  # Get the column name for the current col_index
  col_name <- colnames(rf_df_flow_single_na_omit)[col_index]
  ids <- sample(1:nrow(rf_df_flow_single_na_omit))
  # Create random forest model
  rf.fit <- randomForest(rf_df_flow_single_na_omit[ids, col_index] ~ ., 
                         data = rf_df_flow_single_na_omit[, -c(42:ncol(rf_df_flow_single_na_omit))], 
                         ntree = 500,
                         keep.forest = FALSE, 
                         importance = TRUE)
  
  # Extract variable importance data
  ImpData <- as.data.frame(importance(rf.fit))
  ImpData$Var.Names <- row.names(ImpData)
  ImpData$Var.Names <- factor(ImpData$Var.Names, levels = ImpData$Var.Names[order(ImpData$`%IncMSE`)])
  ImpData$Drug <- col_name
  df_ImpData_single_h0 <- rbind(df_ImpData_single_h0,ImpData)
}


#combo H0----
start_col <- 42
end_col <- ncol(rf_df_flow_combo_na_omit)
df_ImpData_combo_h0 <-c()
for (col_index in start_col:end_col) {
  set.seed(col_index)
  # Get the column name for the current col_index
  col_name <- colnames(rf_df_flow_combo_na_omit)[col_index]
  ids <- sample(1:nrow(rf_df_flow_combo_na_omit))
  # Create random forest model
  rf.fit <- randomForest(rf_df_flow_combo_na_omit[ids, col_index] ~ ., 
                         data = rf_df_flow_combo_na_omit[, -c(42:ncol(rf_df_flow_combo_na_omit))], 
                         ntree = 500,
                         keep.forest = FALSE, 
                         importance = TRUE)
  
  # Extract variable importance data
  ImpData <- as.data.frame(importance(rf.fit))
  ImpData$Var.Names <- row.names(ImpData)
  ImpData$Var.Names <- factor(ImpData$Var.Names, levels = ImpData$Var.Names[order(ImpData$`%IncMSE`)])
  ImpData$Drug <- col_name
  df_ImpData_combo_h0 <- rbind(df_ImpData_combo_h0,ImpData)
}


#plots for single drugs----
df_ImpData_sum_single <- df_ImpData_sum_single[order(df_ImpData_sum_single$Mean),]
df_ImpData_single$Var.Names <- factor(df_ImpData_single$Var.Names, levels=rev(df_ImpData_sum_single$Var.Names) )

ggplot(df_ImpData_single)+
  geom_boxplot(aes(Var.Names, `%IncMSE`), fill = "#93b9ed", outlier.colour = "darkgrey", lwd=0.4)+
  geom_hline(yintercept = quantile(df_ImpData_single_h0$`%IncMSE`, 0.95), lty=2)+
  theme_classic()+
  labs(x="")+
  ggtitle("% Increse MSE per variable real data single agents")+
  ylim(-5, 20)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

df_ImpData_single_h0$Var.Names <- factor(df_ImpData_single_h0$Var.Names, levels=rev(df_ImpData_sum_single$Var.Names) )
ggplot(df_ImpData_single_h0)+
  geom_boxplot(aes(Var.Names, `%IncMSE`), fill = "#93b9ed", outlier.colour = "darkgrey", lwd=0.4)+
  theme_classic()+
  labs(x="")+
  ggtitle("% Increse MSE per variable randomized data single agents")+
  ylim(-5, 20)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

df_ImpData_single$p.value <- NA
for(i in 1:nrow(df_ImpData_single)){
  df_ImpData_single$p.value[i] <-  sum(df_ImpData_single_h0$`%IncMSE` >= df_ImpData_single$`%IncMSE`[i])/c(1+nrow(df_ImpData_single_h0))
}

df_ImpData_single$p.adj <- p.adjust(df_ImpData_single$p.value, method = "BH")

df_ImpData_single_signif <- df_ImpData_single %>% group_by(Var.Names) %>% 
  summarise(pval0.01=sum(p.value<0.01),
            pval0.05=sum(p.value<0.05),
            padj0.01=sum(p.adj<0.01),
            padj0.05=sum(p.adj<0.05))

df_ImpData_single_signif <- df_ImpData_single_signif[order(df_ImpData_single_signif$pval0.01),]
df_ImpData_single_signif$Var.Names <- factor(df_ImpData_single_signif$Var.Names, levels = rev(df_ImpData_single_signif$Var.Names))
ggplot(df_ImpData_single_signif)+
  geom_col(aes(Var.Names, pval0.01), fill="#db9d79")+
  theme_classic()+
  labs(x="", y="#p-value<0.01")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle("Number of single drugs with significant p-value (0.01)")+
  ylim(0,45)

df_ImpData_single_signif <- df_ImpData_single_signif[order(df_ImpData_single_signif$pval0.05),]
df_ImpData_single_signif$Var.Names <- factor(df_ImpData_single_signif$Var.Names, levels = rev(df_ImpData_single_signif$Var.Names))
ggplot(df_ImpData_single_signif)+
  geom_col(aes(Var.Names, pval0.05),  fill="#db9d79")+
  theme_classic()+
  labs(x="", y="#p-value<0.05")+
  ggtitle("Number of single drugs with significant p-value (0.05)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,45)


# plots for combo plots----
df_ImpData_sum_combo <- df_ImpData_sum_combo[order(df_ImpData_sum_combo$Mean),]
df_ImpData_combo$Var.Names <- factor(df_ImpData_combo$Var.Names, levels=rev(df_ImpData_sum_combo$Var.Names) )
ggplot(df_ImpData_combo)+
  geom_boxplot(aes(Var.Names, `%IncMSE`),  fill = "#93b9ed", outlier.colour = "darkgrey", lwd=0.4)+
  geom_hline(yintercept = quantile(df_ImpData_combo_h0$`%IncMSE`, 0.95), lty=2)+
  theme_classic()+
  labs(x="")+
  ylim(-5, 20)+
  ggtitle("% Increse MSE per variable real data drug combinations")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

df_ImpData_combo_h0$Var.Names <- factor(df_ImpData_combo_h0$Var.Names, levels=rev(df_ImpData_sum_combo$Var.Names) )
ggplot(df_ImpData_combo_h0)+
  geom_boxplot(aes(Var.Names, `%IncMSE`), fill = "#93b9ed", outlier.colour = "darkgrey", lwd=0.4)+
  theme_classic()+
  labs(x="")+
  ylim(-5, 20)+
  ggtitle("% Increse MSE per variable randomized data drug combinations")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


df_ImpData_combo$p.value <- NA
for(i in 1:nrow(df_ImpData_combo)){
  df_ImpData_combo$p.value[i] <-  sum(df_ImpData_combo_h0$`%IncMSE` >= df_ImpData_combo$`%IncMSE`[i])/c(1+nrow(df_ImpData_combo_h0))
}

df_ImpData_combo$p.adj <- p.adjust(df_ImpData_combo$p.value, method = "BH")

df_ImpData_combo_signif <- df_ImpData_combo %>% group_by(Var.Names) %>% 
  summarise(pval0.01=sum(p.value<0.01),
            pval0.05=sum(p.value<0.05),
            padj0.01=sum(p.adj<0.01),
            padj0.05=sum(p.adj<0.05))

df_ImpData_combo_signif <- df_ImpData_combo_signif[order(df_ImpData_combo_signif$pval0.01),]
df_ImpData_combo_signif$Var.Names <- factor(df_ImpData_combo_signif$Var.Names, levels = rev(df_ImpData_combo_signif$Var.Names))
ggplot(df_ImpData_combo_signif)+
  geom_col(aes(Var.Names, pval0.01), fill="#db9d79")+
  theme_classic()+
  labs(x="", y="#p-value<0.01")+
  ggtitle("Number of drug combinations with significant p-value (0.01)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,45)

df_ImpData_combo_signif <- df_ImpData_combo_signif[order(df_ImpData_combo_signif$pval0.05),]
df_ImpData_combo_signif$Var.Names <- factor(df_ImpData_combo_signif$Var.Names, levels = rev(df_ImpData_combo_signif$Var.Names))
ggplot(df_ImpData_combo_signif)+
  geom_col(aes(Var.Names, pval0.05), fill="#db9d79")+
  theme_classic()+
  labs(x="", y="#p-value<0.05")+
  ggtitle("Number of drug combinations with significant p-value (0.05)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,45)


#save files----
write.csv(df_ImpData_single_signif, "Rf_significance_single_realdata_vs_randomized_177.csv")
write.csv(df_ImpData_combo_signif, "Rf_significance_combo_realdata_vs_randomized_177.csv")
