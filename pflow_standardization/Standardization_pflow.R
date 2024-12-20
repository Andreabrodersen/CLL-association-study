#library----
library(tidyverse)
library(readxl)
library(janitor)
library(reshape2)
library(dplyr)

#read and format raw pflow file: 
Flow_ctrl_long<-read_csv("Pflow_177_raw.csv")

colnames_flow<-colnames(Flow_ctrl_long)
Flow_ctrl_long<-data.frame(Flow_ctrl_long)
rownames_flow<-as.character(Flow_ctrl_long$...1)
Flow_ctrl_long$...1<-NULL
rownames(Flow_ctrl_long)<-rownames_flow
colnames_flow<-colnames_flow[-1]
colnames(Flow_ctrl_long)<-colnames_flow

#transform data----
#making flow sub IGG
flow_ctrl_nonoise<-Flow_ctrl_long
flow_ctrl_nonoise[, 1:32] <- flow_ctrl_nonoise[, 1:32] - flow_ctrl_nonoise[, 9]
flow_ctrl_nonoise$`IgG Kappa (isotype control)`<-NULL
flow_ctrl_nonoise <- replace(flow_ctrl_nonoise, flow_ctrl_nonoise <0, 0)
flow_ctrl_nonoise_colnames<-colnames(flow_ctrl_nonoise)


#standardizing flow sub IGG
# Standardize across rows
flow_ctrl_nonoise_log <- log10(Flow_ctrl_long) - log10(Flow_ctrl_long$`IgG Kappa (isotype control)`)
flow_ctrl_nonoise_log$`IgG Kappa (isotype control)` <- NULL
flow_ctrl_nonoise_log <- replace(flow_ctrl_nonoise_log, flow_ctrl_nonoise_log <0, 0)
flow_ctrl_nonoise_log[flow_ctrl_nonoise_log == "NaN"] <- NA

flow_ctrl_nonoise_log_scaled <- flow_ctrl_nonoise_log
flow_ctrl_nonoise_log_scaled[,] <- t(apply(flow_ctrl_nonoise_log, 1, scale))



# Plots for inspecting effects of MFI handling----
#raw MFI
flow_ctrl_long_for_plot <- gather(Flow_ctrl_long, key = "variable", value = "value")
p<-ggplot(flow_ctrl_long_for_plot, aes(x = variable, y = value)) +
  geom_boxplot(fill = "white", colour = "#7199A1",outlier.colour = "darkgray") +
  labs(x = "Phospho target", y = "MFI", title = "Boxplot of Phospho-flow Variables") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#noice corrected MFI
flow_ctrl_nonoise_for_plot <- gather(flow_ctrl_nonoise, key = "variable", value = "value")
p<-ggplot(flow_ctrl_nonoise_for_plot, aes(x = variable, y = value)) +
  geom_boxplot(fill = "white", colour = "#7199A1",outlier.colour = "darkgray") +
  labs(x = "Phospho target", y = "noise corrected MFI", title = "Boxplot of Phospho-flow Variables") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#Noise corrected and log transformed MFI
flow_ctrl_nonoise_log_for_plot <- gather(flow_ctrl_nonoise_log, key = "variable", value = "value")
p<-ggplot(flow_ctrl_nonoise_log_for_plot, aes(x = variable, y = value)) +
  geom_boxplot(fill = "white", colour = "#7199A1",outlier.colour = "darkgray") +
  labs(x = "Phospho target", y = "log10(noise corrected MFI)", title = "Boxplot of Phospho-flow Variables") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#noise corrected, log transformed ans scaled MFI 
flow_ctrl_nonoise_log_scaled_for_plot <- gather(flow_ctrl_nonoise_log_scaled, key = "variable", value = "value")
p<-ggplot(flow_ctrl_nonoise_log_scaled_for_plot, aes(x = variable, y = value)) +
  geom_boxplot(fill = "white", colour = "#7199A1",outlier.colour = "darkgray") +
  labs(x = "Phospho target", y = "z-score(log10(noise corrected MFI))", title = "Boxplot of Phospho-flow Variables") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

write.csv(flow_ctrl_nonoise_log_scaled, "Pflow_177_nonoise_log_scaled.csv")

