#DSS Standardization


#library
library(tidyverse)
library(reshape2)

#file import
df.dss<-read_csv("DSS_177_not_standardized.csv")

colnames_dss<-colnames(df.dss)
df.dss<-data.frame(df.dss)
rownames_dss<-as.character(df.dss$...1)
df.dss$...1<-NULL
rownames(df.dss)<-rownames_dss
colnames_dss<-colnames_dss[-1]
colnames(df.dss)<-colnames_dss


#standardizing data----- 
# Extracting the numeric values from the data frame
numeric_data <- as.matrix(df.dss)
numeric_data<-t(numeric_data)

# Standardizing across patients (scaling each column)
standardized_data <- scale(numeric_data, scale = TRUE)
standardized_data<-t(standardized_data)

# Converting the result back to a data frame
df.dss.final.z <- as.data.frame(standardized_data)


#write final dfs for DSSs----
write.csv(df.dss.final.z, "DSS_177_standardized.csv")
