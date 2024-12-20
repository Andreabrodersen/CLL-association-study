library(tidyverse)
library(readxl)
library(janitor)
library(reshape2)
library(dplyr)

#READ DSSs for identifying overlap later----
df.dss.z<-read.csv("DSS_177_standardized.csv")
df.dss.no.z<-read.csv("DSS_177_not_standardized.csv")
df.dss.z$X<-str_pad(df.dss.z$X, 3,side="left","0")
rownames(df.dss.z)<-df.dss.z$X
df.dss.z$X<-NULL


#prepping mut:----

mut<-read_excel("VISION_HO141_cytogenetics_TP53_IGHV_14June2022-kopi.xlsx")


#number IGHV column:
mut<-as.data.frame(mut)
rownames(mut)<-mut$`study nr`
mut[mut=="P"]<-"1"
mut[mut=="NP"]<-"0"
mut$IGHV[mut$IGHV == "U"]<- "0"
mut$IGHV[mut$IGHV == "M 3-21"]<-"1"
mut$IGHV[mut$IGHV == "M"] <- "1"
mut$`not eligible` <-NULL
mut$`study nr`<-NULL
mut_rownames<-rownames(mut)
mut_rownames<-str_pad(mut_rownames, 3, pad = "0")
rownames(mut)<-mut_rownames
mut$IGHV<-as.numeric(mut$IGHV)
mut$all_tp53 <- ifelse(mut$`TP53-mutation` == 1 | mut$`17p13 deletion (TP53)` == 1, 1, 0)

mut$IGVH_TP53<- ifelse(mut$IGHV==0 & mut$all_tp53==1, 1, 0)
mut$all_tp53<-NULL


common_patient_ids <- intersect(rownames(df.dss.z), rownames(mut))

# Subset both data frames based on the common patient IDs
Mut_sub <- mut[common_patient_ids, ]

no_na_mut<-na.omit(Mut_sub)



#prep MRD:---- 

mrd<-read_excel("mrd_14June2022_updateVALUESapril2023.xls")
mrd$pid<-str_pad(mrd$pid, 3,side="left","0")
mrd<-data.frame(mrd)
mrd<-mrd %>% drop_na(pid)
rownames(mrd)<-mrd$pid
mrd<-mrd[,c("mrd15bp","pid")]

mrd[mrd=="low"]<-"0"
mrd[mrd=="high"]<-"1"
mrd[mrd=="intermediate"]<-"1"

mrd_sub <- mrd[common_patient_ids, ]
mrd_sub$pid<-NULL

Mut_sub<-bind_cols(Mut_sub,mrd_sub)



#prep pfs:-----
pfs<-read_csv("177prediced_risk_group.csv")
pfs$ID<-str_pad(pfs$ID, 3,side="left","0")
pfs<-data.frame(pfs)
rownames(pfs)<-pfs$ID
pfs<-pfs[,c("ID","risk")]

pfs[pfs=="low"]<-"0"
pfs[pfs=="high"]<-"1"


pfs_sub <- pfs[common_patient_ids, ]
pfs_sub$ID<-NULL

Mut_sub<-bind_cols(Mut_sub,pfs_sub)


#prep binet----
binet<-read.csv("binet_177patients.csv")
binet$pid<-str_pad(binet$pid, 3,side="left","0")
rownames(binet)<-binet$pid
binet$pid<-NULL
binet[binet=="A"]<-"0"
binet[binet=="B"]<-"0"
binet[binet=="C"]<-"1"
binet$binet<-as.numeric(binet$binet)
Mut_sub<-bind_cols(Mut_sub,binet)

write.csv(Mut_sub, "Mutations_177.csv")

