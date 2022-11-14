# Load Packages :
library(tidyverse)
library(circlize)
library(seriation)
library(ComplexHeatmap)
library(edgeR)
library(org.Hs.eg.db)
library(reshape2)
library(Glimma)
library(gridExtra)
library(CEMiTool)
library(clusterProfiler)
library(openxlsx)
library(gplots)
library(ggplot2)
library(factoextra)
library(biomaRt)
library(AnnotationDbi)
#require(hgu95av2.db) 
library(ggrepel)
library(cola)
library(eulerr)
library(VennDiagram)
library(reshape2)
library(devtools)
library(cluster)
library(RColorBrewer)

pheno <- data.frame(patient = paste("LB", 1:2, sep = ""), 
                    Medio = c(rep("LB", 1), rep("SALES", 1)), 
                    Age = c(rep("<60", 1), rep(">60", 1)), 
                    Treatment = c(rep("Chemotherapy_curative_intention", 1), rep("Paliative_treatment", 1))
                    )
pheno$group <- "A"        #R
pheno$group[2:2] <- "B" #NR

annotations <- read.delim("GCF_000017765.1_ASM1776v1_genomic.gff.gz", h = F, skip = 9)
#annotations <- annotations[-which(annotations$V3=="gene"),]
#annotations <- annotations[-which(annotations$V3=="region"),]
annotations1 <- annotations[(annotations$V2=="tRNAscan-SE"),]
annotations1 <- annotations[(annotations$V3=="tRNA"),]

annotations2 <- annotations[(annotations$V3=="ncRNA"),]
annotations2 <- annotations[(annotations$V2=="cmsearch"),]
annotations2 <- annotations2[-which(annotations2$V3=="exon"),]

annotations_ok <- rbind(annotations1, annotations2)
annotations_ok <- annotations_ok[-which(annotations_ok$V3=="sequence_feature"),]
annotations_ok <- annotations_ok[-which(annotations_ok$V3=="riboswitch"),]





annotations_ID <- str_split(annotations_ok$V9, "ECHS_", simplify = T) #lista de cosas, pero por matriz
annotations_ID <- as.data.frame(annotations_ID)
annotations_ID <- str_split(annotations_ID$V4, ";", simplify = T) #lista de cosas, pero por matriz
annotations_ID <- as.data.frame(annotations_ID )
annotations_ID$V1
annotations_product <- str_split(annotations_ok$V9, "product=", simplify = T) #lista de cosas, pero por matriz
annotations_product <- as.data.frame(annotations_product)
annotations_product$V2


rownames(annotations_ok) <- annotations_ID$V1
annotations_ok$V9 <- annotations_product$V2
annotations_ok$V1 <- annotations_ID$V1
#annotations <- annotations[!duplicated(annotations$annotations),]
#annotations <- as.data.frame(annotations)
#annotations <- annotations[!(annotations$annotations==''),]

######### ARCHIVOS .TXT

A <- read.delim("m2-A.txt", sep = "\t", h = F, stringsAsFactors = F, skip = 0)
B <- read.delim("m2-B.txt", sep = "\t", h = F, stringsAsFactors = F, skip = 0)

A$V1 <- str_split(A$V1, "ECHS_", simplify = T)[,2]
A$V1 <- str_split(A$V1, "-", simplify = T)[,1]
A <- A[!duplicated(A$V1),]
#A <- as.data.frame(A)
B$V1 <- str_split(B$V1, "ECHS_", simplify = T)[,2]
B$V1 <- str_split(B$V1, "-", simplify = T)[,1]
B <- B[!duplicated(B$V1),]


annotations_final_a <- merge(x=annotations_ok,y=A,by="V1")
colnames(annotations_final_a) <- c("id","source","type","first","last","none","strand","none","name","A.count")
annotations_final_a$none <- NULL
annotations_final_a$none <- NULL

annotations_final_b <- merge(x=annotations_ok,y=B,by="V1")
colnames(annotations_final_b) <- c("id","source","type","first","last","none","strand","none","name","B.count")
annotations_final_b$none <- NULL
annotations_final_b$none <- NULL

annotations_match <- merge(x=annotations_final_a,y=annotations_final_b,by="id")
annotations_match <- annotations_match[,c("id","name.x","A.count","B.count")]
colnames(annotations_match)<-c("id","name","A","B")
write.table(annotations_match, "ANNOT-MATCH.txt", sep = "\t", row.names = F, quote = F)

##ALL
annotations_all <- merge(x=annotations_final_a,y=annotations_final_b,by="id",all.y=TRUE)
annotations_all <- annotations_all[,c("id","name.y","A.count","B.count")]
colnames(annotations_all)<-c("id","name","A","B")
annotations_all[is.na(annotations_all)] <- 0


write.table(annotations_all, "ANNOT-ALL.txt", sep = "\t", row.names = F, quote = F)



