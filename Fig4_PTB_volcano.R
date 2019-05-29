#Jen Fettweis, updated 2019
#Plot DESeq2 results for Figure 4, metatranscriptomic data
library(ggplot2)
library(ggpubr)
library(plyr)
library(gridExtra)
library(grid)
library("data.table")
res <- fread("/PATH/Fig4_CDS_tRNA_oneSample_122_wmts_expanded_091918.txt", header=TRUE, sep = "\t",quote="")
#Remove those with NAs from DESeq2
res2<- subset(res,!(res$padj=="NA"))
res2$threshold <- "Other"
res2$threshold = ifelse(res2$stirrups=="TM7_OTU-H1" & res2$padj>0.05,"TM7_OTU-H1, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="TM7_OTU-H1" & res2$padj<=0.05,"TM7_OTU-H1, padj≤0.05",res2$threshold)
res2$threshold = ifelse(is.na(res2$threshold),"Other",res2$threshold)
plot_TM7_OTUH1 <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("TM7 OTU-H1")+
  theme(plot.title=element_text(hjust=0.5,size=12))+
  theme(axis.title = element_text(size=12))+
geom_point(data = res2,aes(colour = threshold),size = 1.5, stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'TM7_OTU-H1, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'TM7_OTU-H1, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("TM7_OTU-H1, padj>0.05"="#f4a582","TM7_OTU-H1, padj≤0.05" = "#b2182b", "Other"="#d6d6d6"))+
  #geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Lachnospiraceae_BVAB1" & res2$padj>0.05,"Lachnospiraceae_BVAB1, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Lachnospiraceae_BVAB1" & res2$padj<=0.05,"Lachnospiraceae_BVAB1, padj≤0.05",res2$threshold)
plot_Lachnospiraceae_BVAB1 <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Lachnospiraceae BVAB1")+
  theme(plot.title=element_text(hjust=0.5,size=12))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1.5, stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Lachnospiraceae_BVAB1, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Lachnospiraceae_BVAB1, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Lachnospiraceae_BVAB1, padj>0.05"="#f4a582","Lachnospiraceae_BVAB1, padj≤0.05" = "#b2182b", "Other"="#d6d6d6"))+
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Lactobacillus_crispatus_cluster" & res2$padj>0.05,"Lactobacillus_crispatus_cluster, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Lactobacillus_crispatus_cluster" & res2$padj<=0.05,"Lactobacillus_crispatus_cluster, padj≤0.05",res2$threshold)
plot_Lactobacillus_crispatus_cluster <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Lactobacillus crispatus cluster")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1,  stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Lactobacillus_crispatus_cluster, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Lactobacillus_crispatus_cluster, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Lactobacillus_crispatus_cluster, padj>0.05"="#92c5de","Lactobacillus_crispatus_cluster, padj≤0.05" = "#2166ac", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Prevotella_cluster2" & res2$padj>0.05,"Prevotella_cluster2, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Prevotella_cluster2" & res2$padj<=0.05,"Prevotella_cluster2, padj≤0.05",res2$threshold)
plot_Prevotella_cluster2 <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Prevotella cluster 2")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1,  stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Prevotella_cluster2, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Prevotella_cluster2, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Prevotella_cluster2, padj>0.05"="#f4a582","Prevotella_cluster2, padj≤0.05" = "#b2182b", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Sneathia_amnii" & res2$padj>0.05,"Sneathia_amnii, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Sneathia_amnii" & res2$padj<=0.05,"Sneathia_amnii, padj≤0.05",res2$threshold)
plot_Sneathia_amnii <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Sneathia amnii")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1,  stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Sneathia_amnii, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Sneathia_amnii, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Sneathia_amnii, padj>0.05"="#f4a582","Sneathia_amnii, padj≤0.05" = "#b2182b", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Dialister_cluster51" & res2$padj>0.05,"Dialister_cluster51, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Dialister_cluster51" & res2$padj<=0.05,"Dialister_cluster51, padj≤0.05",res2$threshold)
plot_Dialister_cluster51 <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Dialister cluster 51")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
    geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1,  stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Dialister_cluster51, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Dialister_cluster51, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Dialister_cluster51, padj>0.05"="#f4a582","Dialister_cluster51, padj≤0.05" = "#b2182b", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Prevotella_amnii" & res2$padj>0.05,"Prevotella_amnii, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Prevotella_amnii" & res2$padj<=0.05,"Prevotella_amnii, padj≤0.05",res2$threshold)
plot_Prevotella_amnii <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Prevotella amnii")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
    geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1,  stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Prevotella_amnii, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Prevotella_amnii, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Prevotella_amnii, padj>0.05"="#f4a582","Prevotella_amnii, padj≤0.05" = "#b2182b", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Sneathia_sanguinegens" & res2$padj>0.05,"Sneathia_sanguinegens, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Sneathia_sanguinegens" & res2$padj<=0.05,"Sneathia_sanguinegens, padj≤0.05",res2$threshold)
plot_Sneathia_sanguinegens <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Sneathia sanguinegens")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
    geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1,  stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Sneathia_sanguinegens, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Sneathia_sanguinegens, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Sneathia_sanguinegens, padj>0.05"="#f4a582","Sneathia_sanguinegens, padj≤0.05" = "#b2182b", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Aerococcus_christensenii" & res2$padj>0.05,"Aerococcus_christensenii, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Aerococcus_christensenii" & res2$padj<=0.05,"Aerococcus_christensenii, padj≤0.05",res2$threshold)
plot_Aerococcus_christensenii <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Aerococcus christensenii")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1,  stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Aerococcus_christensenii, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Aerococcus_christensenii, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Aerococcus_christensenii, padj>0.05"="#f4a582","Aerococcus_christensenii, padj≤0.05" = "#b2182b", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Dialister_micraerophilus" & res2$padj>0.05,"Dialister_micraerophilus, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Dialister_micraerophilus" & res2$padj<=0.05,"Dialister_micraerophilus, padj≤0.05",res2$threshold)
plot_Dialister_micraerophilus <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Dialister micraerophilus")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1,  stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Dialister_micraerophilus, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Dialister_micraerophilus, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Dialister_micraerophilus, padj>0.05"="#f4a582","Dialister_micraerophilus, padj≤0.05" = "#b2182b", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Coriobacteriaceae_OTU27" & res2$padj>0.05,"Coriobacteriaceae_OTU27, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Coriobacteriaceae_OTU27" & res2$padj<=0.05,"Coriobacteriaceae_OTU27, padj≤0.05",res2$threshold)
plot_Coriobacteriaceae_OTU27 <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Coriobacteriaceae OTU27")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1,  stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Coriobacteriaceae_OTU27, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Coriobacteriaceae_OTU27, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Coriobacteriaceae_OTU27, padj>0.05"="#f4a582","Coriobacteriaceae_OTU27, padj≤0.05" = "#b2182b", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Megasphaera_OTU70_type1" & res2$padj>0.05,"Megasphaera_OTU70_type1, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Megasphaera_OTU70_type1" & res2$padj<=0.05,"Megasphaera_OTU70_type1, padj≤0.05",res2$threshold)
plot_Megasphaera_OTU70_type1 <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Megasphaera OTU70 type1")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1.5, stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Megasphaera_OTU70_type1, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Megasphaera_OTU70_type1, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Megasphaera_OTU70_type1, padj>0.05"="#f4a582","Megasphaera_OTU70_type1, padj≤0.05" = "#b2182b", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Lactobacillus_iners" & res2$padj>0.05,"Lactobacillus_iners, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Lactobacillus_iners" & res2$padj<=0.05,"Lactobacillus_iners, padj≤0.05",res2$threshold)
plot_Lactobacillus_iners <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Lactobacillus iners")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1.5, stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Lactobacillus_iners, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Lactobacillus_iners, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Lactobacillus_iners, padj>0.05"="#a1d99b","Lactobacillus_iners, padj≤0.05" = "#238b45", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Lactobacillus_gasseri_cluster" & res2$padj>0.05,"Lactobacillus_gasseri_cluster, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Lactobacillus_gasseri_cluster" & res2$padj<=0.05,"Lactobacillus_gasseri_cluster, padj≤0.05",res2$threshold)
plot_Lactobacillus_gasseri_cluster <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Lactobacillus gasseri cluster")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1.5, stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Lactobacillus_gasseri_cluster, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Lactobacillus_gasseri_cluster, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Lactobacillus_gasseri_cluster, padj>0.05"="#a1d99b","Lactobacillus_gasseri_cluster, padj≤0.05" = "#238b45", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Lactobacillus_jensenii" & res2$padj>0.05,"Lactobacillus_jensenii, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Lactobacillus_jensenii" & res2$padj<=0.05,"Lactobacillus_jensenii, padj≤0.05",res2$threshold)
plot_Lactobacillus_jensenii <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Lactobacillus jensenii")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1.5, stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Lactobacillus_jensenii, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Lactobacillus_jensenii, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Lactobacillus_jensenii, padj>0.05"="#a1d99b","Lactobacillus_jensenii, padj≤0.05" = "#238b45", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

res2$threshold = "Other"
res2$threshold = ifelse(res2$stirrups=="Gardnerella_vaginalis" & res2$padj>0.05,"Gardnerella_vaginalis, padj>0.05",res2$threshold)
res2$threshold = ifelse(res2$stirrups=="Gardnerella_vaginalis" & res2$padj<=0.05,"Gardnerella_vaginalis, padj≤0.05",res2$threshold)
plot_Gardnerella_vaginalis <- ggplot(res2, aes(x=log2FoldChange, y=-log10(padj))) +
  ggtitle("Gardnerella vaginalis")+
  theme(plot.title=element_text(hjust=0.5,size=12,face="italic"))+
  theme(axis.title = element_text(size=12))+
  geom_point(data = subset(res2,threshold !="NA"),aes(colour = threshold),size = 1.5, stroke = 0) +
  xlim(-32,32)+
  ylim(-5,40)+
  geom_point(data = subset(res2, threshold == 'Gardnerella_vaginalis, padj>0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  geom_point(data = subset(res2, threshold == 'Gardnerella_vaginalis, padj≤0.05'),aes(colour = threshold),size = 1.5, stroke = 0)+
  scale_colour_manual(values = c("Gardnerella_vaginalis, padj>0.05"="#a1d99b","Gardnerella_vaginalis, padj≤0.05" = "#238b45", "Other"="#d6d6d6"))+
  geom_text(aes(label = "")) +
  annotate("text", label = "↑Preterm", x = 5, y = -3, size = 5)+
  annotate("text", label = "↑Term", x = -5, y = -3, size = 5)+
  theme(legend.position="none")

merged_plot <- ggarrange(plot_Prevotella_cluster2,plot_Sneathia_amnii,plot_Lachnospiraceae_BVAB1,plot_TM7_OTUH1,plot_Dialister_cluster51,plot_Prevotella_amnii,plot_Sneathia_sanguinegens,plot_Dialister_micraerophilus,plot_Coriobacteriaceae_OTU27,plot_Aerococcus_christensenii,plot_Megasphaera_OTU70_type1,plot_Lactobacillus_crispatus_cluster,ncol=4,nrow=3)
merged_plot
#Export as 2100 x 800 pixels
#Panel_WMTS_091918_2100x800.png
merged_plot_extended <- ggarrange(plot_Prevotella_cluster2,plot_Sneathia_amnii,plot_Lachnospiraceae_BVAB1,plot_TM7_OTUH1,plot_Dialister_cluster51,plot_Prevotella_amnii,plot_Sneathia_sanguinegens,plot_Dialister_micraerophilus,plot_Coriobacteriaceae_OTU27,plot_Aerococcus_christensenii,plot_Megasphaera_OTU70_type1,plot_Lactobacillus_crispatus_cluster,plot_Lactobacillus_jensenii,plot_Lactobacillus_gasseri_cluster,plot_Lactobacillus_iners,plot_Gardnerella_vaginalis,ncol=4,nrow=4)
merged_plot_extended 
#Export as 2100 x 1067 pixels
#Panel_WMTS_extended_091918_2100x1067.png


PTB <- count(subset(res2,as.numeric(padj)<=0.05 & log2FoldChange>0& !(is.na(stirrups))),"stirrups")
PTB <- PTB[order(PTB$freq),]
rownames(PTB) <- PTB$stirrups
PTB$stirrups <-NULL
names(PTB)[names(PTB) == 'freq'] <- '↑ Preterm'

TB <- count(subset(res2,as.numeric(padj)<=0.05 & log2FoldChange<0 & !(is.na(stirrups))),"stirrups")
rownames(TB) <- TB$stirrups
TB$stirrups <-NULL
names(TB)[names(TB) == 'freq'] <- '↑ Term'

all_genes <- merge(PTB,TB,by="row.names",all=TRUE)
all_genes[is.na(all_genes)]<- 0
rownames(all_genes) <-gsub('_',' ',all_genes$Row.names)
all_genes$Row.names <- NULL

all_genes_tablegrob <-tableGrob(all_genes[order(-all_genes$`↑ Preterm`),])
grid.draw(all_genes_tablegrob)


