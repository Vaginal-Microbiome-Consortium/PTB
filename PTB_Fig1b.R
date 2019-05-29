#Jen Fettweis, 2019
#R code for data resource Figure 1b

library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(forcats)


df_ptb <- read.csv("/PATH/Fig1b_ptb.csv",header = TRUE)
df1_ptb <- subset(df_ptb,Trimester=="First Trimester")
df2_ptb <- subset(df_ptb,Trimester=="Second Trimester")
df3_ptb <- subset(df_ptb,Trimester=="Third Trimester")
df4_ptb <- subset(df_ptb,Trimester=="Postpartum")

df_tbs <- read.csv("/PATH/Fig1b_tbs.csv",header = TRUE)
df1_tbs <- subset(df_tbs,Trimester=="First Trimester")
df2_tbs <- subset(df_tbs,Trimester=="Second Trimester")
df3_tbs <- subset(df_tbs,Trimester=="Third Trimester")
df4_tbs <- subset(df_tbs,Trimester=="Postpartum")

df_other <- read.csv("/PATH/Fig1b_other.csv",header = TRUE)
df1_other <- subset(df_other,Trimester=="First Trimester")
df2_other <- subset(df_other,Trimester=="Second Trimester")
df3_other <- subset(df_other,Trimester=="Third Trimester")
df4_other <- subset(df_other,Trimester=="Postpartum")

#Parse PTB file
df1_ptb$Body_Site_Omic <- factor(df1_ptb$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal MGS","maternal vaginal MTS","maternal vaginal cytokines","maternal buccal 16S", "maternal rectal 16S", " ","infant buccal 16S","infant rectal 16S", "infant stool 16S"))
ptb1_plot <- ggplot(data=df1_ptb, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f", "maternal vaginal MGS" = "#ef6548","maternal vaginal MTS" = "#fc8d59","maternal vaginal cytokines"="#fdd49e","maternal buccal 16S" = "#bcbd22","maternal rectal 16S" = "#3690c0", " "="white", "infant buccal 16S"="#bcbd22","infant buccal cytokines" ="#dbdb8d","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,10.58), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

df2_ptb$Body_Site_Omic <- factor(df2_ptb$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal MGS","maternal vaginal MTS","maternal vaginal cytokines","maternal buccal 16S", "maternal rectal 16S", " ","infant buccal 16S","infant rectal 16S", "infant stool 16S"))
ptb2_plot <- ggplot(data=df2_ptb, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f", "maternal vaginal MGS" = "#ef6548","maternal vaginal MTS" = "#fc8d59","maternal vaginal cytokines"="#fdd49e","maternal buccal 16S" = "#bcbd22","maternal rectal 16S" = "#3690c0", " "="white", "infant buccal 16S"="#bcbd22","infant buccal cytokines" ="#dbdb8d","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,10.58), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

df3_ptb$Body_Site_Omic <- factor(df3_ptb$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal MGS","maternal vaginal MTS","maternal vaginal cytokines","maternal buccal 16S", "maternal rectal 16S", " ","infant buccal 16S","infant rectal 16S", "infant stool 16S"))
ptb3_plot <- ggplot(data=df3_ptb, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f", "maternal vaginal MGS" = "#ef6548","maternal vaginal MTS" = "#fc8d59","maternal vaginal cytokines"="#fdd49e","maternal buccal 16S" = "#bcbd22","maternal rectal 16S" = "#3690c0", " "="white", "infant buccal 16S"="#bcbd22","infant buccal cytokines" ="#dbdb8d","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,10.58), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

df4_ptb$Body_Site_Omic <- factor(df4_ptb$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal MGS","maternal vaginal MTS","maternal vaginal cytokines","maternal buccal 16S", "maternal rectal 16S", " ","infant buccal 16S","infant rectal 16S", "infant stool 16S"))
ptb4_plot <- ggplot(data=df4_ptb, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f", "maternal vaginal MGS" = "#ef6548","maternal vaginal MTS" = "#fc8d59","maternal vaginal cytokines"="#fdd49e","maternal buccal 16S" = "#bcbd22","maternal rectal 16S" = "#3690c0", " "="white", "infant buccal 16S"="#bcbd22","infant buccal cytokines" ="#dbdb8d","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,10.58), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

#Parse TB file
df1_tbs$Body_Site_Omic <- factor(df1_tbs$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal MGS","maternal vaginal MTS","maternal vaginal cytokines","maternal buccal 16S","maternal buccal cytokines", "maternal rectal 16S", "maternal plasma cytokines"," ","infant buccal 16S","infant buccal cytokines","infant rectal 16S", "infant stool 16S","infant stool MGS"))
tbs1_plot <- ggplot(data=df1_tbs, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f", "maternal vaginal MGS" = "#ef6548","maternal vaginal MTS" = "#fc8d59","maternal vaginal cytokines"="#fdd49e","maternal buccal 16S" = "#bcbd22","maternal buccal cytokines" = "#dbdb8d","maternal rectal 16S" = "#3690c0", "maternal plasma cytokines" = "#a6bddb", " "="white", "infant buccal 16S"="#bcbd22","infant buccal cytokines" ="#dbdb8d","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd","infant stool MGS"="#c5b0d5"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,14.6), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

df2_tbs$Body_Site_Omic <- factor(df2_tbs$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal MGS","maternal vaginal MTS","maternal vaginal cytokines","maternal buccal 16S","maternal buccal cytokines", "maternal rectal 16S", "maternal plasma cytokines"," ","infant buccal 16S","infant buccal cytokines","infant rectal 16S", "infant stool 16S","infant stool MGS"))
tbs2_plot <- ggplot(data=df2_tbs, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f", "maternal vaginal MGS" = "#ef6548","maternal vaginal MTS" = "#fc8d59","maternal vaginal cytokines"="#fdd49e","maternal buccal 16S" = "#bcbd22","maternal buccal cytokines" = "#dbdb8d","maternal rectal 16S" = "#3690c0", "maternal plasma cytokines" = "#a6bddb", " "="white", "infant buccal 16S"="#bcbd22","infant buccal cytokines" ="#dbdb8d","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd","infant stool MGS"="#c5b0d5"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,14.6), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

df3_tbs$Body_Site_Omic <- factor(df3_tbs$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal MGS","maternal vaginal MTS","maternal vaginal cytokines","maternal buccal 16S","maternal buccal cytokines", "maternal rectal 16S", "maternal plasma cytokines"," ","infant buccal 16S","infant buccal cytokines","infant rectal 16S", "infant stool 16S","infant stool MGS"))
tbs3_plot <- ggplot(data=df3_tbs, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f", "maternal vaginal MGS" = "#ef6548","maternal vaginal MTS" = "#fc8d59","maternal vaginal cytokines"="#fdd49e","maternal buccal 16S" = "#bcbd22","maternal buccal cytokines" = "#dbdb8d","maternal rectal 16S" = "#3690c0", "maternal plasma cytokines" = "#a6bddb", " "="white", "infant buccal 16S"="#bcbd22","infant buccal cytokines" ="#dbdb8d","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd","infant stool MGS"="#c5b0d5"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,14.6), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

df4_tbs$Body_Site_Omic <- factor(df4_tbs$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal MGS","maternal vaginal MTS","maternal vaginal cytokines","maternal buccal 16S","maternal buccal cytokines", "maternal rectal 16S", "maternal plasma cytokines"," ","infant buccal 16S","infant buccal cytokines","infant rectal 16S", "infant stool 16S","infant stool MGS"))
tbs4_plot <- ggplot(data=df4_tbs, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f", "maternal vaginal MGS" = "#ef6548","maternal vaginal MTS" = "#fc8d59","maternal vaginal cytokines"="#fdd49e","maternal buccal 16S" = "#bcbd22","maternal buccal cytokines" = "#dbdb8d","maternal rectal 16S" = "#3690c0", "maternal plasma cytokines" = "#a6bddb", " "="white", "infant buccal 16S"="#bcbd22","infant buccal cytokines" ="#dbdb8d","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd","infant stool MGS"="#c5b0d5"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,14.6), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

#Parse other file
df1_other$Body_Site_Omic <- factor(df1_other$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal cytokines","maternal vaginal lipids","maternal buccal 16S","maternal rectal 16S", "maternal cervical 16S","maternal chest 16S"," ","infant buccal 16S","infant rectal 16S", "infant stool 16S"))
other1_plot <- ggplot(data=df1_other, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f","maternal vaginal cytokines"="#fdd49e","maternal vaginal lipids" = "#fee8c8","maternal buccal 16S" = "#bcbd22","maternal rectal 16S" = "#3690c0", "maternal cervical 16S"="#17becf","maternal chest 16S"= "#98df8a", " "="white", "infant buccal 16S"="#bcbd22","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,11.6), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

df2_other$Body_Site_Omic <- factor(df2_other$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal cytokines","maternal vaginal lipids","maternal buccal 16S","maternal rectal 16S", "maternal cervical 16S","maternal chest 16S"," ","infant buccal 16S","infant rectal 16S", "infant stool 16S"))
other2_plot <- ggplot(data=df2_other, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f","maternal vaginal cytokines"="#fdd49e","maternal vaginal lipids" = "#fee8c8","maternal buccal 16S" = "#bcbd22","maternal rectal 16S" = "#3690c0", "maternal cervical 16S"="#17becf","maternal chest 16S"= "#98df8a", " "="white", "infant buccal 16S"="#bcbd22","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,11.6), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

df3_other$Body_Site_Omic <- factor(df3_other$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal cytokines","maternal vaginal lipids","maternal buccal 16S","maternal rectal 16S", "maternal cervical 16S","maternal chest 16S"," ","infant buccal 16S","infant rectal 16S", "infant stool 16S"))
other3_plot <- ggplot(data=df3_other, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f","maternal vaginal cytokines"="#fdd49e","maternal vaginal lipids" = "#fee8c8","maternal buccal 16S" = "#bcbd22","maternal rectal 16S" = "#3690c0", "maternal cervical 16S"="#17becf","maternal chest 16S"= "#98df8a", " "="white", "infant buccal 16S"="#bcbd22","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,11.6), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

df4_other$Body_Site_Omic <- factor(df4_other$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal cytokines","maternal vaginal lipids","maternal buccal 16S","maternal rectal 16S", "maternal cervical 16S","maternal chest 16S"," ","infant buccal 16S","infant rectal 16S", "infant stool 16S"))
other4_plot <- ggplot(data=df4_other, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f","maternal vaginal cytokines"="#fdd49e","maternal vaginal lipids" = "#fee8c8","maternal buccal 16S" = "#bcbd22","maternal rectal 16S" = "#3690c0", "maternal cervical 16S"="#17becf","maternal chest 16S"= "#98df8a", " "="white", "infant buccal 16S"="#bcbd22","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd"))+
  ylim(0, 700)+
  theme_gray()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  geom_hline(yintercept = c(100,200,300,400,500,600,700), color="white", size=1)+
  geom_vline(xintercept = c(0.43,0.46,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,11.6), color="white", size=1)+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0),breaks = c(200, 400,600))

#Get legend
all <- subset(rbind(df_ptb,df_tbs,df_other),Body_Site_Omic!=" ")
all$Body_Site_Omic <- as.character(all$Body_Site_Omic)
all_sum <- as.data.frame(aggregate(Number_of_Samples~Body_Site_Omic,data=all,sum))
#all_sum <- rbind(all_sum,c(as.character(" "),as.numeric(2500)))
all_sum$Body_Site_Omic <- factor(all_sum$Body_Site_Omic,levels = c("maternal vaginal 16S", "maternal vaginal MGS","maternal vaginal MTS","maternal vaginal cytokines","maternal vaginal lipids","maternal buccal 16S","maternal buccal cytokines","maternal rectal 16S", "maternal plasma cytokines", "maternal cervical 16S","maternal chest 16S","infant buccal 16S","infant buccal cytokines","infant rectal 16S", "infant stool 16S","infant stool MGS"))
all_plot<-ggplot(data=all_sum, aes(x=fct_rev(Body_Site_Omic), y=Number_of_Samples, fill = Body_Site_Omic)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("maternal vaginal 16S" = "#d7301f", "maternal vaginal MGS" = "#ef6548","maternal vaginal MTS" = "#fc8d59","maternal vaginal cytokines"="#fdd49e","maternal vaginal lipids" = "#fee8c8","maternal buccal 16S" = "#bcbd22","maternal buccal cytokines" = "#dbdb8d","maternal rectal 16S" = "#3690c0", "maternal plasma cytokines" = "#a6bddb", "maternal cervical 16S"="#17becf","maternal chest 16S"= "#98df8a"," "="white", "infant buccal 16S"="#bcbd22","infant buccal cytokines" ="#dbdb8d","infant rectal 16S"="#3690c0", "infant stool 16S"="#a58abd","infant stool MGS"="#c5b0d5"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(2.8),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),plot.margin = unit(c(0, 0, 0, 0), "points"),legend.text=element_text(size=rel(1.5)),legend.title=element_blank(),legend.key = element_rect(size = rel(0.8)),legend.box.margin=margin(0,0,0,0),legend.key.size= unit(4,"lines"),legend.key.height = unit(1,"lines"))+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow=6,byrow=FALSE),hjust = 0.5,vjust=0.5)
all_legend<-get_legend(all_plot)
g_legend<-ggdraw(plot_grid(all_legend))

over_fig<-grid.arrange(arrangeGrob(ptb1_plot,ptb2_plot,ptb3_plot,ptb4_plot,nrow=1,ncol=4),arrangeGrob(tbs1_plot,tbs2_plot,tbs3_plot,tbs4_plot,nrow=1,ncol=4),arrangeGrob(other1_plot,other2_plot,other3_plot,other4_plot,nrow=1,ncol=4),heights=c((10+1.6)/10,(13+1.6)/10,(11+1.6)/10))
grid.arrange(over_fig,g_legend,heights=c(1,0.45))
#1270x820