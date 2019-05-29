#Jen Fettweis,
#R code for data resource Figure 1c
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(forcats)
df_ramsreg <- read.csv("/PATH/ramsreg_resources.csv",header = TRUE)
df_ramsreg$Body_Site<-ordered(df_ramsreg$Body_Site,levels=df_ramsreg$Body_Site)
ramsreg_samples_plot<-ggplot(data=df_ramsreg, aes(x=fct_rev(Body_Site), y=Samples, fill = Body_Site))+
  coord_flip()+
  geom_bar(stat="identity")+
  theme_gray()+
  theme(axis.ticks.y = element_blank(),axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  scale_fill_manual(values = c("maternal vagina" = "#d7301f","maternal buccal mucosa" = "#bcbd22","maternal rectum" = "#3690c0", "maternal blood"="#a6bddb","maternal urine"="#9EDAE5","maternal cervix"="#17becf","maternal chest"="#98df8a","maternal dominant palm"="grey","maternal antecubital fossa"="#C49C94","maternal nares"="#2CA02C","cord blood"="#a6bddb","amniotic fluid"="#FF9896","placental membranes"="#FF7F0E","placental tissue"="#FFBB78","umbilical cord"="#8C564B","infant buccal mucosa"="#bcbd22","infant rectum"="#3690c0", "infant meconium/stool"="#a58abd","infant chest"="#98df8a","infant right palm"="grey","infant nares"="#2CA02C","infant respiratory secretions"="#E377C2"))+
  theme(legend.position="none",axis.text.y = element_blank())+
  scale_y_continuous(limits=c(0, 55000), expand = c(0, 0),breaks=(c(20000,40000)))

ramsreg_visits_plot<-ggplot(data=df_ramsreg, aes(x=fct_rev(Body_Site), y=Visits, fill = Body_Site))+
  coord_flip()+
  geom_bar(stat="identity")+
  theme_gray()+
  theme(axis.ticks.y = element_blank(),axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  scale_fill_manual(values = c("maternal vagina" = "#d7301f","maternal buccal mucosa" = "#bcbd22","maternal rectum" = "#3690c0", "maternal blood"="#a6bddb","maternal urine"="#9EDAE5","maternal cervix"="#17becf","maternal chest"="#98df8a","maternal dominant palm"="grey","maternal antecubital fossa"="#C49C94","maternal nares"="#2CA02C","cord blood"="#a6bddb","amniotic fluid"="#FF9896","placental membranes"="#FF7F0E","placental tissue"="#FFBB78","umbilical cord"="#8C564B","infant buccal mucosa"="#bcbd22","infant rectum"="#3690c0", "infant meconium/stool"="#a58abd","infant chest"="#98df8a","infant right palm"="grey","infant nares"="#2CA02C","infant respiratory secretions"="#E377C2"))+
  theme(legend.position = "none",axis.text.y = element_blank())+
  scale_y_continuous(limits=c(0, 7500), expand = c(0, 0),breaks = c(2000, 4000,6000))


ramsreg_preg_plot<-ggplot(data=df_ramsreg, aes(x=fct_rev(Body_Site), y=Pregnancies, fill = Body_Site))+
  coord_flip()+
  geom_bar(stat="identity")+
  theme_gray()+
  theme(axis.ticks.y = element_blank(),axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  scale_fill_manual(values = c("maternal vagina" = "#d7301f","maternal buccal mucosa" = "#bcbd22","maternal rectum" = "#3690c0", "maternal blood"="#a6bddb","maternal urine"="#9EDAE5","maternal cervix"="#17becf","maternal chest"="#98df8a","maternal dominant palm"="grey","maternal antecubital fossa"="#C49C94","maternal nares"="#2CA02C","cord blood"="#a6bddb","amniotic fluid"="#FF9896","placental membranes"="#FF7F0E","placental tissue"="#FFBB78","umbilical cord"="#8C564B","infant buccal mucosa"="#bcbd22","infant rectum"="#3690c0", "infant meconium/stool"="#a58abd","infant chest"="#98df8a","infant right palm"="grey","infant nares"="#2CA02C","infant respiratory secretions"="#E377C2"))+
  theme(legend.position = "none",axis.text.y = element_blank())+
  scale_y_continuous(limits=c(0, 1700), expand = c(0, 0),breaks = c(500, 1000,1500))

ramsreg_preg_plot_leg<-ggplot(data=df_ramsreg, aes(x=fct_rev(Body_Site), y=Pregnancies, fill = Body_Site))+
  coord_flip()+
  geom_bar(stat="identity")+
  theme_gray()+
  theme(axis.ticks = element_blank(),axis.text.x=element_text(size=rel(1.3),margin=margin(t=0,r=0,l=0,b=0)),axis.title.y = element_blank(),legend.position="none",plot.margin = unit(c(2, 2, 2, 2), "points"),legend.text=element_text(size=rel(0.8)),legend.key.size= unit(0.96,"lines"))+
  scale_fill_manual(values = c("maternal vagina" = "#d7301f","maternal buccal mucosa" = "#bcbd22","maternal rectum" = "#3690c0", "maternal blood"="#a6bddb","maternal urine"="#9EDAE5","maternal cervix"="#17becf","maternal chest"="#98df8a","maternal dominant palm"="grey","maternal antecubital fossa"="#C49C94","maternal nares"="#2CA02C","cord blood"="#a6bddb","amniotic fluid"="#FF9896","placental membranes"="#FF7F0E","placental tissue"="#FFBB78","umbilical cord"="#8C564B","infant buccal mucosa"="#bcbd22","infant rectum"="#3690c0", "infant meconium/stool"="#a58abd","infant chest"="#98df8a","infant right palm"="grey","infant nares"="#2CA02C","infant respiratory secretions"="#E377C2"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(size=rel(2.8),margin=margin(t=0,r=0,l=0,b=0)),axis.ticks = element_blank(),axis.title.y = element_blank(),plot.margin = unit(c(0, 0, 0, 0), "points"),legend.text=element_text(size=rel(1.5)),legend.title=element_blank(),legend.key = element_rect(size = rel(0.8)),legend.box.margin=margin(0,0,0,0),legend.key.size= unit(4,"lines"),legend.key.height = unit(1,"lines"))+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow=5,byrow=FALSE),hjust = 0.5,vjust=0.5)+
  scale_y_continuous(limits=c(0, 1700), expand = c(0, 0),breaks = c(500, 1000,1500))


rr_legend<-ggdraw(plot_grid(get_legend(ramsreg_preg_plot_leg)))

rr_over_fig <- grid.arrange(ramsreg_samples_plot,ramsreg_visits_plot,ramsreg_preg_plot, nrow=1,ncol=3, widths= c(1,1,1))
