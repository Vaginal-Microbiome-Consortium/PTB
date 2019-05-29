#Jen Fettweis, 2019
#R code for Fig 5 to make sCCA figure from David Edwards's analysis
library("ggplot2")
library("ggrepel")

df_ptb_cca <- as.data.frame(read.csv('/PATH/Fig5_CCA_coordinates_preterm.csv', sep = ",", header = T))
ptb_cca <- ggplot(df_ptb_cca, aes(x=x, y=y, shape = pch, size = size)) +
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(-1, 1))+
  annotate("path",x=0+1*cos(seq(0,2*pi,length.out=100)),y=0+1*sin(seq(0,2*pi,length.out=100)), alpha = 1, size =0.5, color = "#ffffff")+
  annotate("path",x=0+0.5*cos(seq(0,2*pi,length.out=100)),y=0+0.5*sin(seq(0,2*pi,length.out=100)), alpha = 1, size = 0.5, color = "#ffffff")+
  geom_point(aes(color = Omic), show.legend = FALSE)+ scale_shape_identity()+scale_color_manual(values = c("#053061", "#67001f","#2166ac","#b2182b"))+
  geom_text_repel(aes(label=Shortname, fontface = fontface, color = col), size=4, point.padding=0.1, box.padding = 0.3, show.legend = FALSE)+
  labs(x = "Component 1", y = "Component 2")+
  ggtitle("Preterm") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.y = element_text(size = rel(0.75), angle = 90))+
  theme(axis.title.x = element_text(size = rel(0.75), angle = 0))+
  theme(
    panel.background = element_rect(fill = "#e0e0e0",
                                    colour = "#e0e0e0",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "#f7f7f7"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "#f7f7f7")
  )
ggsave("preterm_cca2-c.png")

df_tb_cca <- as.data.frame(read.csv('/PATH/Fig5_CCA_coordinates_term.csv', sep = ",", header = T))
tb_cca <- ggplot(df_tb_cca, aes(x=x, y=y, shape = pch, size = size))+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(-1, 1))+
  annotate("path",x=0+1*cos(seq(0,2*pi,length.out=100)),y=0+1*sin(seq(0,2*pi,length.out=100)), alpha = 1, size =0.5, color = "#ffffff")+
  annotate("path",x=0+0.5*cos(seq(0,2*pi,length.out=100)),y=0+0.5*sin(seq(0,2*pi,length.out=100)), alpha = 1, size = 0.5, color = "#ffffff")+
  geom_point(aes(color = Omic), show.legend = FALSE)+ scale_shape_identity()+scale_color_manual(values = c("#053061", "#67001f","#2166ac","#b2182b"))+
  geom_text_repel(aes(label=Shortname, fontface = fontface, color = col), size=4, point.padding=0.1, box.padding = 0.3, show.legend = FALSE)+
  labs(x = "Component 1", y = "Component 2")+
  ggtitle("Term") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_reverse( lim=c(1,-1))+
  theme(axis.title.y = element_text(size = rel(0.75), angle = 90))+
  theme(axis.title.x = element_text(size = rel(0.75), angle = 0))+
theme(
  panel.background = element_rect(fill = "#e0e0e0",
                                  colour = "#e0e0e0",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "#f7f7f7"),
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "#f7f7f7")
)

#ggsave("term_cca2-c.png")
ggarrange(tb_cca,ptb_cca,ncol=2)
#ggsave("cca_term_preterm2-w274h135.png", g, dpi=300,units= "mm",width = 274, height = 135)

