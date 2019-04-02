setwd('D:\\Peptide prediction\\Anti-Angiogenic Peptides\\backup')

library(reshape2)
library(ggplot2)
library(scales)

data = read.csv("Heatmap DPC.csv", header = TRUE)

ggplot(data, aes(ntree,mtry,fill = MDGI, label = round(MDGI,2))) + 
  geom_tile() + 
  geom_text() +
scale_fill_gradient2(mid= "pink", high= "red", low= "dodgerblue3", midpoint =2,limit = c(0,4),
space = "Lab",name="MDGI") +
  labs(
    x = "", 
    y = "",
    fill = "MDGI") +
theme(panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=15,face="bold"),
        axis.text.x = element_text(angle = 360, vjust = 1, size = 10,hjust = 1),
	axis.text.y = element_text(angle = 360, vjust = 1, size = 10,hjust = 1))
