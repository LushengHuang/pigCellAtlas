library(ggplot2)
require("ggrepel")
library(ggpubr)

#------------------ Set parameters ------------------#
setwd('Pathway')
Plotdata = read.table(file = 'IRX6_celltype_average_expression.csv', sep = ',',header = T)

#------------------ Preprocessing ------------------#
Plotdata$Rank <- rank(-Plotdata$IRX6, ties.method = "first")  # Use negative IRX6 for descending order
Plotdata$ColorGroup <- ifelse(Plotdata$IRX6 >= sort(Plotdata$IRX6, decreasing = TRUE)[10], "Top10", "Others")
Plotdata$Label <- ifelse(Plotdata$IRX6 >= sort(Plotdata$IRX6, decreasing = TRUE)[5], Plotdata$celltype, NA)

# Create the plot
p1=ggplot(data = Plotdata, aes(x = Rank, y = IRX6)) + 
  theme_bw() +
  xlab("751 tissue celltypes") +
  ylab("IRX6 mean expression") +
  geom_point(aes(color = ColorGroup, size = IRX6)) +
  scale_color_manual(values = c("Top10" = "red", "Others" = "grey")) +
  geom_text_repel(aes(label = Label), size = 5, fontface = 'bold', point.padding = 0.5, min.segment.length = 0) +
  theme(plot.title = element_text(size = 14), 
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_text(size=14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 14))
ggsave(paste0("IRX6_expression_in_tissue_celltypes.pdf"),width=5,height=6)

