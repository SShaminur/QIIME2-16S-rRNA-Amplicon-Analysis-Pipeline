setwd("C:/Users/USER/Desktop/Nazmul-Bhai/ITS-USA/pheatmap")

library(pheatmap)
library("openxlsx")
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(ggthemes)
library("viridis") 
library(grid)
library(gridExtra)
library(tidyverse)

data <- read.xlsx("heatmap.xlsx", sheet = "data", colNames =TRUE, rowNames = TRUE)
my_gene_col <- read.xlsx("heatmap.xlsx", sheet = "data_r", colNames =TRUE, rowNames = TRUE)
my_sample_col <- read.xlsx("heatmap.xlsx", sheet = "data_c", colNames =TRUE, rowNames = TRUE)

#Log transformation
pseudocount <- 1
p1 <- log10(data[, -1] + pseudocount)

#Percentage transformation
to_percent <- function(column) {
  prop.table(column) * 100
}

# Apply the function to the numeric columns
#df_percent <- data
p <- to_percent(data)


#pheatmap(data, border_color="gray80", cellwidth = 60, cellheight = 12, annotation_row = my_gene_col, cluster_cols = T, cluster_rows = T, display_numbers = F, fontsize_number=10, number_format = "%.2f", number_color = "black", color = colorRampPalette(c("white", "blue", "red"))(100))
#pheatmap(data, border_color="gray80", cellwidth = 30, cellheight = 12, annotation_row = my_gene_col, cluster_cols = T, cluster_rows = T, display_numbers = F, fontsize_number=10, number_format = "%.2f", number_color = "black", color = colorRampPalette(c("white", "blue", "red"))(100))
#pheatmap(data, border_color="gray80", cellwidth = 30, cellheight = 12, annotation_row = my_gene_col, annotation_col = my_sample_col, cluster_cols = T, cluster_rows = T, display_numbers = F, fontsize_number=10, number_format = "%.2f", number_color = "black", color = colorRampPalette(c("white", "blue", "red"))(100))
#pheatmap(data1, border_color="gray80", cellwidth = 30, cellheight = 12, annotation_row = my_gene_col, annotation_col = my_sample_col, cluster_cols = T, cluster_rows = F, display_numbers = F, fontsize_number=10, number_format = "%.2f", number_color = "black", color = colorRampPalette(c("white", "blue", "red"))(100))
#pheatmap(data1, border_color="gray80", cellwidth = 30, cellheight = 12, annotation_row = my_gene_col, annotation_col = my_sample_col, cluster_cols = T, cluster_rows = F, display_numbers = F, fontsize_number=10, number_format = "%.2f", number_color = "black", col = brewer.pal(8, 'RdYlGn'))

pheatmap(p1, cutree_rows = 4, cutree_cols = 5, border_color="gray80", cellwidth = 10, cellheight = 10, annotation_row = my_gene_col, annotation_col = my_sample_col, cluster_cols = T, cluster_rows = T, display_numbers = F, fontsize_number=10, number_format = "%.2f", number_color = "black")


###########################################
data <- read.xlsx("merged.xlsx", sheet = "All", colNames =TRUE, rowNames = TRUE)

#https://stackoverflow.com/questions/20041136/avoid-ggplot-sorting-the-x-axis-while-plotting-geom-bar
data$Pathways <- factor(data$Pathways, levels = data$Pathways)

p <- ggplot(data, aes(fill=pathway_class, y=Percentages, x=Pathways)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + scale_x_discrete(labels = data$Pathways) +
  labs(x = "Pathways") + theme_gdocs() + scale_fill_brewer(palette="Set2")


p

####comparison############
data <- read.xlsx("merged.xlsx", sheet = "SW", colNames =TRUE, rowNames = TRUE)

#https://stackoverflow.com/questions/20041136/avoid-ggplot-sorting-the-x-axis-while-plotting-geom-bar
data$Pathways <- factor(data$Pathways, levels = data$Pathways)

#scale_x_discrete(labels = data$Pathways)
#scale_x_discrete(limits=data$Pathways)

p2 <- ggplot(data=data, aes(x= fct_inorder(Pathways), y=percentages, fill=Taste)) +
  geom_bar(stat="identity", position="dodge") + coord_flip() + 
   theme_gdocs() + scale_fill_brewer(palette="Set2") 



p2 



#p + scale_fill_manual(values = c("A" = "blue", "B" = "red", "C" = "green", "D" = "#FF00FF", "E" = "white", "F" = "orange"))
#p + scale_fill_manual(values = c("A" = "blue", "B" = "red"))
#ggplot(data, aes(x = Pathways, y = Percentages, fill = pathway_class)) +
  geom_bar(stat = "identity") +
  labs(x = "Pathways", y = "Percentages", title = "Percentage Distribution by Pathways") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#1f78b4", "#33a02c", "#e31a1c", "#f39e5c","#91b4e1", "#e2ad14")) +  # Adjust colors as needed
  guides(fill = guide_legend(title = "pathway Class")) + scale_x_discrete(breaks = unique(data$Pathways)) +
  coord_flip()

#p <- ggplot(data, aes(fill=pathway_class, y=Percentages, x=Pathways)) + 
  geom_bar(stat="identity") + coord_flip() + scale_x_discrete(labels = data$Pathways) +
  labs(x = "pathway_class") + theme_tufte()


#####################################Percentage-conversion###############
df <- data.frame(
  Category = c("A", "B", "C"),
  Value1 = c(10, 20, 30),
  Value2 = c(40, 50, 60)
)

# Function to convert a column to percentages
to_percent <- function(column) {
  prop.table(column) * 100
}

# Apply the function to the numeric columns
#df_percent <- data
p <- to_percent(data)

# View the DataFrame with percentages
#print(df_percent)

#pheatmap(p, cutree_rows = 2, cutree_cols = 2, border_color="gray80", cellwidth = 10, cellheight = 10, annotation_row = my_gene_col, annotation_col = my_sample_col, cluster_cols = T, cluster_rows = F, display_numbers = F, fontsize_number=10, number_format = "%.2f", number_color = "black")
##############################################################


#install.packages("gridExtra")

F2 <- p + theme(legend.position="FALSE")
F3 <- p2 + theme(legend.position="FALSE")

Bar <- grid.arrange(F2, F3, ncol = 2) 
#grid.arrange(F2, F3, nrow = 2)  

#install.packages("grid")
ph2 <- ph$gtable


grid.arrange(ph2, Bar, ncol = 1)
