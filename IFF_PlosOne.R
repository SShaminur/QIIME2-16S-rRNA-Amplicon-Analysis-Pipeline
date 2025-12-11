setwd("C:/Users/USER/Desktop/IFF-Revision_SR/phyloseq")
##########R Script Written By ###################
#Authore: Md. Shaminur Rahman
#Department of Microbiology
#Jashore University of Science and Technology
#Email: s.rahman@just.edu.bd


library(tidyverse)
library(ggrepel) # for offset labels
library(ggtree) # for visualizing phylogenetic trees
library(ape)
library(microbiome)
library(phyloseq)
library("openxlsx")
library("DESeq2")
library(dplyr)
library(data.table)
library(DT)
library(microbiomeutilities)
library(microbiome)
library(microbial)  #Normalization
library(ggpubr)
library(vegan)
library(ggplot2)
library(ggthemes)
library(gridExtra) 
library(ggtext)
library(RVAideMemoire)
library(MicrobiomeStat)
library(MicrobiomeSurv)


otu_mat<- read.xlsx("data.xlsx", sheet = "organisms", colNames =TRUE, rowNames = TRUE)
tax_mat<- read.xlsx("data.xlsx", sheet = "Taxonomy_table", colNames =TRUE, rowNames = TRUE)
samples_df <- read.xlsx("data.xlsx", sheet = "Metadata", colNames =TRUE, rowNames = TRUE)

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
physeq <- phyloseq(OTU, TAX, samples)

######Normalization######################################
total = sum(sample_sums(physeq)) #TSS

standf = function(x, t=total) round(t * (x / sum(x)))
physeq = transform_sample_counts(physeq, standf)

#Alpha Diversity
comps <- make_pairs(sample_data(physeq)$Type)
comps

rich = estimate_richness(physeq, measures= c("Observed", "Chao1", "Shannon", "InvSimpson"))
write.csv(rich, file = "alpha_diversity_results.csv")
alpha_meas = c("Observed", "Chao1", "Shannon", "InvSimpson")
head(rich)
write.csv(rich, "alpha_diversity_measures.csv")

p <- plot_richness(physeq, "Type", measures=alpha_meas)
p

p5 <- p + geom_boxplot(data=p$data, 
                       aes(x=Type, fill=Type), alpha=0.3,
                       outlier.size = 4, outlier.shape = 2, 
                       outlier.colour = "red") + 
  geom_jitter() + theme_bw() + rremove("x.text")

p5

p7 <- p5 + stat_compare_means(
  comparisons=comps,
  label = "p.signif", method = "wilcox.test", #method = "wilcox.test", "t.test",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")
  )
)

p7 <- p7 + scale_fill_manual(values=c("#EB5406", "#4EE2EC", "#d98880", "#a9cce3", "#d2b4de"))

p7

#Fig.1 exported as PDF

########Kruskal_test################
my_factors = data.frame(sample_data(physeq))
kruskal.test(rich$Chao1 ~ my_factors$Type)

kruskal.sh = kruskal.test(rich$Shannon ~ sample_data(physeq)$Type)
summary(kruskal.sh)
kruskal.sh



##PERMANOVA###################################################
set.seed(999)
abrel_bray <- phyloseq::distance(physeq, method = "bray")

my_factors <- data.frame(sample_data(physeq))
adonis2(abrel_bray ~ Type, data = my_factors)

#PCoA

ord1_1 = ordinate(physeq, method="PCoA", distance = "bray")


beta1_1 <- plot_ordination(physeq, ord1_1, color = "Type", shape= "NA") + 
  geom_point(size= 5, alpha = 0.5) + 
  stat_ellipse(aes(Type=Type)) + theme_bw() + scale_colour_few() + 
  theme(text = element_text(size = 15)) +
  theme(axis.line = element_line(colour = 'black', size = .5)) + theme(legend.position = "right") + 
  annotate("text", x=.8, y=.6, label= "PCoA: Bray", col="black", size=4, parse=TRUE) +
  annotate("text", x=.8, y=.50, label= "PERMANOVA", col="black", size=3, parse=TRUE) +
  annotate("text", x=.8, y=.40, label= "R2=0.473", col="black", size=3, parse=TRUE) +
  annotate("text", x=.8, y=.30, label= "p=0.001", col="black", size=3, parse=TRUE) 

beta1_1 <- beta1_1 + scale_color_manual(values=c("#EB5406", "#4EE2EC", "#d98880", "#a9cce3", "#d2b4de"))

beta1_1

#Fig.2 Exported

####Top genus and phylum statitics#########################

mycols <- c("#EB5406", "#4EE2EC", "#d98880", "#a9cce3", "#d2b4de")



p <- plot_taxa_boxplot(physeq,
                       taxonomic.level = "Genus",
                       top.otu = 26, 
                       group = "Type",
                       add.violin= FALSE,
                       title = "Top twenty five Genus", 
                       keep.other = FALSE,
                       group.order = c("Chicken_gut",	"Droppings",	"Feed",	"Fish_Intestine",	"Sediment"),
                       dot.opacity = .4,
                       box.opacity = .4,
                       group.colors = mycols,
                       dot.size = 2) + theme_biome_utils() + rremove("x.text")

p

physeq.f <- format_to_besthit(physeq)

comps <- make_pairs(sample_data(physeq.f)$Type)
print(comps)

p1 <- p + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test") 


p1

pg <- p + stat_compare_means(
  comparisons = comps, method = "wilcox.test",
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")))

pg 



#Phyla

mycols <- c("#EB5406", "#4EE2EC", "#d98880", "#a9cce3", "#d2b4de")



p <- plot_taxa_boxplot(physeq,
                       taxonomic.level = "Phylum",
                       top.otu = 25, 
                       group = "Type",
                       add.violin= FALSE,
                       title = "Top twenty five Genus", 
                       keep.other = FALSE,
                       group.order = c("Chicken_gut",	"Droppings",	"Feed",	"Fish_Intestine",	"Sediment"),
                       dot.opacity = .4,
                       box.opacity = .4,
                       group.colors = mycols,
                       dot.size = 2) + theme_biome_utils() + rremove("x.text")

p

physeq.f <- format_to_besthit(physeq)

comps <- make_pairs(sample_data(physeq.f)$Type)
print(comps)

p1 <- p + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test") 


p1

py <- p + stat_compare_means(
  comparisons = comps, method = "wilcox.test",
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")))

py 

#Arrange###################
pA <- py + theme(legend.position = "right")
pA

bB <- pg + theme(legend.position = "right")
bB


all <- ggarrange(pA, bB,
                 labels = c("A", "B"),
                 ncol = 2, nrow = 1)


#############Bar_plot############
###########Phylum_top_20#########################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and verify data
# Check your Phylum names first
head(tax_table(physeq)[,"Phylum"])

# Handle NA values
tax_table(physeq)[is.na(tax_table(physeq)[,"Phylum"]), "Phylum"] <- "Unclassified"

# 2. Merge samples by Type and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Type") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("Phylum")

# 3. Check what phyla exist after merging
print("Phylum abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 10 phyla (or fewer if not available)
top_phyla <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(22, ntaxa(physeq_merged))])

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_phyla)) {
  physeq_top <- prune_taxa(top_phyla, physeq_merged)
  other_phyla <- setdiff(taxa_names(physeq_merged), top_phyla)
  physeq_other <- prune_taxa(other_phyla, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_phyla)
  tax_table(merged_other)[, "Phylum"] <- "Others"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Verify final phyloseq object
print("Final Phylum abundances:")
print(sort(taxa_sums(physeq_final), decreasing = TRUE))

# 7. Prepare plotting data - CRUCIAL STEP
plot_data <- psmelt(physeq_final) %>%
  mutate(
    Phylum = as.character(Phylum),
    Abundance = as.numeric(Abundance)  # Ensure numeric values
  ) %>%
  group_by(Sample, Phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 8. Set color palette
my_colors <- c(
  # Vibrant primary colors (10)
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#B2DF8A", "#A65628", "#F781BF", "#999999", "#66C2A5",
  
  # Distinct secondary colors (10)
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#1B9E77",
  "#E5C494", "#B3B3B3", "#8DD3C7", "#FB8072", "#80B1D3"
)



# Apply alpha (0.8 transparency)
my_colors <- adjustcolor(my_colors, alpha.f = 0.8)

# Use only needed colors
n_phyla <- length(unique(plot_data$Phylum))
phyla_colors <- my_colors[1:n_phyla]
names(phyla_colors) <- unique(plot_data$Phylum)

# 9. Create the plot - SIMPLIFIED VERSION
sp <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_col(position = "fill") +  # Using geom_col instead of geom_bar
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Type", y = "Relative Abundance", fill = "Phylum") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 12, face = "italic"),          # Italic legend labels
    legend.title = element_text(size = 11, face = "plain"),         # Italic legend title
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

sp


#################Genus_top_50#####################################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and prepare Genus data
tax_table(physeq)[is.na(tax_table(physeq)[,"Genus"]), "Genus"] <- "Unclassified"

# 2. Merge samples by Type and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Type") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("Genus")

# 3. Select top 50 Genus
top_n <- 80
top_Genus <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(top_n, ntaxa(physeq_merged))])

# 4. Create "Other" category if needed
if(ntaxa(physeq_merged) > length(top_Genus)) {
  physeq_top <- prune_taxa(top_Genus, physeq_merged)
  other_Genus <- setdiff(taxa_names(physeq_merged), top_Genus)
  physeq_other <- prune_taxa(other_Genus, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_Genus)
  tax_table(merged_other)[, "Genus"] <- "Others"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 5. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    Genus = as.character(Genus),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, Genus) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 6. MANUAL CONTRASTING COLOR PALETTE (51 colors)
Genus_colors <- c(
  # Vibrant primary colors (10)
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  
  # Distinct secondary colors (10)
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#8DD3C7", "#FB8072", "#80B1D3",
  
  # Tertiary colors (10)
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02", "#7570B3",
  
  # Additional high-contrast colors (20)
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928", "#8E0152", "#C51B7D", "#DE77AE",
  "#F1B6DA", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41",
  
  # Final contrasting colors (1)
  "#999999"  # Gray for "Other"
)

# Use only needed colors
n_Genus <- length(unique(plot_data$Genus))
Genus_colors <- Genus_colors[1:n_Genus]
names(Genus_colors) <- unique(plot_data$Genus)

# 7. Create optimized plot
sg <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = Genus_colors) +
  labs(x = "Type", y = "Relative Abundance", fill = "Genus") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 9, face = "italic"),          # Italic legend labels
    legend.title = element_text(size = 10, face = "plain"),         # Italic legend title
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

sg


#Arrange###################
#pA <- py + theme(legend.position = "right")
#pA

#bB <- pg + theme(legend.position = "right")
#bB


all <- ggarrange(sp, sg,
                 labels = c("A", "B"),
                 ncol = 2, nrow = 1)
all


###########################################################
#########################################################

###########################################################
#########################################################
###############Phylum##########################
#All_Metadata_Column
# Load required packages
library(phyloseq)
library(tidyverse)

## Step 1: Convert counts to relative abundance (percentages)
physeq_percent <- transform_sample_counts(physeq, function(x) (x / sum(x)) * 100)
otu_table(physeq_percent) <- round(otu_table(physeq_percent), 5)

## Step 2: Merge taxa at Phylum level
physeq_Phylum <- tax_glom(physeq_percent, "Phylum")

## Step 3: Get all metadata columns (excluding sample names)
metadata_cols <- colnames(sample_data(physeq_Phylum))
metadata_cols <- metadata_cols[!metadata_cols %in% c("Sample")] # Exclude sample ID column

## Step 4: Create a function to process each metadata column
generate_abundance_table <- function(physeq, metadata_column) {
  df <- psmelt(physeq)
  
  result <- df %>%
    group_by(Phylum, !!sym(metadata_column)) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(metadata_column), values_from = Abundance) %>%
    mutate(across(where(is.numeric), ~ round(., 5)))
  
  return(result)
}

## Step 5: Process all metadata columns and save results
for (col in metadata_cols) {
  cat("\nProcessing column:", col, "\n")
  
  # Generate the table
  result <- generate_abundance_table(physeq_Phylum, col)
  
  # Print first few rows
  print(head(result))
  
  # Save to CSV
  safe_col <- make.names(col) # Ensure valid filename
  write.csv(result, 
            file = paste0("Phylum_abundance_by_", safe_col, ".csv"), 
            row.names = FALSE)
  
  cat("Saved to:", paste0("Phylum_abundance_by_", safe_col, ".csv"), "\n")
}

## Optional: Create a combined report with all results
all_results <- lapply(metadata_cols, function(col) {
  res <- generate_abundance_table(physeq_Phylum, col)
  res %>% 
    pivot_longer(-Phylum, names_to = "Category", values_to = "Abundance") %>%
    mutate(Metadata_Column = col) %>%
    select(Metadata_Column, Phylum, Category, Abundance)
}) %>% bind_rows()

write.csv(all_results, "combined_Phylum_abundance_all_metadata.csv", row.names = FALSE)



####################################################
########Genus#################################
#############Genus#########################################
#All_Metadata_Column
# Load required packages
library(phyloseq)
library(tidyverse)

## Step 1: Convert counts to relative abundance (percentages)
physeq_percent <- transform_sample_counts(physeq, function(x) (x / sum(x)) * 100)
otu_table(physeq_percent) <- round(otu_table(physeq_percent), 5)

## Step 2: Merge taxa at Genus level
physeq_Genus <- tax_glom(physeq_percent, "Genus")

## Step 3: Get all metadata columns (excluding sample names)
metadata_cols <- colnames(sample_data(physeq_Genus))
metadata_cols <- metadata_cols[!metadata_cols %in% c("Sample")] # Exclude sample ID column

## Step 4: Create a function to process each metadata column
generate_abundance_table <- function(physeq, metadata_column) {
  df <- psmelt(physeq)
  
  result <- df %>%
    group_by(Genus, !!sym(metadata_column)) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(metadata_column), values_from = Abundance) %>%
    mutate(across(where(is.numeric), ~ round(., 5)))
  
  return(result)
}

## Step 5: Process all metadata columns and save results
for (col in metadata_cols) {
  cat("\nProcessing column:", col, "\n")
  
  # Generate the table
  result <- generate_abundance_table(physeq_Genus, col)
  
  # Print first few rows
  print(head(result))
  
  # Save to CSV
  safe_col <- make.names(col) # Ensure valid filename
  write.csv(result, 
            file = paste0("Genus_abundance_by_", safe_col, ".csv"), 
            row.names = FALSE)
  
  cat("Saved to:", paste0("Genus_abundance_by_", safe_col, ".csv"), "\n")
}

## Optional: Create a combined report with all results
all_results <- lapply(metadata_cols, function(col) {
  res <- generate_abundance_table(physeq_Genus, col)
  res %>% 
    pivot_longer(-Genus, names_to = "Category", values_to = "Abundance") %>%
    mutate(Metadata_Column = col) %>%
    select(Metadata_Column, Genus, Category, Abundance)
}) %>% bind_rows()

write.csv(all_results, "combined_Genus_abundance_all_metadata.csv", row.names = FALSE)


