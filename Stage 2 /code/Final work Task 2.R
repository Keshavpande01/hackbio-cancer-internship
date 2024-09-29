# Load necessary libraries
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# Load the dataset
genes_data <- read.csv('https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv', row.names = 1)

# Convert the data frame to a matrix
genes_matrix <- as.matrix(genes_data)

# Sequential color palette
seq_pal <- colorRampPalette(brewer.pal(n=9, name="Blues"))(100)

# Generate Sequential heatmap
heatmap.2(genes_matrix, col=seq_pal, scale="row", trace="none", main="Sequential Heatmap")

# Diverging color palette
div_pal <- bluered(100)

# Generate Diverging heatmap
heatmap.2(genes_matrix, col=div_pal, scale="row", trace="none", main="Diverging Heatmap")

# Clustering heatmaps
# Row clustering only
heatmap.2(genes_matrix, Rowv=TRUE, Colv=FALSE, dendrogram="row", col=div_pal, scale="row", trace="none", main="Gene Clustered Heatmap")

# Column clustering only
heatmap.2(genes_matrix, Rowv=FALSE, Colv=TRUE, dendrogram="col", col=div_pal, scale="row", trace="none", main="Sample Clustered Heatmap")

# Both row and column clustering
heatmap.2(genes_matrix, Rowv=TRUE, Colv=TRUE, dendrogram="both", col=div_pal, scale="row", trace="none", main="Both Gene and Sample Clustered Heatmap")

# 5. Distribution of Fold Change
group_1 <- genes_matrix[, 1:5]  # Select the first 5 columns as Group 1
group_2 <- genes_matrix[, 6:10]  # Select the next 5 columns as Group 2

group_1_mean <- apply(group_1, 1, mean)
group_2_mean <- apply(group_2, 1, mean)

fold_change <- log2(group_1_mean) - log2(group_2_mean)
hist(fold_change, xlab = "log2 Fold Change (Group 1 vs Group 2)",
     main = "Distribution of Fold Change")

# 6. Volcano Plot and Lollipop Chart

# Calculate p-values for each gene
p_values <- numeric(nrow(genes_matrix))
for (i in 1:nrow(genes_matrix)) {
  g1 <- as.numeric(group_1[i, ])
  g2 <- as.numeric(group_2[i, ])
  t_test_result <- t.test(g1, g2, var.equal = FALSE)
  p_values[i] <- t_test_result$p.value
}

# Adjust p-values using the Benjamini-Hochberg method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Set thresholds
p_value_threshold <- 0.05
fold_change_threshold <- 1  # Corresponds to 2-fold change

# Identify significant genes
significant_genes <- data.frame(
  Gene = rownames(genes_matrix),
  log2FoldChange = fold_change,
  pValue = p_values
)

# Add a column for significance
significant_genes$Significant <- "Not Significant"
significant_genes$Significant[significant_genes$pValue < p_value_threshold & 
                                significant_genes$log2FoldChange > fold_change_threshold] <- "Upregulated"
significant_genes$Significant[significant_genes$pValue < p_value_threshold & 
                                significant_genes$log2FoldChange < -fold_change_threshold] <- "Downregulated"

# Volcano Plot
ggplot(significant_genes, aes(x = log2FoldChange, y = -log10(pValue), color = Significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 P-Value") +
  geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-fold_change_threshold, fold_change_threshold), linetype = "dashed", color = "black")

ggsave("volcano_plot.png", width = 10, height = 6)

# Example lollipop plot for top pathways
top_pathways <- data.frame(Pathway=c("Proteolysis", "Glutathione Metabolism", "Cell Cycle"), nGenes=c(10, 8, 6), pValue=c(0.001, 0.005, 0.01))

ggplot(top_pathways, aes(x=reorder(Pathway, -nGenes), y=nGenes)) +
  geom_point(aes(size=-log10(pValue), color=-log10(pValue))) +
  coord_flip() +
  labs(title="Top Pathways", x="Pathway", y="Gene Count", size="-log10(p-value)")



