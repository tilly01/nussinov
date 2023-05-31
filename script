library(gplots)
library(ggplot2)
library(heatmaply)

set.seed(0)

# Define gene names and sample names
genes <- paste0("gene", 1:10)
samples <- c()
for (i in 1:4) {
  for (j in 1:3) {
    sample_name <- paste0("patient_", i, " ", ifelse(j == 1, "_before", "_after"), "_replicate_", j)
    samples <- c(samples, sample_name)
  }
}

# Generate random gene expression data
data <- matrix(rnorm(120, mean = 0, sd = 50), nrow = 10, ncol = 12)
data <- pmax(pmin(data, 100), -100)  # Clip values to -100 and 100

# Create a data frame with the generated data
df <- data.frame(data, row.names = genes, stringsAsFactors = FALSE)

# Calculate the correlation matrix
corr_matrix <- cor(t(df))

# Perform hierarchical clustering
dendro <- hclust(as.dist(1 - corr_matrix), method = "average")

# Interactive heatmap with plotly
heatmaply(df, scale = "none", col = colorRampPalette(c("green", "red"))(256),
          labRow = genes, labCol = samples, xlab = "Samples", ylab = "Genes",
          main = "Gene Expression Heatmap")

# Print the interactive plots
plot(dendro, main = "Hierarchical Clustering of Tissue Samples", xlab = "Samples", ylab = "Distance")

