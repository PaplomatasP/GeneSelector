library(Seurat)
library(GA)  # Genetic algorithm 

# Preprocess before running the pipeline
OptimizeExpressionMatrix <- function(data, min_genes = 200, max_genes = 2500) {
  library(Seurat)
  obj=data
  # Convert data to a sparse matrix if it's a dataframe
  if (is.data.frame(data)) {
    data <- as(as.matrix(data), "dgCMatrix")
  }
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = data)
  
  # Filter cells based on the number of expressed genes
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes )
  
  # Extract the filtered data as a data frame
  seurat_df <- as.data.frame(obj[, rownames(seurat_obj@meta.data)] )
  
  return(seurat_df)
}