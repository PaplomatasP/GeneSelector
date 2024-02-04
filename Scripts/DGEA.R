# Function to run the Differential Gene Expression Analysis (DGEA)
wilcoxFun <- function(data, Labels) {
  
  data <- t(data)
  
  # Prepare the metadata
  MetaData <- data.frame(groups = Labels )
  
  rownames(MetaData) <- colnames(data)  # Match row and column names
  
  # Create a Seurat object with metadata
  suppressWarnings(Input <- Seurat::CreateSeuratObject(counts = data, meta.data = MetaData) )
  suppressMessages(Input <- NormalizeData (Input,verbose = FALSE) )
  
  # Determine the number of levels in Labels
  num_levels <- length(levels(as.factor(Labels)))
  
  # Create a list of ident values based on the number of levels
  ident_values <- levels(as.factor(Labels))[1:num_levels]
  
  # Create a list of ident arguments for Seurat::FindMarkers
  ident_args <- paste0("ident.", 1:num_levels, " = ident_values[", 1:num_levels, "]")
  
  # Combine the ident arguments into a single string
  ident_string <- paste(ident_args, collapse = ", ")
  
  # Use the generated ident_string in Seurat::FindMarkers
  result <- eval(parse(text = paste("Seurat::FindMarkers(object = Input, 
                                    group.by = 'groups', logfc.threshold = -Inf, 
                                    test.use = 'wilcox', only.pos = FALSE, verbose = FALSE,
                                    min.cells.group = 1, ", ident_string, ")")))
  
  # Convert results to a data frame for easier manipulation
  results_Seurat <- data.frame(
    gene_names = rownames(result),
    pvalue = result$p_val,
    FDR = result$p_val_adj,
    logFC = result$avg_log2FC
  )
  
  # Filter for significant and differentially expressed genes
  significant_Statistical_genes <- subset(
    results_Seurat, 
    FDR <= 0.05 & abs(logFC) >= 1
  )
  
  # Get the names of the significant genes
  Sign_Genes <- significant_Statistical_genes
  
  return(results_Seurat)
}