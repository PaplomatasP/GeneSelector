library(Seurat)
library(GA)  # Genetic algorithm 

# Preprocess before run the pipeline
OptimizeExpressionMatrix <- function(data, min_genes = 200, max_genes = 2500) {
  library(Seurat)
  obj=data
  # Convert data to a sparse matrix if it's a dataframe
  if (is.data.frame(data)) {
    data <- as(as.matrix(data), "dgCMatrix")
  }
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = data)
  
  
  # Filter cells
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes )
  
  seurat_df=as.data.frame(obj[, rownames(seurat_obj@meta.data)] )
  
  return(seurat_df)
}



# Function run the DGEA analysis

wilcoxFun <- function(data, Labels) {
  
  data <- t(data)
  
  # Prepare the metadata
  MetaData <- data.frame(groups = Labels )
  
  rownames(MetaData) <- colnames(data)  # Match row and column names
  
  # Create a Seurat object
  suppressWarnings(Input <- Seurat::CreateSeuratObject(counts = data, meta.data = MetaData) )
  suppressMessages(Input <- NormalizeData (Input,verbose = FALSE) )
  #Determine the number of levels in Labels
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

# Function make the normalization

normalize_score <- function(score, score_min, score_max) {
  # # Check if the range of values is zero
  if (score_max == score_min) {
    return(rep(0, length(score)))  # Returns a vector with zero values
  }
  
  normalized_score <- (score - score_min) / (score_max - score_min)
  normalized_score[is.nan(normalized_score)] <- 0  # Zero divide solution
  return(normalized_score)
}

# Costum function for initialization the population of the GA
custom_initialization <- function(object, ...) {
  nBits <- object@nBits  # Number of genes
  popSize <- object@popSize  # Population size
  prob <- 0.01  # 1% chance of selecting each feature
  
  # Initialize population with each gene having a 'prob' chance of being selected
  population <- matrix(rbinom(nBits * popSize, size = 1, prob = prob), nrow = popSize, ncol = nBits)
  return(population)
}


# # Genetic_similarity is a function that performs genetic algorithm-based feature selection
# on single-cell RNA sequencing (scRNA-seq) data to identify genes that exhibit relevant
# biological differences between groups defined by 'Labels'.
# 
# Input:
# - data: A matrix or dataframe containing scRNA-seq data.
# - Labels: A vector specifying group labels or conditions for the data.
# - popSize: Population size for the genetic algorithm.
# - Generations: The number of generations (iterations) for the genetic algorithm.
# 
# The function first normalizes the input data and then evaluates the fitness of gene subsets.
# Fitness is determined based on a combination of statistical measures including Wilcoxon test results
# (p-values and log fold change), gene variance, and normalization scores.
# 
# The genetic algorithm seeks to maximize fitness, ultimately selecting a subset of genes that
# are likely to be biologically relevant for distinguishing between the defined groups.
# 
# The function returns a list of selected genes and their fitness scores.


Genetic_similarity <- function(data, Labels,popSize,Generations) {
  set.seed(999)
  Input <- Seurat::CreateSeuratObject(counts = data)
  invisible(capture.output(Input <- NormalizeData(Input), file = NULL))
  
  data <- GetAssayData(Input, assay = "RNA", layer = "data")
  data <- as.data.frame(data)
  
  
  fitness_function <- function(gene_indices) {
    set.seed(999)
    selected_genes <- data[, gene_indices == 1, drop = FALSE]
    if (ncol(selected_genes) == 0) {
      return(0)
    }
    
    # Wilcoxon Test 
    wilcox_results <- suppressMessages ( wilcoxFun(selected_genes, Labels) )
    
    # Calculate and normalization the  p-values, logFC and variance
    pvalues_adjusted <- -log10(wilcox_results$pvalue + 1e-100)
    
    transformed_pvalues <- 1 - (-log10(wilcox_results$pvalue + 1e-100))
    pvalue_mean <- mean(normalize_score(transformed_pvalues, min(transformed_pvalues), max(transformed_pvalues)))
    
    logFC_mean <- mean(normalize_score(abs(wilcox_results$logFC), min(abs(wilcox_results$logFC)), max(abs(wilcox_results$logFC))))
    
    variance_values <- colSums((selected_genes - colMeans(selected_genes))^2) / (nrow(selected_genes) - 1)
    variance_mean <- mean(normalize_score(variance_values, min(variance_values), max(variance_values)))
    
    
    fitness <- pvalue_mean *0.2 + logFC_mean * 0.3  + variance_mean * 0.5 # * 0.4
    
    
    return(fitness)
  }
 

  # RUN GA
  ga_result <- ga(
    type = "binary",
    fitness = fitness_function,
    nBits = ncol(data), 
    popSize = popSize,
    maxiter = Generations,
    selection = ga_tourSelection,
    pcrossover = runif(1, min = 0.6, max = 0.95) ,#, 0.9,
    pmutation = runif(1, min = 0.001, max = 0.1) , #0.03,
    elitism = base::max(1, round(popSize * 0.5)),
    population = custom_initialization
  )
  
  
  # Results 
  selected_gene_indices <- which(ga_result@solution[1, ] == 1)
  top_genes <- colnames(data)[selected_gene_indices]
  
  return(top_genes)
}

# GeneSelector is a function for identifying and selecting genes with high biological significance
# from single-cell RNA sequencing (scRNA-seq) data using a genetic algorithm-based approach.
# 
# Input:
# - data: A matrix or dataframe containing scRNA-seq data.
# - Generations: The number of generations (iterations) for the genetic algorithm.
# - min_genes: The minimum number of genes to consider in the selection process.
# - max_genes: The maximum number of genes to consider in the selection process.
# - popSize: Population size for the genetic algorithm.
# - Labels: A vector specifying group labels or conditions for the data.
# 
# The function first optimizes the input data by adjusting the number of genes to fall within the
# specified range (min_genes to max_genes). It then employs a genetic algorithm to evaluate the
# fitness of gene subsets, aiming to maximize biological relevance for distinguishing between the
# defined groups represented by 'Labels'.
# 
# The selected genes are returned as 'G_S_genes', and the filtered scRNA-seq data containing
# only these genes is returned as 'data_GA_filtered'.

GeneSelector <- function(data, Generations,min_genes = 200, max_genes = 2500,
                        popSize,Labels) {
  set.seed(999) 
  pb <- txtProgressBar(min = 0, max = 4, style = 3)  # create a progress bar with 5 steps
  
  data <- OptimizeExpressionMatrix(data=data, min_genes = min_genes, max_genes = max_genes) 
  setTxtProgressBar(pb, 2)
  
  # Step 1: Genetic Algorithm (GA)
  
  G_S_genes <- Genetic_similarity(data,Labels,popSize=popSize,Generations=Generations)  # Assuming Genetic_similarity returns top n_GA genes
  setTxtProgressBar(pb, 3)
  # Reduce data to only include G_S_genes
  data_GA_filtered <- as.data.frame(data[, G_S_genes] )
  
  close(pb)  # Close the progress bar
  
  return( list(G_S_genes=G_S_genes ,  data_GA_filtered= data_GA_filtered,Labels=Labels))
}


#### EXAMPLE
# Store the start time
start_time <- Sys.time()

# Use the GeneSelector function on your example data
# Adjust n_GA, n_PCA, and n_NMF based on how many genes you want to keep at each step
selected_genes <- GeneSelector(GBM_Dataset,Generations=400, min_genes = 200, max_genes = 2500,
                               popSize=400, Labels=GBM_Labels)

# Store the end time
end_time <- Sys.time()

# Compute and print the time difference
time_taken <- end_time - start_time

# Output the selected genes
print(time_taken)

