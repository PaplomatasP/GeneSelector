# Custom function for initializing the population of the Genetic Algorithm (GA)
custom_initialization <- function(object, ...) {
  nBits <- object@nBits  # Number of genes
  popSize <- object@popSize  # Population size
  prob <- 0.01  # 1% chance of selecting each feature
  
  # Initialize the population with each gene having a 'prob' chance of being selected
  population <- matrix(rbinom(nBits * popSize, size = 1, prob = prob), nrow = popSize, ncol = nBits)
  return(population)
}

# Genetic_similarity is a function that performs genetic algorithm-based feature selection
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

Genetic_similarity <- function(data, Labels, popSize, Generations) {
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
    
    # Calculate and normalize the p-values, logFC, and variance
    pvalues_adjusted <- -log10(wilcox_results$pvalue + 1e-100)
    
    transformed_pvalues <- 1 - (-log10(wilcox_results$pvalue + 1e-100))
    pvalue_mean <- mean(normalize_score(transformed_pvalues, min(transformed_pvalues), max(transformed_pvalues)))
    
    logFC_mean <- mean(normalize_score(abs(wilcox_results$logFC), min(abs(wilcox_results$logFC)), max(abs(wilcox_results$logFC))))
    
    variance_values <- colSums((selected_genes - colMeans(selected_genes))^2) / (nrow(selected_genes) - 1)
    variance_mean <- mean(normalize_score(variance_values, min(variance_values), max(variance_values)))
    
    fitness <- pvalue_mean *0.2 + logFC_mean * 0.3  + variance_mean * 0.5
    
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
    pcrossover = runif(1, min = 0.6, max = 0.95),
    pmutation = runif(1, min = 0.001, max = 0.1),
    elitism = base::max(1, round(popSize * 0.5)),
    population = custom_initialization
  )
  
  # Results 
  selected_gene_indices <- which(ga_result@solution[1, ] == 1)
  top_genes <- colnames(data)[selected_gene_indices]
  
  return(top_genes)
}