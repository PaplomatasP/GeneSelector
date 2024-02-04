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
selected_genes <- GeneSelector(GBM_Dataset,Generations=500, min_genes = 200, max_genes = 2500,
                               popSize=400, Labels=GBM_Labels)

# Store the end time
end_time <- Sys.time()

# Compute and print the time difference
time_taken <- end_time - start_time

# Output the selected genes
print(time_taken)
