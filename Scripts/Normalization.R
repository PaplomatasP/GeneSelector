# Function to normalize scores
normalize_score <- function(score, score_min, score_max) {
  # Check if the range of values is zero
  if (score_max == score_min) {
    return(rep(0, length(score)))  # Returns a vector with zero values
  }
  
  normalized_score <- (score - score_min) / (score_max - score_min)
  normalized_score[is.nan(normalized_score)] <- 0  # Handle zero divide
  return(normalized_score)
}