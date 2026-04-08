get_factors <- function(n) {
  if (n <= 0 || n != as.integer(n)) {
    stop("Input must be a positive integer.")
  }

  factors <- which(n %% seq_len(n) == 0)
  return(factors)
}
