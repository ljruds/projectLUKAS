
#' Get all factors of a positive integer
#'
#' Computes and returns all positive divisors (factors) of a given integer.
#'
#' @param n A positive integer.
#'
#' @return A numeric vector containing all factors of \code{n}, sorted in ascending order.
#'
#' @details
#' This function finds all integers that divide \code{n} without remainder.
#' The implementation checks all integers from 1 to \code{n}.
#' For improved performance with large \code{n}, see \code{get_factors_fast()}.
#'
#' @examples
#' get_factors(12)
#' # Returns: 1 2 3 4 6 12
#'
#' get_factors(7)
#' # Returns: 1 7
#'
#' @export
get_factors <- function(n) {
  if (n <= 0 || n != as.integer(n)) {
    stop("Input must be a positive integer.")
  }

  factors <- which(n %% seq_len(n) == 0)
  return(factors)
}
