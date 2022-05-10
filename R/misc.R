
# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
matrix_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix
#' @noRd
rcpp_to_matrix <- function(x) {
  ret <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
  return(ret)
}

# -----------------------------------
# takes list format returned from Rcpp and converts to three-dimensional array.
# Array indexing is in the same order as the underlying list, for example
# x[i,j,k] is equivalent to l[[i]][[j]][[k]]
#' @noRd
rcpp_to_array <- function(x) {
  ret <- array(unlist(x), dim = c(length(x[[1]][[1]]), length(x[[1]]), length(x)))
  ret <- aperm(ret, perm = c(3,2,1))
  return(ret)
}

#------------------------------------------------
# update progress bar
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
# close = whether to close when reach end
#' @importFrom utils setTxtProgressBar
#' @noRd
update_progress <- function(pb_list, name, i, max_i, close = TRUE) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i & close) {
    close(pb_list[[name]])
  }
}

#------------------------------------------------
# produce all unique pairwise comparisons from 1 to n in long dataframe format
#' @noRd
pairwise_long <- function(n, include_diagonal = FALSE) {
  if (include_diagonal) {
    ret <- data.frame(x = rep(seq_len(n), times = rev(seq_len(n))),
                      y = unlist(mapply(function(x) x:n, 1:n, SIMPLIFY = FALSE)))
  } else {
    ret <- data.frame(x = rep(seq_len(n - 1), times = rev(seq_len(n - 1))),
                      y = unlist(mapply(function(x) x:n, 2:n, SIMPLIFY = FALSE)))
  }
  return(ret)
}


