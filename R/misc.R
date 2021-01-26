
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
# does interval x intersect interval y
# both integer vectors of length 2
#' @noRd
interval_intersect <- function(x, y) {
  (x[1] < y[2]) & (x[2] > y[1])
}
