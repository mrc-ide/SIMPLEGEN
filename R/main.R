
#------------------------------------------------
#' @title text
#'
#' @description text
#'
#' @param x text
#'
#' @importFrom utils txtProgressBar
#' @export

test1 <- function(x) {
  
  # make progress bars
  pb_sim <- txtProgressBar(min = 0, max = 1e8, initial = NA, style = 3)
  args_progress <- list(pb_sim = pb_sim)
  
  # run cpp function
  sim_indiv_deploy(args_progress)
}


