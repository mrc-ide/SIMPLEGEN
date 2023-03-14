
#------------------------------------------------
#' @title Load system file
#'
#' @description Load a file from within the inst/extdata folder of the SIMPLEGEN
#'   package. File extension must be one of .csv, .txt, or .rds.
#'
#' @param name the name of a file within the inst/extdata folder.
#'
#' @importFrom utils txtProgressBar
#' @export

test1 <- function(x) {
  
  # make progress bars
  pb_sim <- txtProgressBar(min = 0, max = 100, initial = NA, style = 3)
  args_progress <- list(pb_sim = pb_sim)
  
  # run cpp function
  sim_indiv_deploy(args_progress)
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
  if ((i == max_i) && close) {
    close(pb_list[[name]])
  }
}
