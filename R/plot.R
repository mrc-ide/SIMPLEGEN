
#------------------------------------------------
#' @title Plot daily counts of each host state
#'
#' @description TODO.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param deme which deme to plot.
#' @param states which states to plot. Can be any subset of \code{c("Sh", "Eh",
#'   "Ah", "Ch", "Sv", "Ev", "Iv")}.
#'
#' @export

plot_daily_states <- function(project, deme = 1, states = c("Sh", "Eh", "Ah", "Ch")) {
  
  # check inputs
  assert_custom_class(project, "simplegen_project")
  
  # check for epi output
  assert_dataframe(project$epi_output$daily_values)
  
  # check deme within output
  assert_in(deme, project$epi_output$daily_values$deme)
  
  # subset to desired rows and columns
  df_wide <- project$epi_output$daily_values[, c("time", "deme", states)]
  df_wide <- df_wide[df_wide$deme == deme,]
  
  # get to long format
  state <- count <- NULL  # needed to avoid error "no visible global variable"
  df_long <- df_wide
  df_long <- tidyr::gather(df_wide, state, count, states, factor_key = TRUE)
  
  # produce plot
  ggplot2::ggplot(df_long) + ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes_(x = ~time, y = ~count, color = ~state))
  
}
