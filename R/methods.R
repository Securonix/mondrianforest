#' Summary method for mondrianforest.
#'
#' Return a list with a brief summary of the forest and training data.
#'
#' @param mforest a mondrianforest class object.
#' @return A list containing the number of trees, dimensions of the training
#'   data, summary of the data, and parameters.
#' @export

summary.mondrianforest <- function(mforest) {
  out <- list()
  out[["ntree"]] <- length(mforest$forest)
  out[["Dimensions of Data"]] <- list("X_mat"=dim(mforest$data$X_F$X_mat), "F_mat"=dim(mforest$data$X_F$F_mat))
  out[["Data Summary"]] <- purrr::map(mforest$data$X_F, summary)
  out[["Parameters"]] <- mforest$parameters
  out
}
