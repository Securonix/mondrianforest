compute_G <- function(cust, tab, G_parent, d) { # helper function
  if (sum(cust) > 0) {
    return((cust - d*tab + d*sum(tab)*G_parent)/sum(cust))
  }
  if (sum(cust) == 0) {
    return(G_parent)
  }
}

#' @importFrom purrr map2_lgl
predict_mondrian_obs <- function(tree, x, f, gamma, f_scale) {
  H <- vector(mode = "numeric", length = length(tree$cust)) # Zero vector with correct length,
  H[sample(seq_along(H), size = 1)] <- 1 # Sample label uniformly.

  if (!is.null(tree$leaf)) { # If root is leaf.
    return((tree$cust+H)/sum(tree$cust + H))
  }

  G_parent <- H # Base distribution
  s <- 0
  p <- 1
  parent_node <- tree
  current_node <- if ((parent_node$split$split_type=="categorical" & f[parent_node$split$split_dimension]!=parent_node$split$split_point) |
                      (parent_node$split$split_type=="numeric" & x[parent_node$split$split_dimension]<parent_node$split$split_point)) {parent_node$L} else {parent_node$R}

  while (TRUE) {
    time_diff <- current_node$split_time - parent_node$split_time
    eta_x <- sum(pmax.int(x - current_node$B_j$upper, 0)) + sum(pmax.int(current_node$B_j$lower - x, 0))
    eta_f <- sum(map2_lgl(current_node$B_j$factor_uniques, f, function(x, y) {
      return(!(y %in% x))
    })*f_scale)
    eta <- eta_x + eta_f
    p_x <- 1 - exp(-time_diff*eta)
    if (p_x > 0) { # x branches off into two new nodes j_tilde => child(j_tilde)
      d_bar <- (1 + exp(-eta*time_diff))/2 ## Expected discount
      cust <- tab <- pmin(current_node$cust, 1L) # counts for new parent node j_tilde
      G <- compute_G(cust, tab, G_parent, d_bar) # G will be same for child(j_tilde) as j_tilde
      s <- s + p*p_x*G
    }
    if (!is.null(current_node$leaf)) {
      d <- exp(-gamma*time_diff)
      s <- s + p*(1 - p_x)*compute_G(current_node$cust, current_node$tab, G_parent, d)
      return(s)
    }
    else {
      d <- exp(-gamma*time_diff)
      p <- p*(1 - p_x)
      G_parent <- compute_G(current_node$cust, current_node$tab, G_parent, d) # Assign current G to G_parent
      parent_node <- current_node
      if ((current_node$split$split_type=="categorical" & f[current_node$split$split_dimension]!=current_node$split$split_point) |
          (current_node$split$split_type=="numeric" & x[current_node$split$split_dimension]<current_node$split$split_point)) {
        current_node <- current_node$L
      }
      else {
        current_node <- current_node$R
      }
    }
  }
}

#' Predicts with an existing mondrianforest
#'
#' Predict method for object of class 'mondrianforest'.
#'
#' @param mforest A mondrianforest.
#' @param newdata A new dataset. If missing, predicts on training set.
#' @param gamma Hierarchical smoothing parameter, defaults to 0.0001.
#' @return A matrix with class probabilities if type = "prob". If type = "class", returns MAP class estimates.
#' @examples
#' predict(mf, type = "class")
#'
#' @export
#' @importFrom stats predict

predict.mondrianforest <- function(mforest, newdata, gamma = 0.0001, type = "prob") {
  if (missing(newdata)) {
    newdata_list <- mforest$data$X_F
  }
  else {
    if (!any(c("data.frame", "matrix") %in% class(newdata))) {
      stop("newdata must be of class 'data.frame' or 'matrix'.")
    }
    if (names(mforest$parameters$y_col_num) %in% names(newdata)) {
      newdata <- newdata[, -mforest$parameters$y_col_num] # get rid of label if it's there.
    }
    newdata_list <- convert_to_matrices(newdata)
  }
  if (!inherits(mforest, "mondrianforest")) {
    stop("Object not of class mf.")
  }
  if (!(type %in% c("prob", "class"))) {
    stop("Argument 'type' must be either 'prob' or 'class'.")
  }
  num_levels <- length(levels(mforest$data[["Y"]]))

  preds <- matrix(nrow = max(nrow(newdata_list$X_mat), nrow(newdata_list$F_mat)), ncol = num_levels)
  for (i in 1:nrow(preds)) {
    p <- matrix(nrow = length(mforest$forest), ncol = num_levels)
    for (j in 1:length(mforest$forest)) {
      p[j, ] <- predict_mondrian_obs(tree = mforest$forest[[j]], x = newdata_list$X_mat[i, , drop = TRUE], f = newdata_list$F_mat[i, , drop = TRUE], gamma = gamma, f_scale = mforest$parameters[["f_scale"]]) # For binary predictions, only need one element
    }
    preds[i, ] <- colMeans(p)
  }
  colnames(preds) <- levels(mforest$data[["Y"]])
  if (type == "class") {
    preds <- apply(preds, 1, function(x) as.character(which.max(x)))
    for (ind in 1:length(levels(mforest$data[["Y"]]))) {
      preds[preds==as.character(ind)] <- levels(mforest$data[["Y"]])[ind]
    }
    return(preds)
  }
  if (type == "prob") {
    return(preds)
  }
}
