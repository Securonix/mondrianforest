.check_types_mf <- function(X, warn = FALSE) {
  if (is.data.frame(X)) {
    if (any(purrr::map_lgl(X, function(x) !is.numeric(x)))) {
      stop("All covariates must be numeric or binary flags. Run model_matrix() on the factors.")
    }
    var_ranges <- purrr::map_dbl(X, function(x) max(x) - min(x))
    multiple <- max(var_ranges)/min(var_ranges)
    if (max(var_ranges)/min(var_ranges) > 3 & warn == TRUE) {
      warning("Maximum covariate range is ", round(multiple, digits = 1), " times minimum range. Consider standardizing covariates.\n" )
    }
    return(as.matrix(X))
  }
  else {
    return(X)
  }
}

.tree_depth_recurser <- function(root, current_depth = 0) {
  node <- root
  if ("leaf" %in% names(node)) {
    cust <- node$cust
    out <- matrix(c(current_depth, cust), nrow = 1)
    if (current_depth == 0) {
      colnames(out) <- c("Leaf_Depth", names(node$cust))
    }
    return(out)
  }
  current_depth <- current_depth + 1
  out <- rbind(tree_depth_recurser(node$R, current_depth),
               tree_depth_recurser(node$L, current_depth))
  colnames(out) <- c("Leaf_Depth", names(node$cust))
  return(out)
}
