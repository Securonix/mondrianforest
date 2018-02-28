#' Extends a MondrianForest.
#'
#' \code{mondrian_forest} implements a fast (looks at multiple new observations
#' at once) version of Lakshminarayanan et al's Extend Mondrian Block algorithm
#' described in [Lakshminarayanan et al. 2014]
#' (\url{https://arxiv.org/abs/1406.2673}) with a modification allowing for more
#' space-efficient treatment of categorical variables as if they were
#' dummy-encoded.
#'
#' @param mforest A mondrianforest.
#' @param newdata Data (matrix or data frame) containing new data to extend the
#'   forest.
#' @return A mondrianforest.
#' @examples
#' library(mondrianforest); library(dplyr); library(purrr); library(magrittr)
#' set.seed(1)
#' test <- data.frame(x1 = rnorm(1000),
#'                    x2 = runif(1000),
#'                    x3 = rbeta(n = 1000, shape1 = 3, shape2 = 8),
#'                     # x3 is noise
#'                    x4 = rbinom(n = 1000, size = 4, prob = 0.2)) %>%
#'   map_df(function(x) (x - min(x))/(max(x) - min(x))) %>%
#'   mutate(y = x1*x2^2 + exp(x1) - x4,
#'          x1 = cut_number(sin(x1), n = 10),
#'          label = as.factor(case_when(y < 1 ~ "A",
#'                                      y >=1 & y < 1.7 ~ "B",
#'                                      y >= 1.7 ~ "C",
#'          ))) %>%
#'   select(-y)
#' mf <- mondrian_forest(test[1:5, ], y_col_num = 5, lambda = 3)
#' table(test$label[751:1000], predict(mf, test[751:1000, ], type = "class"))
#' mf_ext <- extend_mondrian_forest(mf, test[6:750, ])
#' table(test$label[751:1000], predict(mf_ext, test[751:1000, ], type = "class"))
#' @export

extend_mondrian_forest <- function(mforest, newdata) {
  # Need to update mforest$data$X and Y
  # Convert to X_F_all.
  new_X_F <- convert_to_matrices(newdata[, -mforest$parameters$y_col_num, drop = FALSE]) # Turn to matrix
  X_F_all <- purrr::map2(mforest$data$X_F, new_X_F, function(old_mat, new) {
    all <- rbind(old_mat, new)
    rownames(all) <- 1:nrow(all)
    return(all)
  })

  mforest$data$X_F <- X_F_all
  new_Y <- as.factor(newdata[[mforest$parameters$y_col_num]])
  all_Y <- factor(c(as.character(mforest$data$Y), as.character(new_Y)), levels = unique(c(levels(new_Y), levels(mforest$data$Y))))
  mforest$data$Y <- all_Y

  # X_mat, F_mat, Y are data used to build the tree. x,f,y are new data used to extend.
  for (j in 1:(length(mforest$forest))) {
    mforest$forest[[j]] <- extend_mondrian_tree(lambda = mforest$parameters[["lambda"]],
                                                 f_scale = mforest$parameters[["f_scale"]],
                                                 tree = mforest$forest[[j]],
                                                 X_mat = X_F_all$X_mat[1:(length(all_Y) - nrow(newdata)), , drop = FALSE],
                                                 F_mat = X_F_all$F_mat[1:(length(all_Y) - nrow(newdata)), , drop = FALSE],
                                                 Y = all_Y[1:(length(all_Y) - nrow(newdata))],
                                                 x = X_F_all$X_mat[((length(all_Y) - nrow(newdata) + 1):length(all_Y)), , drop = FALSE],
                                                 f = X_F_all$F_mat[((length(all_Y) - nrow(newdata) + 1):length(all_Y)), , drop = FALSE],
                                                 y = all_Y[((length(all_Y) - nrow(newdata) + 1):length(all_Y))])
  }
  class(mforest) <- "mondrianforest"
  return(mforest)
}

extend_mondrian_tree <- function(lambda, f_scale, tree, X_mat, F_mat, Y, x, f, y) {
  tree <- extend_mondrian_block(lambda, f_scale, tree, X_mat, F_mat, Y, x, f, y, split_time_parent = 0)
  return(tree)
}

#' @importFrom purrr map_int
extend_mondrian_block <- function(lambda, f_scale, node, X_mat, F_mat, Y, x, f, y, split_time_parent) {
  if ("pause" %in% names(node)) {
    node <- update_paused_node(lambda, f_scale, node, X_mat, F_mat, Y, x, f, y, split_time_parent)
    return(node)
  }
  else {
    e_j <- compute_extensions(node$B_j$lower, node$B_j$upper, node$B_j$factor_uniques, x, f) # Create vector of extensions. Make this faster and not a matrix.
    e_factor_lengths <- map_int(e_j$e_factor_uniques, length)
    E <- rexp_or_lambda(sum(e_j$e_lower + e_j$e_upper) + sum(e_factor_lengths*f_scale), lambda)

    if (E + split_time_parent < node$split_time) {
      split <- sample_ext_dimension_split_pt(node$B_j$upper, node$B_j$lower, e_j$e_upper, e_j$e_lower, e_j$e_factor_uniques, e_factor_lengths*f_scale)
      node <- replace_by_parent(lambda, f_scale, node, X_mat, F_mat, Y, x, f, y, split, E + split_time_parent, e_j) # Adds parent above and recurses.
    }
    else { # Failed to Split. Extend Box, update counts if leaf, else recurse down.
      node <- update_node_extent_parallel(node, node$B_j$lower, e_j$e_lower, node$B_j$upper, e_j$e_upper, node$B_j$factor_uniques, e_j$e_factor_uniques)
      if ("leaf" %in% names(node)) {
        node <- update_counts(node, y)
      }
      else {
        node <- recurse_down_tree(lambda, f_scale, node, X_mat, F_mat, Y, x, f, y)
      }
    }
    return(node)
  }
}


# Returns unpaused node via sampling mondrian block
update_paused_node <- function(lambda, f_scale, node, X_mat, F_mat, Y, x, f, y, split_time_parent) {
  X_j <- X_mat[node$obs, , drop = FALSE]
  F_j <- F_mat[node$obs, , drop = FALSE]
  Y_j <- Y[node$obs]
  if (length(unique(Y_j))==1 && length(unique(y))==1 && Y_j[1]==y[1]) { # pausing conditions for new obs
    node[["obs"]] <- c(node$obs, as.integer(rownames(x)))
    node <- update_counts(node, y)
    for (i in seq_len(ncol(x))) {
      # tryCatch(node[["B_j"]]$lower <- min(x[, i, drop = TRUE], node$B_j$lower[i]),error=stop(), warning=broser())
      node[["B_j"]]$lower[i] <- min(x[, i, drop = TRUE], node$B_j$lower[i]) # min(vector, scalar) = min(min(vector), scalar)
      node[["B_j"]]$upper[i] <- max(x[, i, drop = TRUE], node$B_j$upper[i])
    }
    for (i in seq_len(ncol(f))) {
      node[["B_j"]]$factor_uniques[[i]] <- unique(c(node[["B_j"]]$factor_uniques[[i]], f[, i]))
    }
  }
  else {
    node <- sample_mondrian_block(lambda, f_scale, X_j = rbind(X_j, x), F_j = rbind(F_j, f), Y_j = factor(c(as.character(Y_j), as.character(y)), levels = unique(c(levels(Y_j), levels(y)))), split_time_parent = split_time_parent) # unpause.
  }
  return(node)
}

disjoint_runif <- function(lower1, upper1, lower2, upper2) {
  length1 <- upper1 - lower1
  length2 <- upper2 - lower2
  interval <- sample(c("length1", "length2"), size = 1, prob = c(length1, length2))
  X <- switch(interval, "length1" = runif(1, lower1, upper1),
              "length2" = runif(1, lower2, upper2))
  return(X)
}

#' @importFrom purrr map2
compute_current_extent <- function(B_j, e_j) {
  # Updates current extent through
  upper <- B_j$upper + e_j$e_upper
  lower <- B_j$lower - e_j$e_lower
  factor_uniques <- map2(B_j$factor_uniques, e_j$e_factor_uniques, c)
  return(list("lower" = lower,
              "upper" = upper,
              "factor_uniques" = factor_uniques))
}

compute_extensions <- function(lower, upper, factor_uniques, x, f) {
  e_lower <- e_upper <- vector(length = ncol(x))
  e_factor_uniques <- list()
  for (i in seq_len(ncol(f))) {
    e_factor_uniques[[i]] <- setdiff(unique(f[, i]), factor_uniques[[i]])
  }
  e_lower <- vector()
  for (i in seq_len(ncol(x))) {
    e_lower[i] <- max(lower[i] - min(x[, i]), 0)
    e_upper[i] <- max(max(x[, i]) - upper[i], 0)
  }
  return(list("e_lower"=e_lower, "e_upper"=e_upper, "e_factor_uniques"=e_factor_uniques))
}

sample_ext_dimension_split_pt <- function(upper, lower, e_upper, e_lower, e_factor_uniques, e_factor_scales) { # Make sample dimension and split point
  variable_type <- sample(x = c("numeric", "categorical"), size = 1, prob = c(sum(e_upper+e_lower), sum(e_factor_scales)))
  if (variable_type=="numeric") {
    split_dimension <- sample(seq_along(e_upper), size = 1, prob = e_upper+e_lower)
    split_point <- disjoint_runif(lower[split_dimension] - e_lower[split_dimension],
                                   lower[split_dimension],
                                   upper[split_dimension],
                                   upper[split_dimension] + e_upper[split_dimension])
    return(list("split_dimension"=split_dimension, "split_point"=split_point, "split_type"=variable_type))
  }
  else {
    split_dimension <- sample(seq_along(e_factor_uniques), size = 1, prob = e_factor_scales) # Expects a list of the extensions.
    split_category <- sample(e_factor_uniques[[split_dimension]], size = 1)
    return(list("split_dimension"=split_dimension, "split_point"=split_category, "split_type"=variable_type))
  }
}

update_counts <- function(node, Y_j) {
  if (!missing(Y_j)) {
    node[["cust"]] <- node[["cust"]] + get_vector_table(Y_j)              # Cust and Tab are needed for predictions.
    node[["tab"]] <- named_pmin(node$cust, 1L)
    return(node)
  }
  else {
    node[["cust"]] <- node$L$tab + node$R$tab
    node[["tab"]] <- named_pmin(node$cust, 1L)
    return(node)
  }
}

replace_by_parent <- function(lambda, f_scale, node, X_mat, F_mat, Y, x, f, y, split, split_time, e_j) {
  if (split$split_type == "categorical") {
    outside_obs <- split$split_point == f[, split$split_dimension, drop = TRUE]
  }
  else if (split$split_point < node$B_j$lower[split$split_dimension]) {
    outside_obs <- x[, split$split_dimension, drop = TRUE] < split$split_point
  }
  else if (split$split_point >= node$B_j$upper[split$split_dimension]) {
    outside_obs <- x[, split$split_dimension, drop = TRUE] >= split$split_point
  }

  parent <- list("B_j" = compute_current_extent(node$B_j, e_j),
                 "split_time" = split_time,
                 "split" = split)
  if (split$split_type == "numeric" && split$split_point < node$B_j$lower[split$split_dimension]) {
    parent[["L"]] <- sample_mondrian_block(lambda,
                                            f_scale,
                                            X_j = x[outside_obs, , drop = FALSE],
                                            F_j = f[outside_obs, , drop = FALSE],
                                            Y_j = y[outside_obs],
                                            split_time)
    parent[["R"]] <- if (sum(!outside_obs)==0) node else extend_mondrian_block(lambda,
                                                                                f_scale,
                                                                                node,
                                                                                X_mat,
                                                                                F_mat,
                                                                                Y,
                                                                                x[!outside_obs, , drop = FALSE],
                                                                                f[!outside_obs, , drop = FALSE],
                                                                                y[!outside_obs],
                                                                                split_time)
  }
  if (split$split_type == "numeric" && split$split_point >= node$B_j$upper[split$split_dimension]) {
    parent[["R"]] <- sample_mondrian_block(lambda = lambda,
                                            f_scale = f_scale,
                                            X_j = x[outside_obs, , drop = FALSE],
                                            F_j = f[outside_obs, , drop = FALSE],
                                            Y_j = y[outside_obs],
                                            split_time_parent = split_time)
    parent[["L"]] <- if (sum(!outside_obs)==0) node else extend_mondrian_block(lambda,
                                                                                f_scale,
                                                                                node,
                                                                                X_mat,
                                                                                F_mat,
                                                                                Y,
                                                                                x[!outside_obs, , drop = FALSE],
                                                                                f[!outside_obs, , drop = FALSE],
                                                                                y[!outside_obs],
                                                                                split_time)
  }
  if (split$split_type == "categorical") {
    parent[["L"]] <- if (sum(!outside_obs)==0) node else extend_mondrian_block(lambda,
                                                                                f_scale,
                                                                                node,
                                                                                X_mat,
                                                                                F_mat,
                                                                                Y,
                                                                                x[!outside_obs, , drop = FALSE],
                                                                                f[!outside_obs, , drop = FALSE],
                                                                                y[!outside_obs],
                                                                                split_time)
    parent[["R"]] <- sample_mondrian_block(lambda,
                                            f_scale,
                                            x[outside_obs, , drop = FALSE],
                                            f[outside_obs, , drop = FALSE],
                                            y[outside_obs],
                                            split_time)
  }
  parent <- update_counts(parent)
  return(parent)
}

#' @importFrom purrr map2
update_node_extent_parallel <- function(node, lower, e_lower, upper, e_upper, factor_uniques, e_factor_uniques) {
  node[["B_j"]][["lower"]] <- lower - e_lower
  node[["B_j"]][["upper"]] <- upper + e_upper
  node[["B_j"]][["factor_uniques"]] <- map2(factor_uniques, e_factor_uniques, function(x, y) {
    return(c(y,x))
  })
  return(node)
}

recurse_down_tree <- function(lambda, f_scale, node, X_mat, F_mat, Y, x, f, y) {
  if (node$split$split_type == "categorical") {
    left_obs <- f[, node$split$split_dimension] != node$split$split_point
  }
  else {
    left_obs <- x[, node$split$split_dimension] < node$split$split_point
  }
  node$L <- if (sum(left_obs)==0) node$L else extend_mondrian_block(lambda, f_scale, node$L, X_mat, F_mat, Y,
                                                                     x[left_obs, , drop = FALSE],
                                                                     f[left_obs, , drop = FALSE],
                                                                     y[left_obs],
                                                                     node$split_time)
  node$R <- if (sum(!left_obs)==0) node$R else extend_mondrian_block(lambda, f_scale, node$R, X_mat, F_mat, Y,
                                                                      x[!left_obs, , drop = FALSE],
                                                                      f[!left_obs, , drop = FALSE],
                                                                      y[!left_obs],
                                                                      node$split_time)
  node <- update_counts(node)
  return(node)
}
