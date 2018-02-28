#' Creates matrices of covariates.
#'
#' Internal function to separate a dataset into matrices of numeric and factor variables.
#'
#' @param X Data (matrix or data frame).
#' @return A list of matrices for numeric and factor variables (if any).
#' @examples
#' set.seed(1)
#' convert_to_matrices(data.frame(x1=rnorm(100),
#'                                x2=rbeta(100, 1, 2),
#'                                f1=sample(letters, 100, replace = TRUE)))

convert_to_matrices <- function(X) {
  if (is.data.frame(X)) {
    numeric_vars <- which(purrr::map_lgl(X, function(x) is.numeric(x)))
    factor_vars <- setdiff(1:ncol(X), numeric_vars)

    return(list("X_mat"=as.matrix(X[, numeric_vars, drop = FALSE]),
                "F_mat"=as.matrix(X[, factor_vars, drop = FALSE])))
  }
  else {
    if (is.numeric(X)) {return(list("X_mat"=as.matrix(X), "F_mat"=NULL))}
    else {return(list("X_mat"=NULL, "F_mat"=as.matrix(X)))}
  }
}

get_vector_table <- function(factor) {
  mytable <- table(factor)
  counts <- as.integer(mytable)
  names(counts) <- names(mytable)
  return(counts)
}

named_pmin <- function(x_named, y) {
  min <- pmin.int(x_named, y)
  names(min) <- names(x_named)
  return(min)
}

rexp_or_lambda <- function(rate, lambda) {
  if (rate==0) {
    E <- lambda # Leaf
  }
  else {
    E <- rexp(1, rate) # sample from exponential(lambda) distribution
  }
  return(E)
}

#' Builds a MondrianForest.
#'
#' \code{mondrian_forest} implements Lakshminarayanan et al's Mondrian Tree
#' algorithms described in [Lakshminarayanan et al. 2014]
#' (\url{https://arxiv.org/abs/1406.2673}) with a modification allowing for more
#' space-efficient dummy variable treatment of categorical variables.
#'
#' @param X Data (matrix or data frame) containing features and column of labels.
#' @param y_col_num Numeric length 1 vector of column number of label (defaults to last
#'   ncol(X)).
#' @param lambda Budget parameter, see [Lakshminarayanan et al. 2014].
#' @param f_scale Numeric length 1 vector if constant or length equal to the number
#'   of categorical variables. \code{f_scale} represents the implicit range of the
#'   dummy-encoded categorical variables.
#' @param ntree Numeric length 1 vector of number of trees to build in forest.
#' @param verbose Boolean length 1 vector, prints additional information while algorithm is running
#'   (currently just time to build trees).
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
#' mf <- mondrian_forest(test[1:750, ], y_col_num = 5, lambda = 3)
#' table(test$label[751:1000], predict(mf, test[751:1000, ], type = "class"))
#' # Compare to Random Forest in installed:
#' # rf <- randomForest::randomForest(test[1:750, -5],
#'                                    y = test$label[1:750],
#'                                    ntree = 1000)
#' # table(test$label[751:1000], predict(rf, test[751:1000, ]))
#' @export

mondrian_forest <- function(X,
                            y_col_num,
                            lambda,
                            f_scale = 1,
                            ntree = 25,
                            verbose=FALSE) {
  if (missing(y_col_num)) { # Check for given label
    y_col_num <- ncol(X)
    names(y_col_num) <- colnames(X)[y_col_num]
    warning("y_col_num not given, assuming label is the last column in X.\n")
  }
  else {names(y_col_num) <- colnames(X)[y_col_num]}
  if (is.data.frame(X)) {Y <- as.factor(X[[y_col_num]])} # Get vector, turn to factor
  else {Y <- as.factor(X[, y_col_num, drop = TRUE])}
  X <- X[, -y_col_num, drop = FALSE] # Keep as D.F/Matrix

  X_F <- convert_to_matrices(X) # Turn to matrix

  X_F <- purrr::map(X_F, function(mat) {
    rownames(mat) <- 1:nrow(mat)
    return(mat)
  })

  if (!(length(f_scale) %in% c(1, ncol(X_F$F_mat)))) {
    stop("f_scale must have length 1 or length equal to the number of categorical variables.")
  }

  mforest <- list(forest = list(),
                  "data" = list(),
                  parameters = list("lambda"=lambda, "f_scale"=f_scale, "y_col_num"=y_col_num))
  for (i in 1:ntree) {
    t1 <- proc.time()
    mforest$forest[[paste("Tree", i, sep="_")]] <- sample_mondrian_tree(lambda,
                                                                         f_scale,
                                                                         X_F$X_mat,
                                                                         X_F$F_mat,
                                                                         Y)
    if (verbose==TRUE) {
      t_diff <- proc.time() - t1
      cat("Tree ", i, ": ", "built in ", round(unname(t_diff[3]), 1), " seconds", "\n", sep = "")
    }
  }
  mforest[["data"]][["X_F"]] <- X_F # Store the dataset so it can be referenced by the leaves
  mforest[["data"]][["Y"]] <- Y
  class(mforest) <- "mondrianforest"
  return(mforest)
}

sample_dimension_split_pt <- function(upper, lower, factor_uniques, factor_lengths, f_scale) { # Make sample dimension and split point
  variable_type <- sample(x = c("numeric", "categorical"),
                          size = 1,
                          prob = c(sum(upper-lower), sum((factor_lengths-1)*f_scale)))
  if (variable_type=="numeric") {
    split_dimension <- sample(seq_along(upper), size = 1, prob = upper - lower)
    split_point <- runif(n = 1, min = lower[split_dimension], max = upper[split_dimension])
    return(list("split_dimension"=split_dimension,
                "split_point"=split_point,
                "split_type"=variable_type))
  }
  else {
    split_dimension <- sample(seq_along(factor_lengths), size = 1, prob = factor_lengths - 1)
    split_category <- sample(factor_uniques[[split_dimension]], size = 1)
    return(list("split_dimension"=split_dimension,
                "split_point"=split_category,
                "split_type"=variable_type))
  }
}

sample_mondrian_tree <- function(lambda, f_scale, X_mat, F_mat, Y) {
  tree <- sample_mondrian_block(lambda, f_scale, X_mat, F_mat, Y, split_time_parent = 0)
  return(tree)
}

initialize_node_and_extent <- function(X_j, F_j) {
  node <- list()
  lower <- upper <- vector(length = ncol(X_j))
  factor_uniques <- list()
  for (i in seq_len(ncol(X_j))) {
    lower[i] <- min(X_j[, i])            # Dimension-wise min and max
    upper[i] <- max(X_j[, i])
  }
  for (i in seq_len(ncol(F_j))) {
    factor_uniques[[i]] <- as.character(unique(F_j[, i]))   # Dimension-wise number of uniques
  }
  node[["B_j"]] <- list("lower"=lower, "upper"=upper, "factor_uniques"=factor_uniques)
  return(node)
}

pause_leaf <- function(node, X_j) {
  node[["pause"]] <- TRUE
  node[["obs"]] <- as.integer(rownames(X_j))
  return(node)
}

create_split <- function(node, split_time, factor_lengths, f_scale) {
  node[["split_time"]] <- split_time
  node[["split"]] <- sample_dimension_split_pt(node$B_j$upper,
                                               node$B_j$lower,
                                               node$B_j$factor_uniques,
                                               factor_lengths,
                                               f_scale)
  return(node)
}

create_leaf <- function(node, lambda, Y_j) {
  node[["leaf"]] <- TRUE
  node[["split_time"]] <- lambda
  return(node)
}

compute_counts <- function(node, Y_j) {
  if (!missing(Y_j)) {
    node[["cust"]] <- get_vector_table(Y_j)              # Cust and Tab are needed for predictions.
    node[["tab"]] <- named_pmin(node$cust, 1L)
    return(node)
  }
  else {
    node[["cust"]] <- node$L$tab + node$R$tab
    node[["tab"]] <- named_pmin(node$cust, 1L)
    return(node)
  }
}
#' @importFrom purrr map_int
sample_mondrian_block <- function(lambda, f_scale, X_j, F_j, Y_j, split_time_parent) {
  # X_j and Y_j should be interpreted as the data associated with the CURRENT node.
  node <- initialize_node_and_extent(X_j, F_j) # Returns node and B_j
  if (length(unique(Y_j)) == 1){ # Pause conditions
    node <- pause_leaf(node, X_j)
    node <- create_leaf(node, lambda, Y_j)
    node <- compute_counts(node, Y_j)
    return(node)
  }
  else { # Not paused
    factor_lengths <- map_int(node$B_j$factor_uniques, length)
    E <- rexp_or_lambda(sum(node$B_j$upper - node$B_j$lower) + sum((factor_lengths - 1)*f_scale), lambda) # sample from exponential(rate) distribution
    if (split_time_parent + E < lambda) { # If current "time" is less than the budget, recurse and make deeper nodes
      node <- create_split(node, split_time_parent + E, factor_lengths, f_scale)
      if (node$split$split_type=="numeric") {
        left_obs <- X_j[, node$split$split_dimension] < node$split$split_point
      }
      else {
        left_obs <- F_j[, node$split$split_dimension] != node$split$split_point
      }
      node$L <- sample_mondrian_block(lambda,
                                       f_scale,
                                       X_j[left_obs, , drop = FALSE],
                                       F_j[left_obs, , drop = FALSE],
                                       Y_j[left_obs],
                                       split_time_parent = node$split_time)
      node$R <- sample_mondrian_block(lambda,
                                       f_scale,
                                       X_j[!left_obs, , drop = FALSE],
                                       F_j[!left_obs, , drop = FALSE],
                                       Y_j[!left_obs],
                                       split_time_parent = node$split_time)
      node <- compute_counts(node)
    }
    else { # Else, create a leaf
      node <- create_leaf(node, lambda, Y_j)
      node <- compute_counts(node, Y_j)
    }
    return(node)
  }
}
