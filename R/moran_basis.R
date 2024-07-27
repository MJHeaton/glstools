library(RSpectra)

#' Calculate Moran Basis Functions with respect to all X's
#'
#' @param X A model matrix representing the variables in the spatial model.
#' @param A An adjacency matrix representing the spatial relationships between observations.
#' @param tol Tolerance level
#'
#' @return Dataframe of Moran basis functions
#' @export
#'
#' @examples
#' X <- model.matrix(~ x + y + z, data = data)
#' A <- nb2mat(poly2nb(shape_file), style="B")
#' M <- moranBasis(X, A, tol=0.95)
#' moran_data <- bind_cols(data, M)

moranBasis <- function(X, A, tol = 0.95) {
  po <- (diag(nrow(X)) - X %*% chol2inv(chol(t(X) %*% X)) %*% t(X))
  e <- RSpectra::eigs_sym(po %*% A %*% po, which = "LA", k = nrow(po) - 1)
  pos <- which(e$values > 0)
  sel <- which(cumsum(e$values[pos]) / sum(e$values[pos]) < tol)
  M <- e$vectors[, pos[sel]]
  M <- as.data.frame(M)
  names(M) <- paste0("B", 1:ncol(M))
  return(M)
}
