#' Calculate Standardized Residuals for GLS Model
#'
#' This function calculates standardized residuals for a given fitted GLS, accounting for
#' variance and correlation structures that may be present in the model.
#'
#' @param glsobj A fitted GLS model object.
#'
#' @return Vector of standardized residuals.
#' @export

stdres.gls <- function(glsobj) {
  ## If original model has a variance structure then construct the diagonal
  ## matrices accordingly
  if ("varStruct" %in% names(glsobj$modelStruct)) {
    Dinv <- varWeights(glsobj$modelStruct$varStruct) * (1 / sigma(glsobj))
  } else {
    ## No variance model structure then weights are all 1
    norig <- nrow(eval(glsobj$call$data))
    Dinv <- rep(1, norig) * (1 / sigma(glsobj))
  }
  
  if ("corStruct" %in% names(glsobj$modelStruct)) {
    Linv <- corMatrix(glsobj$modelStruct$corStruct, corr = FALSE)
    if (is.null(glsobj$groups) | length(unique(glsobj$group)) == 1) {
      decorr.resid <- as.numeric(Linv %*% glsobj$residuals) * Dinv
    } else {
      resid.list <- base::split(residuals(glsobj), getGroups(glsobj))
      decorr.resid <- sapply(1:length(resid.list), function(gind) {
        as.numeric(Linv[[gind]] %*% resid.list[[gind]])
      })
      decorr.resid <- c(decorr.resid) * Dinv
    }
  } else {
    decorr.resid <- Dinv * glsobj$residuals
  }
  return(decorr.resid)
}
