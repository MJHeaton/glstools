library(RSpectra)

## Get Moran Basis Functions w.r.t all Xs
moranBasis <- function(X, A, tol=0.95){
  po <- (diag(nrow(X)) - X%*%chol2inv(chol(t(X)%*%X))%*%t(X)) 
  e <-  RSpectra::eigs_sym(po%*%A%*%po, which="LA", k=nrow(po)-1)
  pos <- which(e$values>0)
  sel <- which(cumsum(e$values[pos])/sum(e$values[pos])<tol) #Use Tol instead
  M <- e$vectors[,pos[sel]]
  M <- as.data.frame(M)
  names(M) <- paste0("B", 1:ncol(M))
  return(M)
}
