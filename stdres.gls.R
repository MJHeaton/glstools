stdres.gls <- function(glsobj){
  
  ## If original model has a variance structure then construct the diagonal
  ## matrices accordingly
  if("varStruct"%in%names(glsobj$modelStruct)){
    
    Dinv <- varWeights(glsobj$modelStruct$varStruct)*(1/sigma(glsobj))
    
  } else {
    
    ## No variance model structure then weights are all 1
    norig <- nrow(eval(glsobj$call$data))
    Dinv <- rep(1,norig)*(1/sigma(glsobj))

  } ## End if("varStruct"%in%names(glsobj$modelStruct))
  
  if("corStruct"%in%names(glsobj$modelStruct)){
    
    # cor.pars <- coef(glsobj$modelStruct$corStruct,unconstrained=FALSE)
    # cor.call <- deparse(glsobj$call$correlation)
    # cor.call <- paste(substr(cor.call,1,nchar(cor.call)-1),", value = c(",
    #                   paste(as.character(cor.pars),collapse=","),"),fixed=TRUE)")
    # Linv <- corMatrix(Initialize(eval(parse(text=cor.call)),
    #                              data=eval(glsobj$call$data)), 
    #                   corr=FALSE)
    Linv <- corMatrix(glsobj$modelStruct$corStruct, corr=FALSE)
    if(is.null(glsobj$groups) | length(unique(glsobj$group))==1){
      #decorr.resid <- as.numeric(solve(t(chol(R)))%*%glsobj$residuals)*Dinv
      decorr.resid <- as.numeric(Linv%*%glsobj$residuals)*Dinv
    } else {
      
      resid.list <- base::split(residuals(glsobj), getGroups(glsobj))
      decorr.resid <- sapply(1:length(resid.list), function(gind){
        #as.numeric(solve(t(chol(R[[gind]])))%*%resid.list[[gind]])
        as.numeric(Linv[[gind]]%*%resid.list[[gind]])
      })
      decorr.resid <- c(decorr.resid)*Dinv
      
    } ## End groups if statement
    
  } else {
    
    decorr.resid <- Dinv*glsobj$residuals
    
  } ## End if("corStruct"%in%names(glsobj$modelStruct))
  
  return(decorr.resid)
  
} #End stdres.gls function