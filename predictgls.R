library(nlme)

############################################
## Function to predict with gls correctly ##
## by exploiting correlation              ##
############################################
predictgls <- function(glsobj, newdframe=NULL, level=0.95){
  
  ## If no new dataframe provided, used the dataframe from glsobj
  data.char <- as.character(glsobj$call$data)
  if(is.null(newdframe)){
    newdframe <- get(data.char)
  }
  
  ## Get point predictions of new data frame
  ## Need to break apart formula, remove response then rebuild formula
  the.form <- as.formula(glsobj$call$model)
  the.terms <- terms(the.form,data=newdframe)
  the.terms <- delete.response(the.terms)
  the.form <- as.formula(paste0(the.terms))
  Xpred <- model.matrix(the.form,model.frame(the.terms,data=newdframe,na.action="na.pass"))
  predVals <- c(Xpred%*%glsobj$coefficients)
  
  ## Create a joint dataframe of observed and predicted data
  jdframe <- rbind(get(data.char)[names(newdframe)],newdframe)
  
  ## If original model has a variance structure then construct the diagonal
  ## matrices accordingly
  if("varStruct"%in%names(glsobj$modelStruct)){
    
    ## Get Parameters of Variance structure
    var.pars <- coef(glsobj$modelStruct$varStruct,unconstrained=FALSE)
    
    ## If fixed variance structure then there will be no parameters
    if(length(var.pars)==0){
      ## Initialize var matrix using joint data frame
      var.call <- deparse(glsobj$call$weights)
      
    } else {
      ## Initialize var matrix using joint data frame
      var.call <- deparse(glsobj$call$weights)
      var.call <- paste(substr(var.call,1,nchar(var.call)-1),", value = c(",
                      paste(paste(names(var.pars),as.character(var.pars),sep="="),
                            collapse=","),"))")
    }
      
    #Initialize weights
    varMat.i <- Initialize(eval(parse(text=var.call)),data=jdframe)
    
    ## Create the SD weights
    norig <- nrow(get(data.char))
    W <- varWeights(varMat.i)
    W1 <- W[1:norig]
    W2 <- W[-(1:norig)]
    
  } else {
    
    ## No variance model structure then weights are all 1
    norig <- nrow(get(data.char))
    W1 <- rep(1,norig)
    W2 <- rep(1,nrow(jdframe)-norig)
    
  } ## End if("varStruct"%in%names(glsobj$modelStruct))
  
  ## If original model has correlation structure then construct joint
  ## correlation matrix and calculate conditional correlation matrix
  if("corStruct"%in%names(glsobj$modelStruct)){
    ## Get Parameters of correlation structure
    cor.pars <- coef(glsobj$modelStruct$corStruct,unconstrained=FALSE)
    if(as.character(glsobj$call$correlation)[1]=='corSymm'){
      warning(paste("The general correlation structure corSymm() cannot be used",
              "for prediction.  Returning point prediction without correlation."))
      return(cbind(newdframe,data.frame(Prediction=predVals)))
    } else {
      ## Initialize joint correlation matrix using joint data frame
      cor.call <- deparse(glsobj$call$correlation)
      cor.call <- paste(substr(cor.call,1,nchar(cor.call)-1),", value = c(",
                      paste(as.character(cor.pars),collapse=","),"),fixed=TRUE)")
      corMat.i <- Initialize(eval(parse(text=cor.call)),data=jdframe)
      
      ## If there are no groups (all data belong to same group) then create a joint
      ## correlation matrix.  If there are groups (each group independent) then loop
      ## through the groups calculating correlation each time.
      if(is.null(glsobj$groups) | length(unique(glsobj$group))==1){
        norig <- nrow(get(data.char))
        corMats <- corMatrix(corMat.i)
        nr <- nrow(corMats)
        predMat <- corMats[(norig+1):nr,1:norig]%*%chol2inv(chol(corMats[(1:norig),(1:norig)]))
        allpreds <- predVals+(1/W2)*(predMat%*%(glsobj$residuals*W1))
        allpreds.se <- (1/W2)*sqrt((1-rowSums(predMat*corMats[(norig+1):nr,1:norig])))
      } else {
        pred.grps <- sort(getGroups(newdframe,form=as.formula(glsobj$call$correlation$form)))
        ugrps <- as.character(unique(pred.grps))
        norig <- table(glsobj$groups)
        corMats <- corMatrix(corMat.i)[ugrps]
        getpred <- function(x){
          nr <- nrow(corMats[[x]])
          predMat <- corMats[[x]][(norig[x]+1):nr,1:norig[x]]%*%
            chol2inv(chol(corMats[[x]][(1:norig[x]),(1:norig[x])]))
          thepred <- predVals[pred.grps==ugrps[x]]+(1/W2[pred.grps==ugrps[x]])*(predMat%*%(glsobj$residuals[glsobj$groups==ugrps[x]]*W1[glsobj$groups==ugrps[x]]))
          thepred.se <- sqrt((1-rowSums(predMat*corMats[[x]][(norig[x]+1):nr,1:norig[x]])))*(1/W2[pred.grps==ugrps[x]])
          return(list(thepred=thepred,thepred.se=thepred.se))
        }
        tmp <- lapply(1:length(ugrps),getpred)
        allpreds <- rep(0,nrow(newdframe))
        allpreds.se <- rep(0,nrow(newdframe))
        for(gr in 1:length(ugrps)){
          allpreds[pred.grps==ugrps[gr]] <- tmp[[gr]]$thepred
          allpreds.se[pred.grps==ugrps[gr]] <- tmp[[gr]]$thepred.se
        }
      }
    } # End else
  } else {
    allpreds <- predVals
    allpreds.se <- (1/W2)*rep(1,nrow(newdframe))
  }#End if "corStruct"%in%names(glsobj$modelStruct) statement
  
  ## Calculate upper and lower interval limits
  if(level>1){
    level <- level/100
  }
  alpha <- 1-level
  n <- nrow(get(data.char))
  P <- length(coef(glsobj))
  low <- allpreds - qt(1-alpha/2, df=n-P)*glsobj$sigma*allpreds.se
  up <- allpreds + qt(1-alpha/2, df=n-P)*glsobj$sigma*allpreds.se
  
  ## Return the predictions and predictive SE
  return(cbind(newdframe,data.frame(Prediction=allpreds,
                                    SE.pred=glsobj$sigma*allpreds.se,
                                    lwr=low,
                                    upr=up)))
}
