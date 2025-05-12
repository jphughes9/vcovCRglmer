library("sandwich")
library("clubSandwich")
library("lme4")
library(MASS)
vcovCR.glmerMod = function(obj, cluster, type="classic"){
  # Check if obj is a fitted model from lmer or glmer
  if ("merMod" %in% class(obj)) {
    stop("The 'obj' should be an object fitted using lmer or glmer.")
  }
  if (!is.null(obj@call$weights)) 
    stop("Models with prior weights are not currently supported.")
  # Check if cluster is manually input and of class factor
  if (!missing(cluster) && !is.factor(cluster)) {
    stop("If 'cluster' is manually input, it must be of class 'factor'.")
  }
  if (missing(cluster)) 
    cluster <- clubSandwich:::get_outer_group(obj)
  if (!clubSandwich:::is_nested_lmerMod(obj, cluster)) 
    stop("Non-nested random effects detected. Method is not available for such models.")
  if (substr(type,1,2)=="FG") {
    if (nchar(type)==2) {
      r = 0.75
    } else {
      r = as.numeric(substr(type,regexpr("\\(",type)[[1]]+1,regexpr("\\)",type)[[1]]-1))
    }
    type1="FG"
  } else {
    type1=type
  }
  # Check if type is one of the specific allowed values
  allowed_types <- c("classic", "DF", "KC", "MD", "FG")
  if (!(type1 %in% allowed_types)) {
    stop("The 'type' must be one of the following: 'classic', 'DF', 'KC', 'MD', 'FG'.")
  }
#################
# extract information from obj
#################
  n = nobs(obj)
  G = length(unique(cluster))
#
  X = model.matrix(obj,type="fixed")
  beta=matrix(fixef(obj),ncol=1)
  np=dim(beta)[1]
#
  Z = model.matrix(obj,type="random")
  nq=dim(Z)[2]
#
  Y = obj@resp$y
# The following allows processing of binomial data
  if (isLMM(obj)) nden=rep(1,length(Y)) else nden = obj@resp$n
#
  eta = predict(obj,type="link")
  ginv_eta = predict(obj,type="response")
#
  link = family(obj)$link
#
  sigma2 = sigma(obj)^2
  lambda = getME(obj,"Lambda")
  R = as.matrix(lambda%*%t(lambda)*sigma2)
  WB_B <- R
##################
# Robust variance calculation
##################
  XtVX = vcov(obj)
  WB_C1 = solve(XtVX)
  sum=matrix(0,np,np)
  for (g in 1:G){
    grp = (cluster == g & nden>0)
    ng = sum(grp)
    if (link == "identity") {
      delta = diag(ng)  
      deltainv = delta 
    } else if (link == "logit") {
      term = ginv_eta[grp]*(1-ginv_eta[grp])
      delta = diag(term)
      deltainv = diag(1/term)
    } else if (link == "log") {
      term = ginv_eta[grp]
      delta = diag(term)
      deltainv = diag(1/term)
    } else {
      stop("Link ",link," not supported")
    }
    #
    P = deltainv%*%(Y[grp]-ginv_eta[grp]) + eta[grp]
    e = matrix(P - X[grp,]%*%beta,ncol=1)
    #
    Sigma = diag(sigma2*family(obj)$variance(ginv_eta[grp])/nden[grp])
    
    # this is diagonal, which is the first term of WB
    mtx_DA <- function(D,A) {
      matrix(rep(diag(D),ncol(A))*as.numeric(A), ncol=ncol(A))
    }
    mtx_AD <- function(A,D) {
      matrix(rep(diag(D), each=nrow(A))*as.numeric(A), ncol=ncol(A))
    }
    if (link=="identity") {
      WB_A <- diag(1/diag(Sigma))
    } else {
      WB_A <- diag(1/diag(mtx_AD(mtx_DA(deltainv,Sigma),deltainv)))
    }
    
    WB_U <- Z[grp,] 
    WB_Ut <- t(Z[grp,])
    # Compute the inverse
#    WB_AUB <- mtx_AD(mtx_DA(WB_A,WB_U),WB_B)
    WB_AUB <- mtx_DA(WB_A,WB_U)%*%WB_B
    WB_UtA <- mtx_AD(WB_Ut,WB_A)
    Vinv <- WB_A - WB_AUB%*%solve(diag(nq) + WB_Ut%*%WB_AUB)%*%WB_UtA
#    
    if (type1=="MD") {
      WB_A = diag(ng)
      WB_U = -t(Vinv)%*%X[grp,]
      WB_V = t(X[grp,])
# Since WB_A is identity, the following expressions are simplified from general Woodbury
      O = WB_C1 + WB_V%*%WB_U
      FF = WB_A - WB_U%*%ginv(matrix(as.numeric(O),dim(O)))%*%WB_V
      sum = sum + t(X[grp,])%*%Vinv%*%t(FF)%*%e[,,drop=FALSE]%*%t(e[,,drop=FALSE])%*%FF%*%Vinv%*%X[grp,]
    } else {
      if (type1=="KC") {
      WB_A = diag(ng)
      WB_U = -t(Vinv)%*%X[grp,]
      WB_V = t(X[grp,])
# Since WB_A is identity, the following expressions are simplified from general Woodbury
      O = WB_C1 + WB_V%*%WB_U
      FF = WB_A - WB_U%*%ginv(matrix(as.numeric(O),dim(O)))%*%WB_V
      sum = sum + (t(X[grp,])%*%Vinv%*%t(FF)%*%e[,,drop=FALSE]%*%t(e[,,drop=FALSE])%*%Vinv%*%X[grp,] + t(X[grp,])%*%Vinv%*%e[,,drop=FALSE]%*%t(e[,,drop=FALSE])%*%FF%*%Vinv%*%X[grp,])/2
    } else {
      if (type1=="FG") {
      Q = t(X[grp,])%*%Vinv%*%X[grp,]%*%XtVX
      AA = diag(1/sqrt(1-pmin(r,diag(Q))))
      sum = sum + AA%*%t(X[grp,])%*%Vinv%*%e[,,drop=FALSE]%*%t(e[,,drop=FALSE])%*%Vinv%*%X[grp,]%*%AA
    } else {
      sum = sum + t(X[grp,])%*%Vinv%*%e[,,drop=FALSE]%*%t(e[,,drop=FALSE])%*%Vinv%*%X[grp,] 
    }}}
#
  }
  c = 1
  if (type1=="DF") {
    if (n-np>0) c= n/(n-np) else cat("DF not valid because n-np <= 0; defaulting to classic")
  }
  robustVar = c*XtVX%*%sum%*%XtVX
  robustVar
}


