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
  if (link == "identity") {
    delta = diag(n)  
    deltainv = delta 
  } else if (link == "logit") {
    term = ginv_eta*(1-ginv_eta)
    delta = diag(term)
    deltainv = diag(1/term)
  } else if (link == "log") {
    term = ginv_eta
    delta = diag(term)
    deltainv = diag(1/term)
  } else {
    stop("Link ",link," not supported")
  }
#
  P = deltainv%*%(Y-ginv_eta) + eta
  e = matrix(P - X%*%beta,ncol=1)
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
    Sigma = diag(sigma2*family(obj)$variance(ginv_eta[grp])/nden[grp])
    WB_A <- diag(1/diag(deltainv[grp,grp]%*%Sigma%*%deltainv[grp,grp])) # this is diagonal, which is the first term of WB
    WB_U <- Z[grp,] 
    WB_Ut <- t(Z[grp,])
    # Compute the inverse
    Vinv <- WB_A - WB_A%*%WB_U%*%WB_B%*%solve(diag(nq) + WB_Ut%*%WB_A%*%WB_U%*%WB_B)%*%WB_Ut%*%WB_A
#    
    if (type1=="MD") {
      WB_A = diag(ng)
      WB_U = -t(Vinv)%*%X[grp,]
      WB_V = t(X[grp,])
# Since WB_A is identity, the following expressions are simplified from general Woodbury
      O = WB_C1 + WB_V%*%WB_U
      FF = WB_A - WB_U%*%ginv(matrix(as.numeric(O),dim(O)))%*%WB_V
      sum = sum + t(X[grp,])%*%Vinv%*%t(FF)%*%e[grp,,drop=FALSE]%*%t(e[grp,,drop=FALSE])%*%FF%*%Vinv%*%X[grp,]
    } else {
      if (type1=="KC") {
      WB_A = diag(ng)
      WB_U = -t(Vinv)%*%X[grp,]
      WB_V = t(X[grp,])
# Since WB_A is identity, the following expressions are simplified from general Woodbury
      O = WB_C1 + WB_V%*%WB_U
      FF = WB_A - WB_U%*%ginv(matrix(as.numeric(O),dim(O)))%*%WB_V
      sum = sum + (t(X[grp,])%*%Vinv%*%t(FF)%*%e[grp,,drop=FALSE]%*%t(e[grp,,drop=FALSE])%*%Vinv%*%X[grp,] + t(X[grp,])%*%Vinv%*%e[grp,,drop=FALSE]%*%t(e[grp,,drop=FALSE])%*%FF%*%Vinv%*%X[grp,])/2
    } else {
      if (type1=="FG") {
      Q = t(X[grp,])%*%Vinv%*%X[grp,]%*%XtVX
      AA = diag(1/sqrt(1-pmin(r,diag(Q))))
      sum = sum + AA%*%t(X[grp,])%*%Vinv%*%e[grp,,drop=FALSE]%*%t(e[grp,,drop=FALSE])%*%Vinv%*%X[grp,]%*%AA
    } else {
      sum = sum + t(X[grp,])%*%Vinv%*%e[grp,,drop=FALSE]%*%t(e[grp,,drop=FALSE])%*%Vinv%*%X[grp,] 
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


