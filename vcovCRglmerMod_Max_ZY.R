vcovCR.glmerMod_Max_ZY = function(obj, cluster = NULL, type="classic") {
  # Check if obj is a fitted model from lmer or glmer
  if ("merMod" %in% class(obj)) {
    stop("The 'obj' should be an object fitted using lmer or glmer.")
  }
  
  if (is.null(cluster)) {
    cluster <- clubSandwich:::get_outer_group(obj)
  }
  else if (!is.factor(cluster)) {
    stop("If 'cluster' is manually input, it must be of class 'factor'.")
  }
  
  # Initialize default r for FG method
  r <- 3/4
  method <- type
  
  # Check if the type contains parameters (i.e., parentheses for FG)
  if (grepl("\\(", type)) {
    if (!grepl("^FG\\(", type)) {
      stop("Only the 'FG' method need an r value input.")
    }
    
    method <- gsub("\\(.*$", "", type)  # Extract method name before the parenthesis
    r_str <- gsub(".*\\(|\\)", "", type)  # Extract between parentheses
    r <- as.numeric(r_str)  # Convert to numeric
    
    # Validate r
    if (is.na(r) || r < 0 || r >= 1) {
      stop("The 'r' value must be a numeric between 0 (inclusive) and 1 (exclusive).")
    }
  }
  
  # Validate the method name
  allowed_types <- c("classic", "DF", "KC", "MD", "FG")
  if (!(method %in% allowed_types)) {
    stop("The 'type' must be one of the following: 'classic', 'DF', 'KC', 'MD', 'FG'.")
  }
  
  # essential matrices
  beta=matrix(fixef(obj),ncol=1)
  np=dim(beta)[1]

  X = model.matrix(obj,type="fixed")
  Z = model.matrix(obj,type="random")
  Y = obj@resp$y # example of "slots" in R
  # The following allows processing of binomial data
  if (isLMM(obj)) nden=rep(1,length(Y)) else nden = obj@resp$n #BUG fix

  nq = dim(Z)[2]
  
  eta = predict(obj,type="link")
  ginv_eta = predict(obj,type="response")
  
  # checking link function
  link = family(obj)$link
  switch(link,
         "identity" = {
           delta = diag(nobs(obj))
           deltainv = delta
         },
         "logit" = {
           delta = diag(exp(eta) / (1 + exp(eta))^2)
           deltainv = diag((1 + exp(eta))^2 / exp(eta))
         },
         "log" = {
           delta = diag(exp(eta))
           deltainv = diag(1 / exp(eta))
         },
         {
           stop("Model object must have an identity, logit, or log link")
         }
  )
  
  # checks for DF
  c = ifelse(type == "DF", dim(X)[1] / (dim(X)[1] - dim(X)[2]), 1)
  
  theta = as.data.frame(VarCorr(obj))
  G = length(levels(cluster))
  sigma2 = sigma(obj)^2
  sigma_mat = family(obj)$variance(ginv_eta) * sigma2 / nden  # BUG fix
  lambda = getME(obj,"Lambda")
  R = lambda%*%t(lambda)*sigma2
  mbv = vcov(obj)
  P = deltainv %*% (Y - ginv_eta) + eta
  e = P - X %*% beta
  
  # sum product of matrices over over clusters
  running_sum = matrix(0, nrow = dim(X)[2], ncol = dim(X)[2])
  WB_B = R
  for (clus in unique(cluster)) {
    clus_idx = which(cluster == clus & nden > 0)    # BUG fix, in case there are 0 observations
    WB_A = diag(1 / diag(deltainv[clus_idx, clus_idx] %*% diag(sigma_mat[clus_idx]) %*% deltainv[clus_idx, clus_idx]))
    WB_U = Z[clus_idx, ]
    WB_Ut <- t(Z[clus_idx,])
    # Compute the inverse
    Vinv_g <- WB_A - WB_A%*%WB_U%*%WB_B%*%solve(diag(nq) + WB_Ut%*%WB_A%*%WB_U%*%WB_B)%*%WB_Ut%*%WB_A

    X_g = X[clus_idx, ]
    e_g = e[clus_idx, ]
    switch(type,
           "KC" = {
             H_g = X_g %*% mbv %*% t(X_g) %*% Vinv_g
             #Using Woodbury
             WB_A2 <- diag(nrow(H_g))
             WB_C2 <- t(solve(mbv))
             WB_U2 <- t(Vinv_g) %*% X[clus_idx, ]
             WB_V2 <- t(X[clus_idx, ])
             W2 <- WoodburyMatrix(A=WB_A2, B=WB_C2, U=-WB_U2, V=WB_V2)
             
             I_minus_Hg_prime_inv <- solve(W2)
             
             F_g <- sqrtm(I_minus_Hg_prime_inv)
             running_sum = running_sum + t(X_g) %*% Vinv_g %*% t(F_g) %*% e_g %*% t(e_g) %*% F_g %*% Vinv_g %*% X_g
           },
           "MD" = {
             H_g = X_g %*% mbv %*% t(X_g) %*% Vinv_g
             #Using Woodbury
             WB_A2 <- diag(nrow(H_g))
             WB_C2 <- t(solve(mbv))
             WB_U2 <- t(Vinv_g) %*% X[clus_idx, ]
             WB_V2 <- t(X[clus_idx, ])
             W2 <- WoodburyMatrix(A=WB_A2, B=WB_C2, U=-WB_U2, V=WB_V2)
             
             F_g <- solve(W2)
             
             running_sum = running_sum + t(X_g) %*% Vinv_g %*% t(F_g) %*% e_g %*% t(e_g) %*% F_g %*% Vinv_g %*% X_g
           },
           "FG" = {
             Q_g = t(X_g) %*% Vinv_g %*% X_g %*% mbv
             A_g = diag(sqrt(1 - pmin(r, diag(Q_g))))
             running_sum = running_sum + A_g %*% t(X_g) %*% Vinv_g %*% e_g %*% t(e_g) %*% Vinv_g %*% X_g %*% A_g
           },
           {
             running_sum = running_sum + t(X_g) %*% Vinv_g %*% e_g %*% t(e_g) %*% Vinv_g %*% X_g
           }
    )
  }
  
  # final robust variance formula
  robustVar = c * mbv %*% running_sum %*% mbv
  return(robustVar)
}