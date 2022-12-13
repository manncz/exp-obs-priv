# script: 00-helper-funs.R
# author: Charlotte Mann
# date: Dec 13, 2022
# purpose: Functions used in simulation studies

# MSE
mse_fun <- function(pred, out){
  
  mse <- mean((pred-out)^2)
  
  return(mse)
  
}

# Expit and logit functions
expit <- function(x){
  1/(1+exp(-x))
}

logit <- function(p){
  log(p/(1-p))
}

# Calculate scale parameter for Laplace noise in DP Laplace Mechanism
calc_dp_err_lap <- function(X, epsilon = log(3),delta =NULL,
                            statistic){
  
  if(statistic == "corr"){
    delta_f = rep(NA, nrow(X))
  }else{
    delta_f = matrix(rep(NA, nrow(X)*ncol(X)), nrow = nrow(X))
  }
  
  # calculate vector of statistics for desired moment
  if(statistic == "corr"){
    Xs <- scale(X)
    G <- t(Xs) %*% Xs
    G.compare <- G*upper.tri(G, diag = F)
  }else if(statistic == "mean"){
    E.compare <- apply(X, 2, mean)
  }else if(statistic == "sd"){
    V.compare <- apply(X, 2, sd)
  }
    
  delta_f <- foreach(i = 1:nrow(X), .combine = rbind) %dopar% {
    
    X_i <- as.matrix(X[-i,])
    
    if(statistic == "corr"){
      X_i <- scale(X_i)
      G_i <- t(X_i)%*%X_i
      delta_f <- sum(abs(G_i*upper.tri(G, diag = T)-G.compare))
    }else if(statistic == "mean"){
      E_i = apply(X_i, 2, mean)
      delta_f <- abs(E_i - E.compare)
    }else if(statistic == "sd"){
      V_i = apply(X_i, 2, sd)
      delta_f <- abs(V_i - V.compare)
    }
    
    delta_f
  }
  
  if(statistic == "corr"){
    max_delta_f = max(delta_f)
  }else if(statistic %in% c("mean","sd")){
    max_delta_f = apply(delta_f, 2, max)
    epsilon <- epsilon / ncol(X)
  }
  
  lambda <- max_delta_f/epsilon
  
  return(lambda)
}

# Calculate scale parameter for Gaussian noise in DP Gaussian Mechanism
calc_dp_err_gaus <- function(X, epsilon = log(3), delta = 10^-5,
                             statistic = "corr"){
  
  if(statistic == "corr"){
    delta_f = rep(NA, nrow(X))
  }else{
    delta_f = matrix(rep(NA, nrow(X)*ncol(X)), nrow = nrow(X))
  }
  
  
  # calculate vector of statistics for desired moment
  if(statistic == "corr"){
    Xs <- scale(X)
    G <- t(Xs) %*% Xs
    G.compare <- G*upper.tri(G, diag = F)
  }else if(statistic == "mean"){
    E.compare <- apply(X, 2, mean)
  }else if(statistic == "sd"){
    V.compare <- apply(X, 2, sd)
  }
  
  # for each observation in X, drop the row, recalculate the vector of statistics
  # and calculate the sum of squared differences between that vector and the one
  # that uses all of the observations.
  delta_f <- foreach(i = 1:nrow(X), .combine = rbind) %dopar% {
    
    X_i <- as.matrix(X[-i,])
    
    if(statistic == "corr"){
      X_i <- scale(X_i)
      G_i <- t(X_i)%*%X_i
      G_i <- G_i*upper.tri(G_i, diag = F)
      delta_f <- sum((G_i-G.compare)^2)
    }else if(statistic == "mean"){
      E_i = apply(X_i, 2, mean)
      delta_f <- (E_i - E.compare)^2
    }else if(statistic == "sd"){
      V_i = apply(X_i, 2, sd)
      delta_f <- (V_i - V.compare)^2
    }
    
    delta_f
  }
  
  if(statistic == "corr"){
    max_delta_f = max(delta_f)
  }else if(statistic %in% c("mean","sd")){
    max_delta_f = apply(delta_f, 2, max)
    
    #distribute epsilon and delta budget
    epsilon <- epsilon / ncol(X)
    delta <- delta / ncol(X)
  }
    
  lambda <- max_delta_f/epsilon*sqrt(2*log(1.25/delta))
  
  return(lambda)
}

# Differentially private procedure for a gram matrix 
dp_private_gram <- function(X, epsilon = 3, delta = 10^-5,
                            mean.budget=1/3, sd.budget=1/3, corr.budget=1/3,
                            mechanism = "gaussian", seed = 1234,
                            mean.noise = NULL, sd.noise = NULL, corr.noise = NULL){
  
  set.seed(seed)
  
  p <- ncol(X)
  n <- nrow(X)
  if(mechanism == "gaussian"){
    dp_noise_fun = calc_dp_err_gaus
    #distribute delta according to budget delta
    mean.delta <- delta*mean.budget
    sd.delta<- delta*sd.budget
    corr.delta <- delta*corr.budget
    
  }else{
    dp_noise_fun = calc_dp_err_lap
  }
  
  # distribute epsilon according to budget
  mean.epsilon <- epsilon*mean.budget
  sd.epsilon <- epsilon*sd.budget
  corr.epsilon <- epsilon*corr.budget
  
  # calculate gaussian noise for mean and variance vectors and corr matrix
  #noise for mean and sd will be vectors
  if(is.null(mean.noise)){
    mean.noise <- dp_noise_fun(X, epsilon = mean.epsilon, delta = mean.delta,
                               statistic = "mean")
    sd.noise <- dp_noise_fun(X, epsilon = sd.epsilon, delta = sd.delta,
                             statistic = "sd")
    corr.noise <- dp_noise_fun(X, epsilon = corr.epsilon, delta = corr.delta,
                               statistic = "corr")
  }
  
  #calculate true mean, variance and correlation matrices
  true.mean <- apply(X, 2, mean)
  true.sd <- apply(X, 2, sd)
  X_s <- scale(X)
  corr <- t(X_s) %*% X_s
  corr.mat.indx <- upper.tri(corr, diag = F)
  corr <- corr*corr.mat.indx
  
  #create matrix of gaussian noise matching dimenstions of correlation matrix
  corr.noise.mat <- matrix(rep(0,p^2), nrow = p)
  
  if(mechanism == "gaussian"){
    corr.noise.mat[corr.mat.indx] <- rnorm((p-1)*p/2, mean = 0, sd = sqrt(corr.noise))
    mean.noise.vec <- sapply(sqrt(mean.noise), function(x){rnorm(1, mean = 0, sd = x)})
    sd.noise.vec <- sapply(sqrt(sd.noise), function(x){rnorm(1, mean = 0, sd = x)})
  }else{
    corr.noise.mat[corr.mat.indx] <- rlaplace((p-1)*p/2, scale = corr.noise)
    mean.noise.vec <- sapply(mean.noise, function(x){rlaplace(p, scale = x)})
    sd.noise.vec <- sapply(sd.noise, function(x){rlaplace(p, scale = x)})
  }
 
  # gaussian mechanism for differentially private correlation matrix, mean and variance vectors
  dp.corr <- corr + corr.noise.mat
  dp.mean <- true.mean + mean.noise.vec
  dp.sd  <- true.sd + sd.noise.vec
  
  # transform to non-central moments (empirical E[X_1*X_2] matrix, E[X^2] matrix, E[X] matrix)
  x1x2 <- dp.mean %*% t(dp.mean)*corr.mat.indx
  sig1sig2 <- dp.sd %*% t(dp.sd)*corr.mat.indx
  
  dp.EX1X2 <- dp.corr/n*sig1sig2 + x1x2
  dp.EX2 <- dp.sd^2*(n-1)/n + dp.mean^2
  dp.EX <- dp.mean
  
  #check - should give back the original X'X
  # x1x2 <- true.mean %*% t(true.mean)*corr.mat.indx
  # sig1sig2 <- true.sd %*% t(true.sd)*corr.mat.indx
  # dp.EX1X2 <- corr/n*sig1sig2 + x1x2
  # dp.EX2 <- true.sd^2*(n-1)/n + true.mean^2
  # dp.EX <- true.mean
  
  #reconstruct X'X/n, now differentially private with non-central empirical moments
  dp.G <- dp.EX1X2 + t(dp.EX1X2)
  diag(dp.G) <- dp.EX2
  dp.G <- cbind(dp.EX, dp.G)
  dp.G <- rbind(c(1, dp.EX), dp.G)
  colnames(dp.G)[1] <- "intercept"
  rownames(dp.G)[1] <- "intercept"
  
  # save noise and epsilon values for gaussian mechanism
  names(sd.noise) <- names(mean.noise) <- colnames(X)
  dp.noise <- data.frame(mean.noise, sd.noise,
                         corr = corr.noise)
  
  ep.vec <- c(mean.epsilon, sd.epsilon, corr.epsilon)
  names(ep.vec) <- c("mean", "sd", "corr")
  
  # return differentially private X'X/n, and noise and epsilon parameters
  return(list(G = dp.G, dp.noise = dp.noise, epsilon.budget = ep.vec))
}