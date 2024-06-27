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

# Generate differentially private gram matrix
dp_private_gram <- function(X, epsilon = 3, delta = 10^-5,
                            mean.budget=1/3, e2.budget=1/3, e12.budget=1/3,
                            lbound, ubound,
                            mechanism = "gaussian", seed = 1234){
  
  p <- ncol(X)
  n <- nrow(X)
  
  # divide privacy budget
  mean.ep = epsilon*mean.budget
  var.ep = epsilon*e2.budget
  cov.ep = epsilon*e12.budget/((p+1)*p/2)
  
  mean.delta = delta*mean.budget
  var.delta = delta*e2.budget
  cov.delta = delta*e12.budget/((p+1)*p/2)
  
  covepMat <-  rep(1,p)%*%t(rep(1/cov.ep,p))
  diag(covepMat) <- 1/var.ep
  
  covdeltaMat <-  rep(1,p)%*%t(rep(1/cov.delta,p))
  diag(covdeltaMat) <- 1/var.delta
  
  ########## calculate scale parameter for Gaussian mechanism ###########
  
  # sensitivities based on given vectors of upper and lower bounds for the columns (assumed to have been inputted in order)
  mean.sens <- pmax(abs(ubound),abs(lbound))/n
  cov.sens <- (pmax(abs(ubound),abs(lbound)))%*%t(pmax(abs(ubound),abs(lbound)))/n
  
  mean.scaleparam <-  mean.sens/mean.ep*sqrt(2*log(1.25/mean.delta))
  cov.scaleparam <-  cov.sens*covepMat*sqrt(2*log(1.25*covdeltaMat))*upper.tri(cov.sens, diag = T)
  
  ##########                    inject noise                   ###########
  
  # calculate original covariance matrix and mean based on X, truncated by upper and lower bounds
  lboundMat <- rep(1,n)%*%t(lbound)
  uboundMat <- rep(1,n)%*%t(ubound)
  
  X_trunc <- X*(X >= lboundMat) + lboundMat*(X < lboundMat)
  X_trunc <- X_trunc*(X_trunc <= uboundMat) + uboundMat*(X_trunc > uboundMat)
  
  base.EX1X2 <- t(X_trunc)%*%X_trunc/n
  base.EX <- rep(1,n)%*%X_trunc/n
  
  # generate random noise
  mean.noisevec <- sapply(mean.scaleparam, function(x){rnorm(1, mean = 0, sd = x)})
  
  cov.noisemat <- apply(cov.scaleparam, 1:2,  function(x){rnorm(1, mean = 0, sd = x)})
  cov.noisemat <- cov.noisemat + t(cov.noisemat*upper.tri(cov.noisemat, diag = F))
  rownames(cov.noisemat) <- colnames(cov.noisemat) <- NULL
  #isSymmetric(cov.noisemat)
  
  # add random noise to mean and covariance
  dp.EX <- as.vector(base.EX + mean.noisevec)
  dp.EX1X2 <- base.EX1X2 + cov.noisemat
  #isSymmetric(dp.EX1X2)
  
  ##########                reconstruct matrix                ###########
  dp.G <- cbind(dp.EX, dp.EX1X2)
  dp.G <- rbind(c(1, dp.EX), dp.G)
  colnames(dp.G)[1] <- "intercept"
  rownames(dp.G)[1] <- "intercept"  
  
  ##########              post-process positive-definite                  ###########
  
  #check if positive-definite
  es <- eigen(dp.G, symmetric = T)
  
  #if any negative eigenvalues, reconstruct a positive definite matrix
  if(sum(es$values < 0) > 0){
    
    print("Original DP Gram matrix not positive-definite. Post processing...")
    
    medval <- median(es$values[es$values > 0])
    
    # reset any negative eigenvalues to 0
    poseigvals <- es$values
    poseigvals[es$values < 0] <- 0
    
    # ridge with median positive eigenvalue
    poseigvals <- poseigvals + medval
    
    # reconstruct matrix
    V <- es$vectors
    lambda <- diag(poseigvals)
    
    dp.G.recon <- V%*%lambda%*%solve(V)
    colnames(dp.G.recon) <- rownames(dp.G.recon) <-  colnames(dp.G)
  }else{
    dp.G.recon <- dp.G
  }
  
  
  #returning a list to match the local sensitivity function, just to make integrating the 
  #new function easier to integrate into the simulation functions
  return(list(G = dp.G.recon, dp.noise = NULL))
  
}



################################################################################
 #                                OLD FUNCTIONS                              #
################################################################################

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
      delta_f <- sqrt(sum((G_i-G.compare)^2))
    }else if(statistic == "mean"){
      E_i = apply(X_i, 2, mean)
      delta_f <- sqrt((E_i - E.compare)^2)
    }else if(statistic == "sd"){
      V_i = apply(X_i, 2, sd)
      delta_f <- sqrt((V_i - V.compare)^2)
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
    
  lambda <- (max_delta_f/epsilon)*sqrt(2*log(1.25/delta))
  
  return(lambda)
}


# Differentially private procedure for a gram matrix 
dp_private_gram_local <- function(X, epsilon = 3, delta = 10^-5,
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
  
  # calculate scale parameter for DP mechanism for mean and variance vectors and corr matrix
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
  
  #create matrix of gaussian noise matching dimensions of correlation matrix
  corr.noise.mat <- matrix(rep(0,p^2), nrow = p)
  
  if(mechanism == "gaussian"){
    corr.noise.mat[corr.mat.indx] <- rnorm((p-1)*p/2, mean = 0, sd = corr.noise)
    mean.noise.vec <- sapply(mean.noise, function(x){rnorm(1, mean = 0, sd = x)})
    sd.noise.vec <- sapply(sd.noise, function(x){rnorm(1, mean = 0, sd = x)})
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
  # mean(dp.G - G)
  
  #reconstruct X'X/n, now differentially private with non-central empirical moments
  dp.G <- dp.EX1X2 + t(dp.EX1X2)
  diag(dp.G) <- dp.EX2
  dp.G <- cbind(dp.EX, dp.G)
  dp.G <- rbind(c(1, dp.EX), dp.G)
  colnames(dp.G)[1] <- "intercept"
  rownames(dp.G)[1] <- "intercept"
  
  #check if positive-definite
  es <- eigen(dp.G, symmetric = T)
  
  #if any negative eigenvalues, reconstruct a positive definite matrix
  if(sum(es$values < 0) > 0){
    
    print("Original DP Gram matrix not positive-definite. Post processing...")
    
    medval <- median(es$values[es$values > 0])
    
    # reset any negative eigenvalues to 0
    poseigvals <- es$values
    poseigvals[es$values < 0] <- 0
    
    # ridge with median positive eigenvalue
    poseigvals <- poseigvals + medval
    
    # reconstruct matrix
    V <- es$vectors
    lambda <- diag(poseigvals)
    
    dp.G.recon <- V%*%lambda%*%solve(V)
    colnames(dp.G.recon) <- rownames(dp.G.recon) <-  colnames(dp.G)
  }else{
    dp.G.recon <- dp.G
  }
  
  # save noise and epsilon values for gaussian mechanism
  names(sd.noise) <- names(mean.noise) <- colnames(X)
  dp.noise <- data.frame(mean.noise, sd.noise,
                         corr = corr.noise)
  
  ep.vec <- c(mean.epsilon, sd.epsilon, corr.epsilon)
  names(ep.vec) <- c("mean", "sd", "corr")
  
  # return differentially private X'X/n, and noise and epsilon parameters
  return(list(G = dp.G.recon, dp.noise = dp.noise, epsilon.budget = ep.vec))
}
