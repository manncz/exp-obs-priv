# script: 00-sim-funs.R
# author: Charlotte Mann
# date: Dec 13, 2022
# purpose: Functions for variance simulations

#================================== GENERATE SIMULATED DATA ====================================#

gen_obs_rct_dat <- function(n.obs = 1000, n.rct = 50,
                            P = 10, prop.sig = .8, prop.norm = .6,
                            alpha = 10, beta.bank = 10, unex.sig=2,
                            seed = 12345){
  set.seed(seed)
  
  P.sig <- floor(P*prop.sig) #number of non-null covariates
  P.norm <- floor(P*prop.norm) #number of continuous covariates
  P.cat <- P-P.norm #number of categorical covariates
  
  #generate beta
  beta <- rep(0, P)
  beta[sample(1:(P), P.sig)] <- sqrt(beta.bank/P.sig)
  beta <- c(alpha, beta)
  
  error <- rnorm(n.obs, 0, unex.sig)
  
  #******************** generate obs covariate data *********************#
  
  # generate continuous covariates
  X <- data.frame(matrix(rnorm(n.obs*P.norm,0,1), ncol=P.norm))
  
  # generate categorical covariates
  
  p <- NULL
  if(P.cat != 0){
    p <- runif(P.cat)
    for(i in 1:P.cat){
      temp <- sample(c(0,1), size = n.obs, replace = T, prob = c(1-p[i], p[i]))
      X <- cbind(X,temp)
    }
  }
  
  colnames(X) <- paste0("X", 1:P)
  
  #************************ generate obs full data ************************#
  
  Xm <- model.matrix(~., data = X)
  y <- Xm %*% beta + error
  dat <- data.frame(y, X)
  

  #******************** generate rct covariate data *********************#
  # generate continuous covariates
  X.rct <- data.frame(matrix(rnorm(n.rct*P.norm,0,1), ncol=P.norm))
    
  # generate categorical covariates
  if(P.cat != 0){
    for(i in 1:P.cat){
      temp <- sample(c(0,1), size = n.rct, replace = T, prob = c(1-p[i], p[i]))
      X.rct<- cbind(X.rct,temp)
    }
  }

  colnames(X.rct) <- paste0("X", 1:P)
  
  #************************ generate rct full data ************************#
  error <- rnorm(n.rct, 0, unex.sig)
  
  Xm <- model.matrix(~., data = X.rct)
  y.rct <- Xm %*% beta + error
  dat.rct <- data.frame(y = y.rct, X.rct)

  #save relevant data
  return(list(dat = dat, dat.rct = dat.rct, 
              p = p, p.sig = P.sig, beta = beta))
  
}

#================================== AUXILIARY PREDICTIONS ====================================#

# epsilon is a vector of the epsilons used for DP private gram matrix. The epsilon
# used for the DP private synthetic data is fixed in 00-data-synth.py
gen_aux_predictions <- function(dat, dat.rct, lambda, 
                                epsilon = c(1,3,6), delta = 10^-5,
                                seed = 2795, mechanism = "gaussian", 
                                syn.dp = T, syn.epsilon = 3){
  
  set.seed(seed)
  
  if(!dir.exists("datasynth")){
    dir.create("datasynth")
  }
  write.csv(dat, file = paste0("datasynth/sim_dat",seed,".csv"))

  #fix weird error with synthpop by changing column name of y column
  syn.dat.feed <- dat
  colnames(syn.dat.feed)[1] <- "X0"
  
#------------------------ synthetic data - non DP ------------------------#
  syn.dat <- synthpop::syn(data = syn.dat.feed, seed=seed, minnumlevels = 2)$syn
  
  lm.syn <- lm(X0~., data = syn.dat)
  beta.syn <- lm.syn$coefficients
  
  print("synthetic data complete")
  
#------------------------ synthetic data - DP  ------------------------#
  if(syn.dp){
    
    py_run_string(paste0('input_dat = "datasynth/sim_dat',seed,'.csv"'))
    py_run_string(paste0('desc_file = "datasynth/description',seed,'.json"'))
    py_run_string(paste0('synth_dat = "datasynth/synth_dat',seed,'.csv"'))
    py_run_string(paste0('syn_dp_ep = ', syn.epsilon))
    
    py_run_file("00-data-synth.py", local = T)
    syn.dat.dp <- read.csv(file = paste0("datasynth/synth_dat",seed,".csv"))[,-1]
    lm.syn.dp <- lm(y~., data = syn.dat.dp)
    
    beta.syn.dp <- lm.syn.dp$coefficients
    
    #clean up
    file.remove(paste0("datasynth/synth_dat",seed,".csv"))
    file.remove(paste0("datasynth/description",seed,".json"))
    file.remove(paste0("datasynth/sim_dat",seed,".csv"))
    
    print("DP synthetic data complete")
  }
  # gram matrix
  lm.g <- lm(y~., data = dat)
  beta.g <- lm.g$coefficients
  
  # true gram matrix
  X <- as.matrix(dat)
  n <- nrow(dat)
  p <- ncol(dat)
  X1 <- cbind(rep(1, nrow(dat)),X)
  G <- t(X1)%*%X1/n
  upper.tri.ind <- upper.tri(G, diag = T)
  
#------------------------ Differentially private gram matrix (g-star) ------------------------#
  
  # determine epsilon budget for each type of statistic
  budget <- c(p,p, (p-1)*p/2)
  budget <- budget/sum(budget)
  names(budget) <- c("mean","sd","corr")
  
  # generate differentially private gram matrices. 
  # start with first epsilon and calculate the noise (which involves a loo calculation)
  dp.ep1 <- dp_private_gram(X, epsilon = epsilon[1], mean.budget = budget["mean"],
                            sd.budget =  budget["sd"], corr.budget =  budget["corr"],
                            delta = delta, mechanism = mechanism, seed = seed)
  
  # now all of hte other noise amounts are just based on the first calculation
  mean.noise <- dp.ep1[["dp.noise"]]$mean*epsilon[1]
  sd.noise <- dp.ep1[["dp.noise"]]$sd*epsilon[1]
  corr.noise <- dp.ep1[["dp.noise"]]$corr[1]*epsilon[1]
  
  #generate the rest of the epsilons
  if(length(epsilon > 1)){
    dp.out <- lapply(epsilon[-1], FUN = function(x){dp_private_gram(X, epsilon = x, mean.budget = budget["mean"],
                              sd.budget =  budget["sd"], corr.budget =  budget["corr"],
                              delta = delta, mechanism = mechanism, seed = seed,
                              mean.noise = mean.noise/x, sd.noise = sd.noise/x, corr.noise = corr.noise/x)})
    names(dp.out) <- as.character(epsilon[-1])
  }
  
  dp.out[[as.character(epsilon[1])]] <- dp.ep1
  
  dp.noise <- list()
  
  G_dp <- dp.out[[as.character(epsilon[1])]]$G
  y_idx <- which(colnames(G_dp) == "y")
  beta.g.dp <- (solve(G_dp[-c(y_idx),-c(y_idx)])%*%G_dp[-c(y_idx),y_idx])
  dp.noise[[paste0("epsilon", epsilon[1])]] <- dp.out[[as.character(epsilon[1])]]$dp.noise
  G.mse <- sum((G_dp*upper.tri.ind - G*upper.tri.ind)^2)/(sum(upper.tri.ind)-1)
  
  if(length(epsilon) >1){
    for(i in 2:length(epsilon)){
      G_dp <- dp.out[[as.character(epsilon[i])]]$G
      y_idx <- which(colnames(G_dp) == "y")
      beta.g.dp <- cbind(beta.g.dp, (solve(G_dp[-c(y_idx),-c(y_idx)])%*%G_dp[-c(y_idx),y_idx]))
      dp.noise[[paste0("epsilon", epsilon[i])]] <- dp.out[[as.character(epsilon[i])]]$dp.noise
      G.mse <- c(G.mse, sum((G_dp*upper.tri.ind - G*upper.tri.ind)^2)/(sum(upper.tri.ind)-1))
    }
  }
  
  colnames(beta.g.dp) <- paste0("epsilon", epsilon)
  names(G.mse) <- paste0("epsilon", epsilon)
  
  print("DP gram matrix complete")
  
#------------------------ Noisy Gram Matrix G-tilde  ------------------------#
  set.seed(seed)
  
  E <- matrix(rnorm(nrow(dat)*ncol(dat),sd = lambda), nrow = nrow(dat))
  X <- cbind(rep(1, nrow(dat)), as.matrix(dat)+E)
  G_tilde <- t(X)%*%X
  
  y_idx <- which(colnames(G_tilde) == "y")
  
  beta.gtilde <- (solve(G_tilde[-c(y_idx),-c(y_idx)])%*%G_tilde[-c(y_idx),y_idx])[,1]
  
#**************************** calculate G mse ********************************# 
  if(syn.dp == T){
    Xsyn <- as.matrix(syn.dat.dp)
    X1 <- cbind(rep(1, nrow(dat)),Xsyn)
    G_dpsyn <- t(X1)%*%X1/n
    
    G.mse <- c(G.mse, sum((G_dpsyn/n*upper.tri.ind - G*upper.tri.ind)^2)/(sum(upper.tri.ind)-1))
    names(G.mse)[length(G.mse)] <- "g.dpsyn"
  }
  
  G.mse <- c(G.mse, sum((G_tilde/n*upper.tri.ind - G*upper.tri.ind)^2)/(sum(upper.tri.ind)-1))
  names(G.mse)[length(G.mse)] <- "g.tilde"
  
#**************************** apply to RCT ********************************#
  
  Mm <- model.matrix(y~., data = dat.rct)
  
  dat.rct$p.syn <- Mm %*% beta.syn
  if(syn.dp){dat.rct$p.syn.dp <- Mm %*% beta.syn.dp}
  dat.rct$p.gram <- Mm %*% beta.g
  dat.rct$p.gtilde <- Mm %*% beta.gtilde
  
  g.dp.preds <- Mm %*% beta.g.dp %>%
    data.frame() %>%
    rename_with(~paste0("p.gdp.",.x))
  
  dat.rct <- cbind(dat.rct, g.dp.preds)
  
  preds <- dat.rct %>%
    select(starts_with("p."))
  
  beta.pred <- data.frame(gram = beta.g,
                          syn = beta.syn,
                          #syn.dp = beta.syn.dp,
                          gtilde = beta.gtilde,
                          beta.g.dp)
  
  return(list(preds = preds, beta.pred = beta.pred, dp.noise=dp.noise, G.mse = G.mse))
  
}

#================================== ESTIMATION SIMULATION ====================================#

################### PARALLEL IMPLEMENTATION ########################


run_sims <- function(n.sim = 1000, tau = 2, dat.rct, preds, seed = 2795,
                          beta, syn.dp = T, epsilon, comparison = "reloop"){
  
  set.seed(seed)
  
  start.time <- Sys.time()
  
  cov_idx <- which(colnames(dat.rct) != "y")
  n_cov <- floor(nrow(dat.rct)/5)
  if(n_cov < length(cov_idx)){
    nonzero.beta<- which(beta[-1] != 0)
    if(length(nonzero.beta) < n_cov){
      samp_cov_idx <- nonzero.beta + 1
    }else{
      samp_cov_idx <- sample(which(beta[-1] != 0), n_cov) + 1
    }
  }else{
    samp_cov_idx <- cov_idx
  }
  
  result.dat <- foreach(i = 1:n.sim, .combine=rbind, .errorhandling = "remove") %dopar% {
    
    #print(paste0("iteration:", i))
    
    it.ob <- c()

    Tr <- sample(c(0,1), nrow(dat.rct), replace = T, prob = c(.5,.5))
    Y <- dat.rct$y + Tr*tau
    reg.dat <- data.frame(Y = Y, Z = Tr, dat.rct[,samp_cov_idx])
    
    it.ob["dif.mean"] <- mean(Y[Tr==1]) - mean(Y[Tr==0])
    
    reg.mod <- lm(Y ~., data = reg.dat)
    it.ob["r2"] <- summary(reg.mod)$r.squared
    it.ob["reg.est"] <- reg.mod$coefficients["Z"]
    
    if(comparison == "reloop"){
      
    Z.reloop <- as.matrix(dat.rct[,cov_idx])
    colnames(Z.reloop) <- NULL
      
    it.ob["reloop.syn"] <- loop(Y = Y, Tr= Tr, Z =  Z.reloop, yhat = preds$p.syn,
                                  pred = reloop)[1]
      
    if(syn.dp){
      it.ob["reloop.syn.dp"] <- loop(Y = Y, Tr= Tr, Z =  Z.reloop, yhat = preds$p.syn.dp,
                                    pred = reloop)[1]
    }
      
    it.ob["reloop.g"] <- loop(Y = Y, Tr= Tr, Z =  Z.reloop, yhat = preds$p.gram,
                              pred = reloop)[1]
      
    it.ob["reloop.gtilde"] <- loop(Y = Y, Tr= Tr, Z =  Z.reloop, yhat = preds$p.gtilde,
                                    pred = reloop)[1]
      
    #allow for any number of epsilons
    for(x in epsilon){
      pred.idx <- which(colnames(preds) == paste0("p.gdp.epsilon", x))
      it.ob[paste0("reloop.gdp",x)] <- loop(Y = Y, Tr= Tr, Z =  Z.reloop, yhat = preds[,pred.idx],
                                            pred = reloop)[1]
    }
      
    }else{
      
      reg.dat <- data.frame(Y = Y, Z = Tr, preds)
      
      it.ob["re.syn"] <- lm(Y ~ Z + p.syn, data = reg.dat)$coefficients["Z"]
      
      if(syn.dp){
        it.ob["re.syn.dp"] <-  lm(Y ~ Z + p.syn.dp, data = reg.dat)$coefficients["Z"]
      }
      
      it.ob["re.g"] <- lm(Y ~ Z + p.gram, data = reg.dat)$coefficients["Z"]
      
      it.ob["re.gtilde"] <- lm(Y ~ Z + p.gtilde, data = reg.dat)$coefficients["Z"]
      
      #allow for any number of epsilons
      for(x in epsilon){
        form <- as.formula(paste0("Y ~ Z + p.gdp.epsilon", x))
        it.ob[paste0("re.gdp",x)] <- lm(form, data = reg.dat)$coefficients["Z"]
      }
      
    }
   
    it.ob
  }
  
  
  end.time <- Sys.time()
  
  elapsed <- end.time - start.time
  
  return(list(result.dat = result.dat, samp_cov_idx = samp_cov_idx, time = elapsed))
  
}

################### SIM WRAPPER ########################

sim_wrapper <- function(n.obs, n.rct,                                   #sample sizes
                        P, prop.norm, prop.sig,                         #n variables and types of vaes
                        alpha, beta.bank, unex.sig,                     #y parameters
                        lambda,                                         #privacy parameters for normal noise
                        mechanism = "gaussian", epsilon, delta = 10^-5, #privacy parameters for DP
                        n.sim.dat, n.sim.treat, tau,                    #specifications for simulations
                        comparison = "reloop",                          #whether to use regression or reloop estimators
                        syn.dp = T                                      #whether or not to run synthetic DP 
){
  
  sim.list <- foreach(i = 1:n.sim.dat) %do% {
    print(paste0("simulated data: ", i))
    
    # generate and save simulated data
    dat.list <- gen_obs_rct_dat(n.obs = n.obs, n.rct = n.rct,
                                P=P, prop.norm = prop.norm, 
                                prop.sig = prop.sig,
                                alpha = alpha,
                                beta.bank = beta.bank,
                                unex.sig = unex.sig,
                                seed = i)
    
    dat <- dat.list$dat
    dat.rct <- dat.list$dat.rct
    
    # calculate auxiliary predictions
    aux.out <- gen_aux_predictions(dat, dat.rct, lambda = lambda, mechanism=mechanism, 
                                   epsilon = epsilon, delta = delta,
                                   seed = i, syn.dp = syn.dp)
    
    # run estimation simulations
    sim.out <- run_sims(n.sim = n.sim.treat, tau = tau, dat.rct = dat.rct, preds = aux.out$preds,
                        beta = dat.list$beta, seed = i, syn.dp = syn.dp, epsilon = epsilon,
                        comparison = comparison)
    
    #save DP noise parameters
    sim.out$dp.noise <- aux.out$dp.noise
    
    #save gram matrix mse
    sim.out$G.mse <- aux.out$G.mse
    
    sim.out$obs.y = dat$y
    sim.out$rct.y = dat.rct$y
    sim.out$obs.r2 <- summary(lm(y~., data = dat))$r.squared
    sim.out$rct.r2 <- summary(lm(y~., data = dat.rct))$r.squared
    
    sim.out
    
  }
  
  return(sim.list)
}



