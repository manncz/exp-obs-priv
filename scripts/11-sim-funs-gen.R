# script: 11-sim-funs-gen.R
# author: Charlotte Mann
# date: Nov 17, 2022
# purpose: Functions to support simulation demonstrating esitmating PATE
#           when there is a heterogenous treatment effect and covariate shift
#           between aux data and rct data

#==========================================================================================#

#================================== GENERATE SIMULATED DATA ====================================#
gen_obs_rct_dat_shift <- function(n.obs = 1000, n.rct = 50,
                            P = 10, prop.sig = .8, prop.norm = .6,
                            beta=NULL, select.beta=NULL,
                            alpha = 10, beta.bank = 10, unex.sig=2,
                            select.alpha = -2.5, select.gamma = .5,
                            select.beta.bank = 0, select.prop.sig = .5,
                            seed = 12345){
  set.seed(seed)
  
  P.sig <- floor(P*prop.sig) #number of non-null covariates
  P.norm <- floor(P*prop.norm) #number of continuous covariates
  P.cat <- P-P.norm #number of categorical covariates
  
  #if a beta vector isn't provided, generate a random beta
  if(is.null(beta)){
    #generate beta
    beta <- rep(0, P)
    beta[sample(1:(P), P.sig)] <- sqrt(beta.bank/P.sig)
    #last coefficient is 0 since S only interacts w outcome if treated
    beta <- c(alpha, beta, 0)
  }

  error <- rnorm(n.obs, 0, unex.sig)
  
  #******************** generate obs covariate data *********************#
  
  # generate continuous covariates
  X <- data.frame(matrix(rnorm(n.obs*P.norm,1,1), ncol=P.norm))
  
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
  
  #also generate extra variable S, which interacts with treatment
  # so in the observational data, it does not effect y
  X$S <- rnorm(n.obs,1,1)
  
  #************************ generate obs full data ************************#
  
  Xm <- model.matrix(~., data = X)
  y <- Xm %*% beta + error
  
  dat <- data.frame(y, X)
  
  #******************** generate rct covariate data *********************#
  
  # if a selection beta vector isn't provided, generate a random selection beta vector
  if(is.null(select.beta)){
    #generate select.beta
    P.sig <- floor(P*select.prop.sig) #number of non-null covariates
    select.beta <- rep(0, P)
    select.beta[sample(1:(P), P.sig)] <- select.beta.bank/P.sig
    select.beta <- c(select.alpha, select.beta, select.gamma)
  }
 
    
  #determine how many samples to take from to get approximately n.rct in the actual RCT data
  inflate.factor <- floor(1/expit(sum(select.beta)))
    
  # generate continuous covariates
  X.rct <- data.frame(matrix(rnorm(n.rct*inflate.factor*P.norm,1,1), ncol=P.norm))
  
  # generate categorical covariates
  if(P.cat != 0){
    for(i in 1:P.cat){
      temp <- sample(c(0,1), size = n.rct, replace = T, prob = c(1-p[i], p[i]))
      X.rct<- cbind(X.rct,temp)
    }
  }
  
  colnames(X.rct) <- paste0("X", 1:P)
  X.rct$S <- rnorm(n.rct*inflate.factor,1,1)
  
  Xm <- model.matrix(~., data = X.rct)
  
  p.select <- expit(Xm %*% select.beta)
  select <- rbinom(length(p.select),1,p.select)
  
  # now only select participants based on the selection vector
  X.rct <- X.rct[select ==1,]
  
  #************************ generate rct full data ************************#
  error <- rnorm(nrow(X.rct), 0, unex.sig)
  
  Xm <- model.matrix(~., data = X.rct)
  y.rct <- Xm %*% beta + error
  dat.rct <- data.frame(y = y.rct, X.rct)
  
  #save relevant data
  return(list(dat = dat, dat.rct = dat.rct, 
              p = p, beta = beta,
              select.beta = select.beta))
  
}

#================================== AUXILIARY MEANS ====================================#

gen_aux_means <- function(dat, dat.rct, 
                          lambda, epsilon = c(1,3,6), delta = 10^-5,
                          mechanism = "gaussian", syn.dp = T, syn.epsilon=3,
                          dp.upbound = NULL, dp.lbound = NULL,
                          seed = 2795){
  
  set.seed(seed)
  
  if(!dir.exists("datasynth")){
    dir.create("datasynth")
  }
  write.csv(dat, file = paste0("datasynth/sim_dat",seed,".csv"))
  
  #------------------------       true data        ------------------------#
  
  mu.aux <- colMeans(dat)
  
  #------------------------ synthetic data - non DP ------------------------#
  #fix weird error with synthpop by changing column name of y column
  syn.dat.feed <- dat
  colnames(syn.dat.feed)[1] <- "X0"

  syn.dat <- synthpop::syn(data = syn.dat.feed, seed=seed, minnumlevels = 2)$syn
  
  mu.syn <- colMeans(syn.dat)
  names(mu.syn)[1] <- "y"
  
  #------------------------ synthetic data - DP  ------------------------#
  if(syn.dp){
    
    py_run_string(paste0('input_dat = "datasynth/sim_dat',seed,'.csv"'))
    py_run_string(paste0('desc_file = "datasynth/description',seed,'.json"'))
    py_run_string(paste0('synth_dat = "datasynth/synth_dat',seed,'.csv"'))
    py_run_string(paste0('syn_dp_ep = ', syn.epsilon))
    
    py_run_file("00-data-synth.py", local = T)
    syn.dat.dp <- read.csv(file = paste0("datasynth/synth_dat",seed,".csv"))[,-1]
    
    mu.syn.dp <- colMeans(syn.dat.dp)
    
    #clean up
    file.remove(paste0("datasynth/synth_dat",seed,".csv"))
    file.remove(paste0("datasynth/description",seed,".json"))
    file.remove(paste0("datasynth/sim_dat",seed,".csv"))
  }
  
  print("synthetic data complete")
  
  #------------------------ Differentially private G ------------------------#
  X <- as.matrix(dat)
  n <- nrow(dat)
  p <- ncol(dat)
  
  budget <- c(p, p, (p-1)*p/2)
  budget <- budget/sum(budget)
  names(budget) <- c("mean","e2","e12")
  
  mean.epsilon <- epsilon[1]*budget["mean"]
  
  if(!is.null(dp.upbound)){
    
    mu.dp <- foreach(i = 1:length(epsilon), .combine = cbind) %do% {
      
      dp.g <- dp_private_gram(X, epsilon = epsilon[i], 
                                mean.budget = budget["mean"], e2.budget =  budget["e2"], e12.budget = budget["e12"],
                                lbound = dp.lbound, ubound = dp.upbound,
                                delta = delta, mechanism = mechanism, seed = seed)
      
      it.mu <- dp.g[["G"]][1,-1]
      
    }
    
    mean.noise = NULL
    
  }else{
    if(mechanism == "gaussian"){
      mean.delta <- delta*budget
      mean.noise <- calc_dp_err_gaus(X, epsilon = mean.epsilon, delta = mean.delta,
                                 statistic = "mean")
    }else{
      mean.noise <- calc_dp_err_lap(X, epsilon = mean.epsilon,statistic = "mean")
    }
  
    mu.dp <- foreach(i = 1:length(epsilon),.combine =cbind) %do% {
      set.seed(seed)
      
      mean.noise <- mean.noise*epsilon[1]/epsilon[i]
     
      if(mechanism=="gaussian"){
        mean.noise.vec <- sapply(mean.noise, function(x){rnorm(1, mean = 0, sd = x)})
      }else{
        mean.noise.vec <- sapply(mean.noise, function(x){rlaplace(p, scale = x)})
      }
      
      it.mu <- mu.aux + mean.noise.vec
      
      it.mu
    }
    
    names(mean.noise) <- colnames(dat)
  }
  
  colnames(mu.dp) <- paste0("m.epsilon", epsilon)
  
  print("DP means complete")
  
  #------------------------ Noisy Gram Matrix G-tilde  ------------------------#
  
  E <- matrix(rnorm(nrow(dat)*ncol(dat),sd = lambda), nrow = nrow(dat))
  mu.gtilde <- colMeans(X+E)
  
  #**************************** bring means together ********************************# 
  aux.mus <- data.frame(m.aux = mu.aux,
                        m.gtilde = mu.gtilde,
                        m.syn = mu.syn,
                        mu.dp
                        )
  if(syn.dp){aux.mus$m.syn.dp <- mu.syn.dp}
  
  #**************************** return dat ********************************# 
  return(list(aux.means = aux.mus,
              mean.noise = mean.noise))
  
}

#================================== ESTIMATION SIMULATION ====================================#


run_sims_gen <- function(n.sim, tau, dat.rct, aux.means, beta,
                         level, cw.inference, nboot.inf,  
                         seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  cov_idx <- which(colnames(dat.rct) != "y")
  
  #for regression estimator, only include a number of covariates that is at most
  #1/5 the sample size
  n_cov <- floor(nrow(dat.rct)/5)
  if(n_cov < length(cov_idx)){
    nonzero.beta<- which(beta[-1] != 0)
    if(length(nonzero.beta) < n_cov){
      samp_cov_idx <- c(nonzero.beta + 1,  which(colnames(dat.rct) == "S"))
    }else{
      samp_cov_idx <- sample(which(beta[-1] != 0), n_cov) + 1
      #make sure S is included
      samp_cov_idx <- unique(c(samp_cov_idx, which(colnames(dat.rct) == "S")))
    }
  }else{
    samp_cov_idx <- cov_idx
  }
  
  result.dat <- foreach(i = 1:n.sim, .combine=rbind) %do% {
    
    Tr <- sample(c(0,1), nrow(dat.rct), replace = T, prob = c(.5,.5))
    Y <- dat.rct$y + Tr*tau*dat.rct$S
    X <- dat.rct %>%
      select(-y) %>%
      as.matrix()
    
    it.ob <- matrix(rep(NA, 4),nrow=2)
    colnames(it.ob) <- c("dif.mean", "reg.est")
    rownames(it.ob) <- c("est", "coverage")
    z <- qnorm(p = 1-(level/2))
    
    #difference in means
    it.ob["est","dif.mean"] <- mean(Y[Tr==1]) - mean(Y[Tr==0])
    
    se <- sqrt(var(Y[Tr==1])/sum(Tr==1) + var(Y[Tr==0])/sum(Tr==0))
    it.ob["coverage","dif.mean"] <- as.numeric(tau >= it.ob["est","dif.mean"]-(z*se) & tau <= it.ob["est","dif.mean"]+(z*se))
    
    #regression estimator
    reg.dat <- data.frame(Y = Y, Z = Tr, dat.rct[,samp_cov_idx])
    
    reg.mod <- lm(Y ~., data = reg.dat)
    it.ob["est", "reg.est"] <- reg.mod$coefficients["Z"]
    it.ob["coverage", "reg.est"] <- as.numeric(tau >= confint(reg.mod, "Z")[1] & tau <= confint(reg.mod, "Z")[2])
    
    #calibaration weighted esimators and bootstrap inference
    n.do<- ncol(aux.means)
    aux.res <- foreach(j = 1:n.do, .combine = cbind) %do% {
      est <- cw_estimator(X.trial = X, Y.trial = Y, Tr = Tr,
                           mean.aux = aux.means[-1,j], adj=T)
      
      if(cw.inference){
        ci <- cw_bootstrap_inference(X.trial = X, Y.trial = Y, Tr = Tr, mean.aux = aux.means[-1,j],
                                     adj = T, nboot = nboot.inf, level = level)
        
        coverage <- as.numeric(tau >= ci["lb",] & tau <= ci["ub",])
      }else{
        coverage <- c(NA, NA)  
        names(coverage) <- names(est)
      }
      
      print(paste("cw", j, "complete"))
      
      res <- rbind(est, coverage)
      
      colnames(res) <- paste(str_replace_all(str_to_lower(colnames(res)), "-", ""),
                            str_replace(colnames(aux.means)[j], "m.", ""), sep = ".")
      res
    }
    
    cbind(it.ob, aux.res)
  
  }
  
  result.dat <- result.dat %>%
    data.frame %>%
    mutate(sate = mean(dat.rct$S*tau))
  
  return(result.dat)
  
}

#================================== FULL SIMULATION WRAPPER ====================================#

sim_wrapper_gen <- function(n.obs, n.rct,                               #sample sizes
                        P, prop.norm, prop.sig,                         #n variables and types of vars
                        alpha, beta.bank, unex.sig,                     #y parameters
                        select.alpha, select.gamma,                     #rct selection parameters
                        select.beta.bank, select.prop.sig,
                        lambda,                                         #privacy parameters for normal noise
                        mechanism = "gaussian", epsilon, delta = 10^-5, #privacy parameters for DP
                        lbounddp = NULL, ubounddp = NULL,                   #bounds on columns for DP
                        n.sim.dat, n.sim.treat, tau,                    #specifications for simulations
                        syn.dp = T,  syn.epsilon =3,                    #whether or not to run synthetic DP 
                        level = .05, cw.inference = F, nboot.inf = NULL, #parameters for inference
                        seed = 123){
  
  set.seed(seed)
  
  #generate beta
  P.sig <- floor(P*prop.sig) #number of non-null covariates
  beta <- rep(0, P)
  beta[sample(1:(P), P.sig)] <- sqrt(beta.bank/P.sig)
  #last coefficient is 0 since S only interacts w outcome if treated
  beta <- c(alpha, beta, 0)
  
  #generate select.beta
  P.sig <- floor(P*select.prop.sig) #number of non-null covariates
  select.beta <- rep(0, P)
  select.beta[sample(1:(P), P.sig)] <- select.beta.bank/P.sig
  select.beta <- c(select.alpha, select.beta, select.gamma)
  
  sim.list <- foreach(i = 1:n.sim.dat, .errorhandling = "remove") %do% {
    print(paste0("simulated data: ", i))
    
    dat.list <- gen_obs_rct_dat_shift(n.obs = n.obs, n.rct = n.rct,
                                      P = P, prop.norm = prop.norm,
                                      beta=beta, select.beta=select.beta,
                                      unex.sig = unex.sig,
                                      seed = seed + i)
    
    aux.out <- gen_aux_means(dat=dat.list$dat, dat.rct=dat.list$dat.rct, lambda=lambda,
                             epsilon = epsilon, delta = delta, mechanism = mechanism, 
                             dp.lbound = lbounddp, dp.upbound = ubounddp,
                             syn.dp = syn.dp, syn.epsilon = syn.epsilon,
                             seed = seed + i)
    
    print(paste0("auxiliary means ", i, "done"))
    
    #since we will be only running one treatment assignment vector for each dataset, set 
    #the seed outside of the function, so the treatment assignment is actually random
    set.seed(seed + i)
    
    sim.res <- run_sims_gen(n.sim = n.sim.treat, tau = tau, dat.rct=dat.list$dat.rct,
                            aux.means=aux.out$aux.means, beta=beta,
                            cw.inference = cw.inference, nboot.inf = nboot.inf, level = level)
    
    sim.out <- list(result.dat = sim.res)
    
    sim.out$dp.noise <- aux.out$mean.noise
    sim.out$aux.means <- aux.out$aux.means %>%
      data.frame() %>%
      mutate(m.true = c(sum(beta), rep(1, P+1)))
    
    sim.out
    
  }
  
  sim.list$beta <- beta
  sim.list$select.beta <- select.beta

  return(sim.list)
  
}
