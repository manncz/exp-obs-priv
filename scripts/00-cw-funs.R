# script: 00-cw-funs.R
# author: Charlotte Mann adapted from genRCT package (https://github.com/idasomm/genRCT)
# date: Dec 13, 2022
# purpose: Functions for calibration weighted estimator


# load required packages
library(nleqslv)

### from utilities.R of genRCT package
lamFun <- function(lam, moments, moments.bar) { ## vector lam and x
  qi <- (exp(moments %*% lam) / sum(exp(moments %*% lam)))[, 1]
  colSums(qi  * moments) - moments.bar
  
  
}
####


### adapted from genRCT.estimators.R of genRCT package

cw_estimator <- function(X.trial, Y.trial, Tr, mean.aux,
                         adj = F, seed = NULL, family = "gaussian"){
  
    if (!is.null(seed)) set.seed(seed)
    
    dat.trial <- data.frame(Tr = Tr, X.trial)
    dat.trial$Y <- Y.trial
    h.trial <- cbind(1, as.matrix(X.trial))
    
    ## Estimate adjustment models for ACW-t adjusted estimator
    if(adj == T){
      
      if (family == 'gaussian') {
        fit.t0 <- lm(Y ~.-Tr, data = dat.trial[dat.trial$Tr == 0,])
        fit.t1 <- lm(Y ~.-Tr, data = dat.trial[dat.trial$Tr == 1,])
        
        beta.t0 <- as.numeric(coef(fit.t0))
        beta.t1 <- as.numeric(coef(fit.t1))
        
        dat.trial$Y0.t <- h.trial %*% beta.t0
        dat.trial$Y1.t <- h.trial %*% beta.t1
        
        mu0.aux <- c(1,mean.aux) %*% beta.t0
        mu1.aux <- c(1,mean.aux) %*% beta.t1
        
      } else if (family == 'binomial') {
        fit.t0 <- glm(Y ~.-Tr, data = dat.trial[dat.trial$Tr == 0], family = binomial("logit"))
        fit.t1 <- glm(Y ~.-Tr, data = dat.trial[dat.trial$Tr == 1], family = binomial("logit"))
        
        beta.t0 <- as.numeric(coef(fit.t0))
        beta.t1 <- as.numeric(coef(fit.t1))
        
        dat.trial$Y0.t <- expit(h.trial %*% beta.t0)
        dat.trial$Y1.t <- expit(h.trial %*% beta.t1)
        
        mu0.aux <- expit(c(1,mean.aux) %*% beta.t0)
        mu1.aux <- expit(c(1,mean.aux) %*% beta.t1)
      }
    }
    
    ## estimate weights
    cnt <- 0
    while (cnt <= 50) {
      cnt <- cnt + 1
      lam.hat <- searchZeros(matrix(rnorm(length(mean.aux) * 10, 0, 0.25), nrow = 10), lamFun, moments = X.trial, moments.bar =mean.aux)$x[1,]
      if (!is.null(lam.hat)) break
    }
    
    if (is.null(lam.hat)) {
      warning('No lam.hat solutions')
      lam.hat <- rep(NA, ncol(X.trial))
    }
    q.score <- exp(X.trial %*% lam.hat) / sum(exp(X.trial %*% lam.hat))
    dat.trial$q <- q.score
    
    ## calculate CW estimator
    tau.cw <- dat.trial %>%
      summarize(est = sum(2*q*(Tr*Y - (1-Tr)*Y))/sum(q))
    tau.cw <- tau.cw$est
    
    ## calculate ACW-t estimator
    if(adj == T){
      
      tau.acw.t <- dat.trial %>%
        summarize(est = sum(2*q*(Tr*(Y-Y1.t) - (1-Tr)*(Y-Y0.t)))/sum(q))
      
      tau.acw.t <- tau.acw.t$est + mu1.aux - mu0.aux
      
        tau.cw <- c(tau.cw,tau.acw.t)
        names(tau.cw) <- c("CW", "ACW-t")
    }
    
    return(tau.cw)
    
}

## adapted from genRCT.estimators.R from genRCT packages

cw_bootstrap_inference <- function(X.trial, Y.trial, Tr, mean.aux,
                             adj = F, seed = NULL, family = "gaussian",
                             nboot = 100, level = .05){
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(X.trial)
  
  boot.dist <- foreach(i = 1:nboot, .combine = rbind) %dopar% {
    
    if(i %% 10 == 0) print(paste("completing bootstrap", i, "of", nboot))
    
    #get bootstrap sample
    bsamp <- sample(n, replace = TRUE)
    
    X.trialb <- X.trial[bsamp,]
    Y.trialb <- Y.trial[bsamp]
    Trb <- Tr[bsamp]
  
    boot.est <- cw_estimator(X.trial = X.trialb, Y.trial = Y.trialb, Tr = Trb, mean.aux = mean.aux, adj=adj)
  
    boot.est
    
  }
  
  se <- apply(boot.dist, 2, sd, na.rm = T)
  
  ci <- apply(boot.dist, 2, function(x) quantile(x, probs = c(level/2, 1-(level/2)), na.rm = T))
  rownames(ci) <- c("lb", "ub")
  
  inf.res <- rbind(se, ci)
  
  return(inf.res)
  
}

####
