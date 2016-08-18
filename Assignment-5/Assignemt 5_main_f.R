library(maxLik)
library(msm)
library(ggplot2)
library(stargazer)
library(glmnet)


main_f <- function(X, y, beta.1, eps1 = 1e-6, eps2 = 1e-7, maxit = 50) {
  # Fisher's scoring routine for estimation of L3 model (with line search)
  # Input:
  # X = n-by-(r+1) design matrix
  #y - n-by-l vector of success counts
  #beta.1  (r+1)-by-1 vector of starting values for regression est
  #Iteration controlled by:
  #epsl - absolute convergence criterion for beta
  #eps2 - absolute convergence criterion for log-likelihood
  #maxit 8 maximum allowable number of iterations
  
  #Output:
  #out = list containing:
  #beta.MLE = beta MLE
  #NR.hist = iteration history of convergence differences
  #beta.hist = iteration history of beta
  #beta.cov = beta covariance matrix (inverse Fisher's information matrix at MLE)
  #note . convergence note
  
  beta.2 <- rep(-Inf,length(beta.1)) # init beta.2
  diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
  llike.1 <- LL_poisson(y,X,beta.1) # update loglikelihood
  llike.2 <- LL_poisson(y,X,beta.2) # update loglikelihood
  diff.like <- abs(llike.1 - llike.2) # diff
  if (is.nan(diff.like)) { diff.like <- 1e9 }
  i <- 1 # initial iteration index
  alpha.step <- seq(-1, 2, by = 0.1)[-11] # line search step sizes, excluding 0
  NR.hist <- data.frame(i, diff.beta, diff.like, llike.1, step.size = 1) # iteration history
  
  beta.hist <- matrix(beta.1, nrow = 1)
  while ((i <= maxit) & (diff.beta > eps1) & (diff.like > eps2)) {
    i <- i + 1 # increment iteration
    # update beta
    beta.2 <- beta.1 # old guess is current guess
    mu.2 <- exp(X%*%beta.2) # m*p is mean
   # variance matrix
    v.2 <- diag(as.vector(mu.2))
    #v.2<-mean(y)
    score.2 <- t(X) %*% (y - mu.2) # score function
    
    increm <- solve(t(X) %*% v.2 %*% X, score.2) # solve for increment
    
    # line search for improved step size
    llike.alpha.step <- rep(NA, length(alpha.step)) # init llike for line search
      for (i.alpha.step in 1:length(alpha.step)) {
      llike.alpha.step[i.alpha.step] <- LL_poisson(y, X
                                                   , (beta.2 + alpha.step[i.alpha.step] * increm))
      #llike.alpha.step[i.alpha.step] <- LL_poisson(y,X,beta.2)
    }
    
    # step size index for max increase in log-likelihood (if tie, [1] takes first)
    ind.max.alpha.step <- which(llike.alpha.step == max(llike.alpha.step))[1]
    beta.1 <- beta.2 + alpha.step[ind.max.alpha.step] * increm # update beta
    diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
    llike.2 <- llike.1 # age likelihood value
    llike.1 <- LL_poisson(y,X,beta.1) # update loglikelihood
    diff.like <- abs(llike.1 - llike.2) # diff
    
    # iteration history
    NR.hist <- rbind(NR.hist, c(i, diff.beta, diff.like, llike.1, alpha.step[ind.max.alpha.step]))
    beta.hist <- rbind(beta.hist, matrix(beta.1, nrow = 1))
  }
  #prepare output
  out<-list()
  out$beta.MLE<-beta.1
  out$iter<-i-1
  out$NR.hist<-NR.hist
  out$beta.hist<-beta.hist
  v.1<-diag(as.vector(exp(X%*%beta.2)))
  Iinv.1<-solve(t(X)%*%v.1%*%X) #inverse information matrix
  out$beta.cov<-Iinv.1
  if(!(diff.beta>eps1)&!(diff.like>eps2)){
    out$note<-paste("Absolute convergence of", eps1,"for betas
                    and",eps1, "for log likelihood satisfied")
  }
  if(i>maxit){
    out$note<-paste("Exceeded max iterations of ",maxit)
  }
  
  

  if(!(diff.beta>eps1)&!(diff.like>eps2)){
    out$note<-paste("Absolute convergence of", eps1,"for betas
                    and",eps1, "for log likelihood satisfied")
  }
  return(out)
}

LL_poisson<-function(y,X,beta){
  lambda<-mean(y)
  l<-t(y)%*%(X%*%beta) - sum(exp(X%*%beta))
  return(l)
  }

