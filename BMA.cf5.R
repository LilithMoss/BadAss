###############################
# Admixture Analysis using BMA
###############################
library(Rmpfr)

BMA.cf <- function(dat=NULL,alpha.SE=0.05,beta.SE=0.05,cc=0.5,co=0.5,stand=F,control.sd=NULL){
  #Run simple CO and CC models on standardized deltaA
  reg <- as.list(rep(0, numModels))
  Y<-dat$Y; Z<-dat$Z; Q<-dat$Q; L<-dat$L;
  #Model weights
  pmw <- c(co,cc)
  n <- length(Y)
  #Create SE value - 1/25/2018 if NULL
  if(is.null(control.sd)){
    con <- dat[dat$Y==0,] #Controls only in dat
    control.sd <- sd(as.numeric(con$Q-(con$L/2)))
  }
  alpha.SE <- alpha.SE*control.sd
  beta.SE <- beta.SE*control.sd
  #Create outcome variable
  if(stand==F){
    P <- (as.numeric(Q-(L/2))) #produces more interpretable results
  } else if(stand==T){
    P <- (as.numeric(Q-(L/2)))/(sd(as.numeric(Q-(L/2))))
  }
  #Run standard regression models
  #Case-Only
  reg[[1]] <- lm(P ~ -1 + Y) #Case-only
  #Case-Control
  reg[[2]] <- lm(P ~ 1 + Y)  #Case - Control
  ####################
  #Run BMA analysis
  ####################
  #Calculate BMA Result
  betas <- (unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Estimate"] })))
  # betas.se <- sdA*(unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] }))) #11.8.2017
  betas.se <- (unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] })))
  #Set alpha.SE to what it actually is
  # alpha.SE <- 100*summary(reg[[2]])$coefficients[1,2]
  #Prior Covariance matrix for the case-control model
  V.cc <- matrix(c((as.numeric(alpha.SE)^2/sigma.1),0,0,as.numeric(beta.SE)^2/sigma.1),ncol=2,byrow=T)
  #Prior Variance of Case-only matrix using SE(Beta) for the intercept
  cov1 <- sigma.1*matrix(as.numeric(beta.SE)^2/sigma.1) 
  #Prior Variance of Case-control matrix
  cov2 <- sigma.1*V.cc 
  #Prior covariance matricies for both CO and CC models
  Cov_0 <- list(cov1,cov2)
  #Invert the Prior covariance matrix to get prior precision matrix
  lambda_0 <- list(solve(Cov_0[[1]]),solve(Cov_0[[2]]) )
  #Specify Design Matrix (list of design matrices)
  X <- list( matrix(Y),matrix(c(rep(1,length(Y)),Y),ncol=2) )
  #Specify Posterior Precision Matrix from Prior precision for each model
  lambda_n <- list(t(X[[1]])%*%X[[1]]+lambda_0[[1]],t(X[[2]])%*%X[[2]]+lambda_0[[2]])
  #Specify prior mean vector for Betas and Alpha
  mu_0 <- list(matrix(0),matrix(c(0,0)))
  #Specify posterior mean vector for betas
  mu_n <- list( solve(t(X[[1]])%*%X[[1]]+lambda_0[[1]]) * ( (t(X[[1]])%*%X[[1]]*betas[1]) + (lambda_0[[1]]*mu_0[[1]]) ),
                solve(t(X[[2]])%*%X[[2]]+lambda_0[[2]]) %*% ( (t(X[[2]])%*%X[[2]] %*% matrix(reg[[2]]$coefficients)) + (lambda_0[[2]] %*% mu_0[[2]]) ) )
  b0 = 1
  a0 = 1/(sigma.1) + 1
  #Specify Posterior hyperparameters for sigma^2
  an <- a0+(n/2)
  # bn <- list( b0+(1/2)*(t(Q)%*%Q + t(mu_0[[1]])%*%lambda_0[[1]]%*%mu_0[[1]] - t(mu_n[[1]])%*%lambda_n[[1]]%*%mu_n[[1]]),
  #             b0+(1/2)*(t(Q)%*%Q + t(mu_0[[2]])%*%lambda_0[[2]]%*%mu_0[[2]] - t(mu_n[[2]])%*%lambda_n[[2]]%*%mu_n[[2]]) )
  bn <- list( b0+(1/2)*(t(P)%*%P + t(mu_0[[1]])%*%lambda_0[[1]]%*%mu_0[[1]] - t(mu_n[[1]])%*%lambda_n[[1]]%*%mu_n[[1]]),
              b0+(1/2)*(t(P)%*%P + t(mu_0[[2]])%*%lambda_0[[2]]%*%mu_0[[2]] - t(mu_n[[2]])%*%lambda_n[[2]]%*%mu_n[[2]]) )
  
  #Calculate large values using multiple precision package (Rmpfr)
  lterm1 <- exp(as(((-n/2)*log(2*pi)),"mpfr")) #1/(2pi)^(n/2)
  lterm2 <- list( exp(as(an*log(bn[[1]]),"mpfr")), exp(as(an*log(bn[[2]]),"mpfr")) )
  lterm3 <- gamma(as(an,"mpfr"))
  lterm4 <- gamma(as(a0,"mpfr"))
  #Calculate Marginal Likelihood
  PrDGivenM1 <-lterm1*sqrt(det(lambda_0[[1]])/det(lambda_n[[1]]))*((b0^a0)/lterm2[[1]])*(lterm3/lterm4) 
  PrDGivenM2 <-lterm1*sqrt(det(lambda_0[[2]])/det(lambda_n[[2]]))*((b0^a0)/lterm2[[2]])*(lterm3/lterm4) 
  PrDGivenM <- c(PrDGivenM1,PrDGivenM2)
  #Calculate Posterior Model Probabilities
  PrMGivenD.new <- (PrDGivenM*pmw)/sum( PrDGivenM*pmw )
  post.beta <- sum(betas*PrMGivenD.new)
  post.se <- sqrt(sum(((betas.se^2)+(betas^2))*PrMGivenD.new) - (post.beta^2))
  z.score <- post.beta/post.se
  p.value <- 2*pnorm(-abs(z.score))
  #BMA.cf <- as.numeric( c(post.beta, post.se, z.score, p.value))
  BMA.cf.res <- as.numeric( c(post.beta, post.se, z.score, p.value,
                          as.numeric(PrMGivenD.new[1]),as.numeric(PrMGivenD.new[2])))
  names(BMA.cf.res) <- c("post.beta", "post.se", "z.score", "p.value",
                     "PrDGivenCO","PrDGivenCC")
  return(BMA.cf.res)
}


