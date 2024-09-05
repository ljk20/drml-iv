late.comparison <- function(Y, A, Z, X, late){
  naive.est <- mean(Y[A==1])-mean(Y[A==0])
  naive.se <- sqrt(var(Y[A==1])/(sum(A==1)) + var(Y[A==0])/(sum(A==0)))
  
  res.1 <- ivreg(Y~A+.-Z|Z+.-A, data=data.frame(Y=Y, A=A, Z=Z, X))
  tsls.est <- summary(res.1)$coefficient[2,1]
  tsls.se <- summary(res.1)$coefficient[2,2]
  
  res.2 <- ivreg(Y~A+.-Z|Z+.-A, data=data.frame(Y=Y, A=A, Z=Z))
  tsls.un.est <- summary(res.2)$coefficient[2,1]
  tsls.un.se <- summary(res.2)$coefficient[2,2]
  
  late.est <- late$res$est[1]
  late.se <- late$res$se[1]
  
  #res.df <- data.frame(est=naive.est, se=naive.se, 
  #                       estimator="1. Mean difference")
  res.df <- data.frame(est=tsls.un.est, se=tsls.un.se, 
                                     estimator="1. TSLS Unadjusted")
  res.df <- rbind(res.df, data.frame(est=tsls.est, se=tsls.se,
                                     estimator="2. TSLS Adjusted"))                                    
  res.df <- rbind(res.df, data.frame(est=late.est, se=late.se, 
                                     estimator="3. DRML"))                        
  res.df
}

estimate.ite <- function(ate.res, itt.res, X.1){
  cov.df <- data.frame(X.1)
  cov.df$if.a <- ate.res$ifvals[,2]-ate.res$ifvals[,1]
  cov.df$if.y <- itt.res$ifvals[,2]-itt.res$ifvals[,1]
  ite.ifa.rf <- ranger(if.a~., cov.df)
  ite.ify.rf <- ranger(if.y~., cov.df)
  ite.rf <- predictions(ite.ify.rf, data=cov.df)/
    predictions(ite.ifa.rf, data=cov.df)
  
  truncate.range <- function(x, max.val=1, min.val=-1){
    x[x <= min.val] <- min.val; x[x >= max.val] <- max.val
    x
  }
  list(est=truncate.range(ite.rf))
}


iv.profiling.density <- function(V, v0, iv.strength, bw=NULL){
  n <- length(V)
  kern <- function(x){
    dnorm(x)
  }
  if(is.null(bw)){
    bw <- density(V)$bw
  }
  phi <- 1/bw*sapply((v0-V)/bw, kern)
  
  num.if.c <- phi*(iv.strength$ifvals[,2]-iv.strength$ifvals[,1])
  denom.if.c <- iv.strength$ifvals[,2]-iv.strength$ifvals[,1]
  psi <- mean(num.if.c)/mean(denom.if.c)
  psi.if <- (num.if.c-psi*denom.if.c)/mean(denom.if.c)
  psi.se <- sd(psi.if)/sqrt(n)
  ci.ll <- psi-1.96*psi.se; ci.ul <- psi+1.96*psi.se
  res1 <- data.frame(parameter=paste("Complier"), est=psi, se=psi.se,
                     ci.ll,ci.ul)
  
  num.if.at <- phi*iv.strength$ifvals[,1]
  denom.if.at <- iv.strength$ifvals[,1]
  psi <- mean(num.if.at)/mean(denom.if.at)
  psi.if <- (num.if.at-psi*denom.if.at)/mean(denom.if.at)
  psi.se <- sd(psi.if)/sqrt(n)
  ci.ll <- psi-1.96*psi.se; ci.ul <- psi+1.96*psi.se
  res2 <- data.frame(parameter=paste("Always Taker"), est=psi, se=psi.se,
                     ci.ll,ci.ul)
  
  num.if.nt <- phi*(1-iv.strength$ifvals[,2])
  denom.if.nt <- (1-iv.strength$ifvals[,2])
  psi <- mean(num.if.nt)/mean(denom.if.nt)
  psi.if <- (num.if.nt-psi*denom.if.nt)/mean(denom.if.nt)
  psi.se <- sd(psi.if)/sqrt(n)
  ci.ll <- psi-1.96*psi.se; ci.ul <- psi+1.96*psi.se
  res3 <- data.frame(parameter=paste("Never Taker"), est=psi, se=psi.se,
                     ci.ll,ci.ul)
  
  res <- rbind(res1,res2,res3); rownames(res) <- NULL
  return(res)
}

iv.profiling.prob <- function(V, v0, iv.strength){
  phi <- as.numeric(V==v0)
  n <- length(V)
  num.if <- cbind(phi*(iv.strength$ifvals[,2]-iv.strength$ifvals[,1]),
                  phi*(iv.strength$ifvals[,1]),
                  phi*(1-iv.strength$ifvals[,2]))
  denom.if <- cbind(iv.strength$ifvals[,2]-iv.strength$ifvals[,1],
                    iv.strength$ifvals[,1],
                    1-iv.strength$ifvals[,2])
  psi <- apply(num.if, 2, mean)/apply(denom.if, 2, mean)
  psi.if <- t(t(num.if-t(psi*t(denom.if)))/apply(denom.if, 2, mean))
  psi.se <- apply(psi.if,2,sd)/sqrt(n)
  ci.ll <- psi-1.96*psi.se; ci.ul <- psi+1.96*psi.se
  data.frame(parameter=c("Complier", "Always Taker", "Never Taker"),
             est=psi, se=psi.se, ci.ll, ci.ul)
}
