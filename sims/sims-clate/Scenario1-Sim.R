
library( tidyverse )
# options(list(dplyr.summarise.inform = FALSE))
library(SuperLearner)
library(AER)
library(ranger)
library(gridExtra)

rm(list=ls())


n <- 10000000
ATE <- 0
X1 <- runif(n,-1,1)
X2 <- rbinom(n, 1, 0.3)
U <- runif(n,-1.5,1.5)
X <- cbind.data.frame(X1, X2)
Z <- rbinom(n, 1, plogis(0.4*X[,1] - 0.8*X[,2] + 0.4*I(X[,1] > 0)))

## enforce monotonicity so CLATE = CATE = r.X
p0 <- plogis(-0.3 - 0.4*X[,1] - 0.14*X[,2] - 0.7*I(X[,1] > 0) + 0.7 * U)
p1 <- plogis(-0.3 - 0.4*X[,1] - 0.14*X[,2] + 1.1 - 0.55*X[,1] - 0.7*I(X[,1] > 0) + 0.7 * U)
A0 <- rbinom(n, 1, p0)
A1 <- A0 + (1-A0)*rbinom(n, 1, (p1 - p0) / (1 - p0))
A <- (1 - Z)*A0 + Z*A1

# CATE function
r.X <- (ATE - 4 * X[,1] + 6 * (X[,2] - 0.3) - 4 * (I(X[,1] > 0) - 0.5))
ite.true <- quantile(r.X, probs = seq(.1, .9, by = .1))
clate.true <- (ATE - 4 * .5 + 6 * (1 - 0.3) - 4 * ((.5 > 0) - 0.5))

rm(n, ATE, X1, X2, U, X, Z, A)

### Sim Funcs

make_data <- function(n, ATE){
  X1 <- runif(n,-1,1)
  X2 <- rbinom(n, 1, 0.3)
  U <- runif(n,-1.5,1.5)
  X <- cbind.data.frame(X1, X2)
  Z <- rbinom(n, 1, plogis(0.4*X[,1] - 0.8*X[,2] + 0.4*I(X[,1] > 0)))
  
  ## enforce monotonicity so CLATE = CATE = r.X
  p0 <- plogis(-0.3 - 0.4*X[,1] - 0.14*X[,2] - 0.7*I(X[,1] > 0) + 0.7 * U)
  p1 <- plogis(-0.3 - 0.4*X[,1] - 0.14*X[,2] + 1.1 - 0.55*X[,1] -
                 0.7*I(X[,1] > 0) + 0.7 * U)
  A0 <- rbinom(n, 1, p0)
  A1 <- A0 + (1-A0)*rbinom(n, 1, (p1 - p0) / (1 - p0))
  A <- (1 - Z)*A0 + Z*A1
  
  # CATE function
  r.X <- (ATE - 4 * X[,1] + 6 * (X[,2] - 0.3) - 4 * (I(X[,1] > 0) - 0.5))
  
  # baseline outcome mean
  gamma.X <- (40 - 7*X[,1] - 8*X[,2] + 10*I(X[,1] > 0))
  
  # random noise
  e.Y <- rnorm(n, sd = 0.2)
  
    Y <- r.X * A + gamma.X + 1.5 * U + e.Y
	simDat <- data.frame(cbind(Z, Y, A, X1, X2))
	return(simDat)
}

dat  = make_data( 5000, 0 )

lm(A ~ Z)


eval_data = function(dat, ite.true, clate.true){

	# DRML
	X <- cbind.data.frame(dat$X1, dat$X2)
	  colnames(X) <- c("X1", "X2")
	Y <- dat$Y
	A <- dat$A
	Z <- dat$Z
	
	pi.mod <- SuperLearner(Z, X, family = binomial(),
                         SL.library = c("SL.rpart", "SL.glm", "SL.rpartPrune"))
	lambda1.mod <- SuperLearner(A[Z == 1], X[Z == 1, ], family = binomial(),
	                            SL.library = c("SL.rpart", "SL.glm", "SL.rpartPrune"))
	lambda0.mod <- SuperLearner(A[Z == 0], X[Z == 0, ], family = binomial(),
	                            SL.library = c("SL.rpart", "SL.glm", "SL.rpartPrune"))
	mu1.mod <- SuperLearner(Y[Z == 1], X[Z == 1, ], SL.library = c("SL.rpart", "SL.glm", "SL.rpartPrune"))
	mu0.mod <- SuperLearner(Y[Z == 0], X[Z == 0, ], SL.library = c("SL.rpart", "SL.glm", "SL.rpartPrune"))
  
	pi.pred <- predict.SuperLearner(pi.mod, type="response")$pred
	lambda.pred.1 <- predict.SuperLearner(lambda1.mod, type="response",
	                                      newdata = X)$pred
	lambda.pred.0 <- predict.SuperLearner(lambda0.mod, type="response",
	                                       newdata = X)$pred
	mu.pred.1 <- predict.SuperLearner(mu1.mod,
                                    newdata = X)$pred
	mu.pred.0 <- predict.SuperLearner(mu0.mod,
	                                  newdata = X)$pred

# standard DRML
gamma.hat <- mu.pred.1 - mu.pred.0
delta.hat <- lambda.pred.1 - lambda.pred.0
Psi.hat <- gamma.hat / delta.hat
if.num <- (2*Z - 1) * (Y - ifelse(Z == 1, mu.pred.1, mu.pred.0)) / 
  ifelse(Z == 1, pi.pred, 1 - pi.pred) + gamma.hat
if.denom <-  (2*Z - 1) * (A - ifelse(Z == 1, lambda.pred.1, lambda.pred.0)) / 
  ifelse(Z == 1, pi.pred, 1 - pi.pred) + delta.hat

# using smoothing.spline as the second-step regressor 
#reg.num <- smooth.spline(X$X1, if.num)
#reg.denom <- smooth.spline(X$X1, if.denom)

## CLATE Estimate
cov.df <- data.frame(X)
cov.df$if.a <- if.denom
cov.df$if.y <- if.num
ite.ifa.rf <- ranger(if.a~., cov.df)
ite.ify.rf <- ranger(if.y~., cov.df)
ite.rf <- predictions(ite.ify.rf, data=cov.df)/
    predictions(ite.ifa.rf, data=cov.df)
ite.est <- quantile(ite.rf, probs = seq(.1, .9, by = .1)) 

clate.df <- data.frame(X1 = .5, X2 = 1)
clate.est <- predictions(ite.ify.rf, data=clate.df)/
    predictions(ite.ifa.rf, data=clate.df)    
     
DR = tibble( est = sqrt(mean((ite.est - ite.true)^2)), clate = clate.est - clate.true )
	
# Frauen and Feuerriegel
tau.init <- gamma.hat/delta.hat
num.pseudo <- Y-A*tau.init-mu.pred.0+lambda.pred.0*tau.init
cov.df$pseudo.outcome <-   (2*Z - 1)/delta.hat * num.pseudo /
  ifelse(Z == 1, pi.pred, 1 - pi.pred)+tau.init
  
reg.clate <- ranger(pseudo.outcome ~ ., data=cov.df)
ite.ff <- predictions(reg.clate, data = cov.df)	
clate.ff <- predictions(reg.clate, data = clate.df)	

ite.ff <- quantile(ite.ff, probs = seq(.1, .9, by = .1))

FF = tibble( est = sqrt(mean((ite.ff - ite.true)^2)), clate = clate.ff - clate.true )
          
bind_rows( FF = FF, DR = DR, .id = "method")

}



### Sim Run
sim_reps = 1000
set.seed(234968)

scenarios = expand_grid( n = c( 1000, 2000, 5000, 10000))

run_scenario = function( n ) {
    
    cat( "Simulation: ", n, "\n" )
    
      
    # Run the Simulation              
    reps_qs0 = map_df( 1:sim_reps, function( id ) {
        bdat  = make_data( n, 0 )
        edat = eval_data( bdat, ite.true, clate.true )
        edat$id = id
        edat
    } )
    
    #cat("Sim Done")
    reps_qs0
}


 # timing
t_fexact <- system.time({ 
scenarios$res = pmap( scenarios, run_scenario )

})

# time required for computation, in minutes
t_fexact[['elapsed']]/60

save(scenarios, file="./output/scenario1_results.RData")




