
library(tidyverse)
library(npcausal)
library(ranger)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(foreign)

rm(list=ls())

data <- read.dta("iv_match.dta")

data.1 <- data %>% filter(shafi17==1)
data.1$outcome <- as.numeric(data.1$died + data.1$prolonged > 0)

rm(data)
nrow(data.1)

vars <- c("comorb", "angus", "age","ynel1","ynel2","ynel3","ynel4", "ynel5",
            "ynel6","ynel7","ynel8","ynel9","ynel10","ynel11","ynel12","ynel13",
            "ynel14","ynel15","ynel16","ynel17","ynel18","ynel19","ynel20","ynel21",
            "ynel22","ynel23","ynel24","ynel25","ynel26","ynel27","ynel28",
            "ynel29","ynel30", "ynel31", "p_cat1", "p_cat2", "p_cat3","p_cat4",
            "p_cat5", "female", "hisp", "afam", "other", "disability_bin")
            
  Y <- data.1$outcome; A <- data.1$surg; Z <- data.1$pref_bin
  X <- as.matrix(data.1[,vars])
  hosp.mat <- matrix(model.matrix(~factor(data.1$hospid)-1), nrow=length(Y))
  X <- cbind(X, hosp.mat)


source('functions.R')
	
sl.list =c("SL.gam","SL.glm","SL.ranger")
nsplits = 5

t_fexact <- system.time({ 
	
y.on.z <- ate(Y, Z, X, nsplits = nsplits, sl.lib = sl.list)
a.on.z <- ate(A, Z, X, nsplits = nsplits, sl.lib = sl.list)
late.five <- ivlate(Y, A, Z, X, nsplits = nsplits, sl.lib = sl.list)
  
})

# time required for computation, in minutes
t_fexact[['elapsed']]/60

# plot distribution of ITE
ite.res <- estimate.ite(ate.res = a.on.z, itt.res = y.on.z, X.1 = X)
 
save(y.on.z, a.on.z, late.five, ite.res, file = "nuis.out")

