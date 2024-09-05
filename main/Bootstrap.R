

library(tidyverse)
library(npcausal)
library(ranger)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(AER)
library(foreign)
library(gridExtra)
library(gam)

rm(list=ls())


data <- read.dta("iv_match.dta")

data.1 <- data %>% filter(shafi17==1)
data.1$outcome <- as.numeric(data.1$died + data.1$prolonged > 0)

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
load("nuis.out")


# estimate bootstrap quantile
resample <- function(x) {
  return(sample(x, size=length(x), replace=TRUE))
}

resample.data.frame <- function(df) {
  return(df[resample(1:nrow(df)), ])
}

# bootstrap sample 
age.range <- 20:80 # grid for estimating cate
comorb.range <- 0:10  # grid for estimating cate
boot.sample <- 1000 # has to be much larger but takes time (I used 5000 for the report) 

cov.df <- data.frame(X)
cov.df$if.a <- a.on.z$ifvals[,2]-a.on.z$ifvals[,1]
cov.df$if.y <- y.on.z$ifvals[,2]-y.on.z$ifvals[,1]
gam.ifa.comorb.sepsis <- gam(if.a~s(comorb), cov.df, family=gaussian)
gam.ify.comorb.sepsis <- gam(if.y~s(comorb), cov.df, family=gaussian)

pred.df.1 <- data.frame(angus=c(rep(0, length(comorb.range)), rep(1, length(comorb.range))), 
                      comorb=rep(comorb.range), 2)


cate.comorb.sepsis <- predict(gam.ify.comorb.sepsis, pred.df.1)/
  predict(gam.ifa.comorb.sepsis, pred.df.1)

## Bootstrap
one.sim.sepsis.comorb <- function(){
  df <- resample.data.frame(cov.df[, c("age", "angus", "comorb","if.a", "if.y")])
  gam.ifa.comorb.sepsis <- gam(if.a~s(comorb), df, family=gaussian)
  gam.ify.comorb.sepsis <- gam(if.y~s(comorb), df, family=gaussian)
  predict(gam.ify.comorb.sepsis, pred.df.1)/predict(gam.ifa.comorb.sepsis, pred.df.1)
}

boots.sample.comorb <- replicate(boot.sample, one.sim.sepsis.comorb())
boots.CI.comorb <- t(apply(boots.sample.comorb, 1, quantile, c(0.025, 0.975)))

## Final Results
pred.df.1$cate <- cate.comorb.sepsis
pred.df.1$upper.ci <- boots.CI.comorb[,2]
pred.df.1$lower.ci <- boots.CI.comorb[,1]

  
save(pred.df.1, file = "boot.out")