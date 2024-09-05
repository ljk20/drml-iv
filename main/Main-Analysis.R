
library(npcausal)
library(ranger)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(AER)
library(gridExtra)
library(gam)


library(tidyverse)
library(foreign)

rm(list=ls())

data <- read.dta("iv_match.dta")

data.1 <- data %>% filter(shafi17==1)
data.1$outcome <- as.numeric(data.1$died + data.1$prolonged > 0)

df.tto <- data %>% group_by(researchid) %>%
           summarize(pref = mean(pref)*100)
           
summary(df.tto$pref)           
                   
ggplot(df.tto, aes(x=pref)) + 
        geom_histogram(fill="#69b3a2", color="#e9ecef") +
        xlab("TTO: Percentage of Operations Peformed") +
        ylab("") +
        theme_minimal() +
        theme(axis.title = element_text(size = 15)) +
        theme(axis.text = element_text(size = 15)) 

          
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

## IV Strength
names <- c("pref_bin","comorb", "angus", "age","ynel1","ynel2","ynel3","ynel4", "ynel5",
            "ynel6","ynel7","ynel8","ynel9","ynel10","ynel11","ynel12","ynel13",
            "ynel14","ynel15","ynel16","ynel17","ynel18","ynel19","ynel20","ynel21",
            "ynel22","ynel23","ynel24","ynel25","ynel26","ynel27","ynel28",
            "ynel29","ynel30", "p_cat1", "p_cat2", "p_cat3","p_cat4",
            "female", "hisp", "afam", "other", "disability_bin")
            
frmla <- as.formula(paste("surg ~", paste(names, collapse = "+"), sep=""))
fs <- lm(frmla, data=data.1)

names <- c("comorb", "angus", "age","ynel1","ynel2","ynel3","ynel4", "ynel5",
            "ynel6","ynel7","ynel8","ynel9","ynel10","ynel11","ynel12","ynel13",
            "ynel14","ynel15","ynel16","ynel17","ynel18","ynel19","ynel20","ynel21",
            "ynel22","ynel23","ynel24","ynel25","ynel26","ynel27","ynel28",
            "ynel29","ynel30", "p_cat1", "p_cat2", "p_cat3","p_cat4",
            "female", "hisp", "afam", "other", "disability_bin")

frmla <- as.formula(paste("surg ~", paste(names, collapse = "+"), sep=""))
fn <- lm(frmla, data=data.1)

a <- waldtest(fs, fn)$F[2]
b <- summary(fs)$coef[2,]
round(b[1], 3)
a.on.z$res[3,]


# Compare estimates of LATE
## Takes 2-3 minutes
plot.df1 <- late.comparison(Y, A, Z, X, late.five)
plot.df1$category <- c("UNADJUSTED TSLS", "ADJUSTED TSLS", "DRML")

p1 <- ggplot(plot.df1, aes(x=category, y=est, ymax=est+1.96*se, ymin=est-1.96*se))+
         geom_point() + geom_errorbar(width=0.3) + 
         geom_hline(yintercept = 0, linetype="dashed") +
         theme_minimal()+
         ylab("Est. Treatment Effect") + xlab("") +
        theme(axis.title = element_text(size = 15)) +
        theme(axis.text = element_text(size = 15)) 
  
p1


# monotonicity sensitivity
defier.max <- (1-a.on.z$res$est[3])/2
defier.seq <- seq(0, defier.max, length.out=50)
eps.seq <- seq(0.0, 0.50, by=0.01)

defier.seq <- rep(defier.seq, length(eps.seq))
eps.seq <- rep(eps.seq, each=length(defier.seq))
est <- late.five$res$est[1] + defier.seq*eps.seq/a.on.z$res$est[3]

plot.data <- data.frame(delta=defier.seq, eps=eps.seq, level=est)

p2 <- ggplot(plot.data, aes(delta, eps)) +
  geom_raster(aes(fill=level)) +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  geom_contour(aes(z = level), 
               breaks = 0, col = 'black', linetype="dashed") +
  theme_minimal() + theme(legend.position = "bottom", legend.key.width= unit(2, 'cm')) +    
  labs(x="P(Defier)",
       y="Average Effect Difference Btw Compliers and Defiers",
       title = "",
       fill = "LATE Upper Bound") +
        theme(axis.title = element_text(size = 15)) +
        theme(axis.text = element_text(size = 15)) 
p2 
                     
 
## Profiling Complier Types

# bootstrap sample 
age.range <- 20:80 # grid for estimating cate
comorb.range <- 0:10  # grid for estimating cate
boot.sample <- 10 # has to be much larger but takes time (I used 5000 for the report) 


# Age
age.profiling <- lapply(age.range, 
                        function(v0){
                        iv.profiling.density(data.frame(X)$age, v0, a.on.z)})
                        
                        
comp.est <- unlist(lapply(age.profiling, function(res){res$est[1]}))
comp.se <- unlist(lapply(age.profiling, function(res){res$se[1]}))
at.est <- unlist(lapply(age.profiling, function(res){res$est[2]}))
at.se <- unlist(lapply(age.profiling, function(res){res$se[2]}))
nt.est <- unlist(lapply(age.profiling, function(res){res$est[3]}))
nt.se <- unlist(lapply(age.profiling, function(res){res$se[3]}))

plot.df <- data.frame(x=age.range, ymin=comp.est-1.96*comp.se, 
                      ymax=comp.est+1.96*comp.se,y=comp.est,
                      type="Complier")
ggplot(plot.df) + geom_errorbar(aes(x=x, ymin=ymin, ymax=ymax)) +
  geom_line(aes(x=x, y=y), size=1) +
  facet_grid(~type) + theme(legend.position = "none") +
  xlab("Age") + ylab("Estimated smoothed density") +
  theme_minimal() + ylim(c(0, .02))
  
                          
plot.df <- rbind(plot.df, data.frame(x=age.range, ymin=at.est-1.96*at.se, 
                                     ymax=at.est+1.96*at.se, y=at.est,
                                     type="Always Taker"))
ggplot(plot.df) + geom_errorbar(aes(x=x, ymin=ymin, ymax=ymax)) +
  geom_line(aes(x=x, y=y), size=1) +
  facet_grid(~type) + theme(legend.position = "none") +
  xlab("Age") + ylab("Estimated smoothed density") +
  theme_minimal() + ylim(c(0, .02))
  
                                       
plot.df <- rbind(plot.df, data.frame(x=age.range, ymin=nt.est-1.96*nt.se, 
                                     ymax=nt.est+1.96*nt.se, 
                                     y=nt.est,
                                     type="Never Taker"))
                                     
p3 <- ggplot(plot.df) + geom_errorbar(aes(x=x, ymin=ymin, ymax=ymax)) +
  geom_line(aes(x=x, y=y), size=1) +
  facet_grid(~type) + theme(legend.position = "none") +
  xlab("Age") + ylab("Estimated Smoothed Density") +
  theme_minimal() + ylim(c(0, .022))                                     

p3

                                     
## Sepsis
sepsis.profile <- lapply(1, function(v0){iv.profiling.prob(data.frame(X)$angus, v0, a.on.z)})
disability.profile <- lapply(1, function(v0){iv.profiling.prob(data.frame(X)$disability_bin, v0, a.on.z)})
afam.profile <- lapply(1, function(v0){iv.profiling.prob(data.frame(X)$afam, v0, a.on.z)})
## Obesity
obese.profile <- lapply(1, function(v0){iv.profiling.prob(data.frame(X)$ynel22, v0, a.on.z)})
## Anemia
anemia.profile <- lapply(1, function(v0){iv.profiling.prob(data.frame(X)$ynel26, v0, a.on.z)})
medicare.profile <- lapply(1, function(v0){iv.profiling.prob(data.frame(X)$p_cat1, v0, a.on.z)})
comorb.profile <- lapply(median(data.1$comorb), function(v0){iv.profiling.density(data.frame(X)$comorb, v0, a.on.z)})

profile.format <- function(x){
	raw <- unlist(x)
	raw <- round(as.numeric(raw[-c(1:3)]), 3)
	comp.est <- raw[1]
	at.est <- raw[2]
	nt.est <- raw[3]
    comp <- data.frame(y = comp.est, ymin = raw[7], ymax = raw[10], type = "Complier")
    at <- data.frame(y = at.est, ymin = raw[8], ymax = raw[11], type = "Always Taker")
    nt <- data.frame(y = nt.est, ymin = raw[9], ymax = raw[12], type = "Never Taker")
    plot.dat <- rbind(comp, at, nt)
 return(plot.dat)
}

plot.sep <- profile.format(sepsis.profile)
p4.1 <- ggplot(plot.sep) + geom_errorbar(aes(x=type, ymin=ymin, ymax=ymax), width=0.3) +
           geom_point(aes(x=type, y=y)) + 
           theme_minimal()+
           theme(legend.position = "none") +
           ylab("Estimated probability of Sepsis") + 
           xlab ('') + coord_cartesian() +
        theme(axis.title = element_text(size = 10)) +
        theme(axis.text = element_text(size = 10)) 
p4.1
  
plot.disab <- profile.format(disability.profile)
p4.2 <- ggplot(plot.disab) + geom_errorbar(aes(x=type, ymin=ymin, ymax=ymax), width=0.3) +
           geom_point(aes(x=type, y=y)) + 
           theme_minimal()+
           theme(legend.position = "none") +
           ylab("Estimated probability of Disability") + 
           xlab ('')  + coord_cartesian() +
        theme(axis.title = element_text(size = 10)) +
        theme(axis.text = element_text(size = 10))
p4.2

plot.afam <- profile.format(afam.profile)
p4.3 <- ggplot(plot.afam) + geom_errorbar(aes(x=type, ymin=ymin, ymax=ymax), width=0.3) +
           geom_point(aes(x=type, y=y)) + 
           theme_minimal()+
           theme(legend.position = "none") +
           ylab("Estimated probability Black") + 
           xlab ('')  + coord_cartesian() +
        theme(axis.title = element_text(size = 10)) +
        theme(axis.text = element_text(size = 10))
p4.3  

plot.obese <- profile.format(obese.profile)
p4.4 <- ggplot(plot.obese) + geom_errorbar(aes(x=type, ymin=ymin, ymax=ymax), width=0.3) +
           geom_point(aes(x=type, y=y)) + 
           theme_minimal()+
           theme(legend.position = "none") +
           ylab("Estimated probability Obese") + 
           xlab ('') + coord_cartesian() +
        theme(axis.title = element_text(size = 10)) +
        theme(axis.text = element_text(size = 10))
p4.4  

plot.anemia <- profile.format(anemia.profile)
p4.5 <- ggplot(plot.anemia) + geom_errorbar(aes(x=type, ymin=ymin, ymax=ymax), width=0.3) +
           geom_point(aes(x=type, y=y)) + 
           theme_minimal()+
           theme(legend.position = "none") +
           ylab("Estimated probability Anemia") + 
           xlab ('') + coord_cartesian() +
        theme(axis.title = element_text(size = 10)) +
        theme(axis.text = element_text(size = 10))
p4.5  

plot.medicare <- profile.format(medicare.profile)
p4.6 <- ggplot(plot.medicare) + geom_errorbar(aes(x=type, ymin=ymin, ymax=ymax), width=0.3) +
           geom_point(aes(x=type, y=y)) + 
           theme_minimal()+
           theme(legend.position = "none") +
           ylab("Estimated probability Medicare Patient") + 
           xlab ('')  + coord_cartesian() +
        theme(axis.title = element_text(size = 10)) +
        theme(axis.text = element_text(size = 10))
p4.6  

pdf("profile.pdf", width=10, height=7, onefile=FALSE, paper="special")
grid.arrange(p4.1, p4.2, p4.3, p4.4, p4.5, p4.6, nrow=2, ncol=3)
dev.off() 


## CATEs 
p7 <- ggplot(data.frame(ite=ite.res$est), aes(x=ite)) + 
  geom_histogram(bins = 50, color="black", fill="lightblue") + 
  xlim(c(-0.7,0.7)) +
  geom_vline(aes(xintercept = late.five$res[1,2], color="LATE"), linetype="dashed", size=1) +
  xlab("Estimated. Individual Treatment Effect") + 
  theme_minimal()+
  theme(legend.position = "none", 
  legend.title= element_blank()) +
  ylab("Frequency") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15))

p7
  

  
# Comorb x sepsis  
load("boot.out")

p8 <- ggplot(pred.df.1, aes(x=comorb, y=cate)) + 
  geom_point(size=1) + geom_line() +
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci))+
  geom_hline(aes(color="red",yintercept=late.five$res$est[1]), linetype="dashed")+
  theme_minimal() + xlab("No. of Comorbidities") + ylab("Conditional Risk of Adverse Event") +
  scale_x_continuous(breaks = seq(0, 10, by=2))+
  scale_colour_manual(name = "", values=c("red"), labels=c("Estimated LATE")) +
  theme(legend.position = "bottom", panel.spacing.x = unit(2, "lines"),
        strip.text.x = element_text(size = 20)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15))
p8


# heat map (no bootstrap)
cov.df <- data.frame(X)
cov.df$if.a <- a.on.z$ifvals[,2]-a.on.z$ifvals[,1]
cov.df$if.y <- y.on.z$ifvals[,2]-y.on.z$ifvals[,1]


gam.ifa.comorb.age.sepsis <- gam(if.a~s(comorb)+s(age)+angus+s(comorb*age)+
                                   s(comorb*angus)+s(age*angus)+s(age*angus*comorb), 
                                 cov.df, family=gaussian)
                                                                 
gam.ify.comorb.age.sepsis <- gam(if.y~s(comorb)+s(age)+angus+s(comorb*age)+
                                   s(comorb*angus) +s(age*angus)+s(age*angus*comorb), 
                                 cov.df, family=gaussian)

pred.df <- data.frame(age=rep(age.range, each = length(comorb.range)), 
                      comorb=rep(comorb.range, length(age.range)), angus=1)
pred.df <- rbind(pred.df, data.frame(age=rep(age.range, each = length(comorb.range)), 
                            comorb=rep(comorb.range, length(age.range)), 
                            angus=0))

cate.comorb.age.sepsis <- predict(gam.ify.comorb.age.sepsis, pred.df)/ predict(gam.ifa.comorb.age.sepsis, pred.df)

pred.df$cate <- cate.comorb.age.sepsis

p9 <- ggplot(pred.df) + 
  geom_tile(aes(x=age,y=factor(comorb), fill=cate)) + 
  theme_minimal()+
  theme(legend.position = "bottom") + 
  facet_grid(~factor(angus), labeller = as_labeller(c(`0`="Not Septic", `1`="Septic")))+ 
  scale_fill_gradient2("Est. CATE",low="navy", high="coral")+guides(fill = guide_colourbar(barwidth =10))+
  ylab("# of comorbidities") + xlab("Age") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15))


p9

vars <- c("age", "comorb", "age_gr6", "p_cat1", "p_cat3", "ynel1",
          "ynel4", "ynel5", "ynel13", "ynel14", "ynel15", "ynel27", "ynel31", 
          "shafi1", "shafi3", "shafi20", "shafi23", "shafi40", "shafi41", "shafi49")

xmat <- match.data[vars]
xmat <- as.matrix(xmat)
y=match.data$died; # Mortality
d=match.data$surg;
z=match.data$pref_bin;
modely.z=glm(y~z+xmat,family=binomial)
modeld.z=glm(d~z+xmat,family=binomial)
modelz.x=glm(z~xmat,family=binomial)

z1mat=cbind(rep(1,nrow(xmat)), rep(1,nrow(xmat)), xmat);
z0mat=cbind(rep(1,nrow(xmat)), rep(0,nrow(xmat)), xmat);

evy.z1=exp(z1mat%*%matrix(coef(modely.z),ncol=1))/(1+exp(z1mat%*%matrix(coef(modely.z),ncol=1)))
evy.z0=exp(z0mat%*%matrix(coef(modely.z),ncol=1))/(1+exp(z0mat%*%matrix(coef(modely.z),ncol=1)))
evd.z1=exp(z1mat%*%matrix(coef(modeld.z),ncol=1))/(1+exp(z1mat%*%matrix(coef(modeld.z),ncol=1)))
evd.z0=exp(z0mat%*%matrix(coef(modeld.z),ncol=1))/(1+exp(z0mat%*%matrix(coef(modeld.z),ncol=1)))


# Estimate proportion of "compliers"
pr.comp <- mean(evd.z1)-mean(evd.z0)

modeld.z1=glm(d~xmat,subset=(z==1),family=binomial)
modeld.z0=glm(d~xmat,subset=(z==0),family=binomial)
z1mat=cbind(rep(1,nrow(xmat)), xmat);
z0mat=cbind(rep(1,nrow(xmat)), xmat);
evd.z1=exp(z1mat%*%matrix(coef(modeld.z1),ncol=1))/(1+exp(z1mat%*%matrix(coef(modeld.z1),ncol=1)))
evd.z0=exp(z0mat%*%matrix(coef(modeld.z0),ncol=1))/(1+exp(z0mat%*%matrix(coef(modeld.z0),ncol=1)))

# Estimated Proportion of Defierss
mean(evd.z0 > evd.z1)

