
library(foreign)
library(dplyr)
library(xtable)

rm(list=ls())

data <- read.dta("iv_match.dta")

        
## Balance
all.cov <- c("comorb", "angus", "age","ynel1","ynel2","ynel3","ynel4", "ynel5",
            "ynel6","ynel7","ynel8","ynel9","ynel10","ynel11","ynel12","ynel13",
            "ynel14","ynel15","ynel16","ynel17","ynel18","ynel19","ynel20","ynel21",
            "ynel22","ynel23","ynel24","ynel25","ynel26","ynel27","ynel28",
            "ynel29","ynel30", "ynel31", "p_cat1", "p_cat2", "p_cat3","p_cat4",
            "p_cat5", "female", "hisp", "afam", "other", "disability_bin")
             
data.var <- data %>% group_by(surg) %>% 
   summarize(across(all.cov, ~var(.x))) %>% as.data.frame()

c.var <- as.numeric(data.var[1,])
t.var <- as.numeric(data.var[2,])
c.var <- c.var[-1]
t.var <- t.var[-1]
pooled.sd <- sqrt((t.var + c.var/2))


um.wt <- data %>% group_by(surg) %>% 
   summarize(across(all.cov, ~mean(.x))) %>% as.data.frame()

           
um.wt.tab <- matrix(NA, length(all.cov), 3)
um.wt.tab[,1] <- unlist(um.wt[1,-1]) 
um.wt.tab[,2] <- unlist(um.wt[2,-1])                        
um.wt.tab[,3] <- (unlist(um.wt[2,-1]) - unlist(um.wt[1,-1]))/pooled.sd


### Table
n_covs <- length(all.cov)

var_names <- c("No. Comorbidities", "Sepsis", "Age",
               "Congestive Heart Failure",  "Cardiac Arrhythmias","Valvular Disease",
               "Pulmonary Circulation Disorders", "Peripheral Vascular Disorders",
               "Hypertension, Uncomplicated","Paralysis", "Other Neurological Disorders", 
               "Chronic Pulmonary Disease","Diabetes, Uncomplicated",
               "Diabetes, Complicated","Hypothyroidism", "Renal Failure","Liver Disease",
               "Peptic Ulcer Disease Excluding Bleeding","AIDS/HIV", "Lymphoma","Metastatic Cancer",
               "Solid Tumor Without Metastasis","Rheumatoid Arthritis/Collagen Vascular",
               "Coagulopathy","Obesity","Weight Loss","Fluid and Electrolyte Disorders",
               "Blood Loss Anemia","Deficiency Anemia","Alcohol Abuse","Drug Abuse",
               "Psychoses","Depression","Hypertension, Complicated",
               "Medicare", "Medicaid", "Private Insurance", "Self Insurance", "Other",
               "Female", "Hispanic", "Black", "Other Race","Disability")
              
rownames(um.wt.tab) <- var_names
colnames(um.wt.tab) <- c("Operative", "Non-Operative", "Std. Diff.")

tab.1 <- um.wt.tab[which(abs(um.wt.tab[,3]) > .15),]
tab.1[,1] <- c(round(tab.1[1,1], 1), round((tab.1[2:9,1])*100, 1))
tab.1[,2] <- c(round(tab.1[1,2], 1), round((tab.1[2:9,2])*100, 1))

tab.1
xtable(tab.1, digits = c(1,1,1,2))            
            
tab.2 <- um.wt.tab
nrow(um.wt.tab)
tab.2[,1] <- c(round(tab.2[1,1], 1), round((tab.2[2:94,1])*100, 1))
tab.2[,2] <- c(round(tab.2[1,2], 1), round((tab.2[2:94,2])*100, 1))

tab.2
xtable(tab.2, digits = c(1,1,1,2))               
            