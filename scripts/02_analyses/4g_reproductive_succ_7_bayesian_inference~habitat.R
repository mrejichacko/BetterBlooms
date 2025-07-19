#####################################################################
###
### 4g. Pollinator experiment: 
### Modelling plant reproductive success: bayesian inference
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
### 
#####################################################################
#####################################################################
#####################################################################
###
### Prep

###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
### load 3 a and 4a-f! it's not enough to just load 3a!!!!!!!!!!!!###
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###

#Draw Baysian inference
library(arm)
nsim <- 10000

#Carrot: Seed set
bsim.mod.1_Carrot <- sim(mod.1_Carrot, n.sim=nsim)
coeff.mod.1_Carrot <- as.data.frame(t(as.data.frame(apply(bsim.mod.1_Carrot@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.1_Carrot$fixeff <- rownames(coeff.mod.1_Carrot)
row.names(coeff.mod.1_Carrot)<-NULL
coeff.mod.1_Carrot$Model <- "mod.1_Carrot"
coeff.mod.1_Carrot$Phytometer <- "Carrot"
coeff.mod.1_Carrot$Response <- "Seed_set"
coeff.mod.1_Carrot$Effect <- "GardenLandscape"
colnames(coeff.mod.1_Carrot)[1] <- "lower"
colnames(coeff.mod.1_Carrot)[2] <- "mean"
colnames(coeff.mod.1_Carrot)[3] <- "upper"
coeff.mod.1_Carrot

#Radish: fruit set
bsim.mod.1_Radish <- sim(mod.1_Radish, n.sim=nsim) #simulate from the posterior distribution
coeff.mod.1_Radish <- as.data.frame(t(as.data.frame(apply(bsim.mod.1_Radish@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.1_Radish$fixeff <- rownames(coeff.mod.1_Radish)
row.names(coeff.mod.1_Radish)<-NULL
coeff.mod.1_Radish$Model <- "mod.1_Radish"
coeff.mod.1_Radish$Phytometer <- "Radish"
coeff.mod.1_Radish$Response <- "Fruit_set"
coeff.mod.1_Radish$Effect <- "GardenLandscape"
colnames(coeff.mod.1_Radish)[1] <- "lower"
colnames(coeff.mod.1_Radish)[2] <- "mean"
colnames(coeff.mod.1_Radish)[3] <- "upper"
coeff.mod.1_Radish

#Radish: seed set
bsim.mod.1_Radish_seedset <- sim(mod.1_Radish_seedset, n.sim=nsim)
coeff.mod.1_Radish_seedset <- as.data.frame(t(as.data.frame(apply(bsim.mod.1_Radish_seedset@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.1_Radish_seedset$fixeff <- rownames(coeff.mod.1_Radish_seedset)
row.names(coeff.mod.1_Radish_seedset)<-NULL
coeff.mod.1_Radish_seedset$Model <- "mod.1_Radish_seedset"
coeff.mod.1_Radish_seedset$Phytometer <- "Radish"
coeff.mod.1_Radish_seedset$Response <- "Seed_set"
coeff.mod.1_Radish_seedset$Effect <- "GardenLandscape"
colnames(coeff.mod.1_Radish_seedset)[1] <- "lower"
colnames(coeff.mod.1_Radish_seedset)[2] <- "mean"
colnames(coeff.mod.1_Radish_seedset)[3] <- "upper"
coeff.mod.1_Radish_seedset

#Sainfoin: Fruit set: full model
bsim.mod.1_Sainfoin <- sim(mod.1_Sainfoin, n.sim=nsim)
coeff.mod.1_Sainfoin <- as.data.frame(t(as.data.frame(apply(bsim.mod.1_Sainfoin@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.1_Sainfoin$fixeff <- rownames(coeff.mod.1_Sainfoin)
row.names(coeff.mod.1_Sainfoin)<-NULL
coeff.mod.1_Sainfoin$Model <- "mod.1_Sainfoin"
coeff.mod.1_Sainfoin$Phytometer <- "Sainfoin"
coeff.mod.1_Sainfoin$Response <- "Fruit_set"
coeff.mod.1_Sainfoin$Effect <- "GardenLandscape"
colnames(coeff.mod.1_Sainfoin)[1] <- "lower"
colnames(coeff.mod.1_Sainfoin)[2] <- "mean"
colnames(coeff.mod.1_Sainfoin)[3] <- "upper"
coeff.mod.1_Sainfoin

#Comfrey: fruit set
bsim.mod.1_Comfrey <- sim(mod.1_Comfrey, n.sim=nsim)
coeff.mod.1_Comfrey <- as.data.frame(t(as.data.frame(apply(bsim.mod.1_Comfrey@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.1_Comfrey
coeff.mod.1_Comfrey$fixeff <- rownames(coeff.mod.1_Comfrey)
row.names(coeff.mod.1_Comfrey)<-NULL
coeff.mod.1_Comfrey
coeff.mod.1_Comfrey$Model <- "mod.1_Comfrey"
coeff.mod.1_Comfrey$Phytometer <- "Comfrey"
coeff.mod.1_Comfrey$Response <- "Fruit_set"
coeff.mod.1_Comfrey$Effect <- "GardenLandscape"
colnames(coeff.mod.1_Comfrey)[1] <- "lower"
colnames(coeff.mod.1_Comfrey)[2] <- "mean"
colnames(coeff.mod.1_Comfrey)[3] <- "upper"
coeff.mod.1_Comfrey

#Comfrey: seed set
bsim.mod.1_Comfrey_seedset <- sim(mod.1_Comfrey_seedset, n.sim=nsim)
coeff.mod.1_Comfrey_seedset <- as.data.frame(t(as.data.frame(apply(bsim.mod.1_Comfrey_seedset@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.1_Comfrey_seedset
coeff.mod.1_Comfrey_seedset$fixeff <- rownames(coeff.mod.1_Comfrey_seedset)
row.names(coeff.mod.1_Comfrey_seedset)<-NULL
coeff.mod.1_Comfrey_seedset
coeff.mod.1_Comfrey_seedset$Model <- "mod.1_Comfrey_seedset"
coeff.mod.1_Comfrey_seedset$Phytometer <- "Comfrey"
coeff.mod.1_Comfrey_seedset$Response <- "Seed_set"
coeff.mod.1_Comfrey_seedset$Effect <- "GardenLandscape"
colnames(coeff.mod.1_Comfrey_seedset)[1] <- "lower"
colnames(coeff.mod.1_Comfrey_seedset)[2] <- "mean"
colnames(coeff.mod.1_Comfrey_seedset)[3] <- "upper"
coeff.mod.1_Comfrey_seedset


#Combine them
coeff.mods.ReprSucc.GxL <- rbind(coeff.mod.1_Carrot, coeff.mod.1_Radish, coeff.mod.1_Radish_seedset, coeff.mod.1_Sainfoin, coeff.mod.1_Comfrey, coeff.mod.1_Comfrey_seedset)
head(coeff.mods.ReprSucc.GxL)
tail(coeff.mods.ReprSucc.GxL)

coeff.mods.ReprSucc.GxL$fixeff <- as.factor(coeff.mods.ReprSucc.GxL$fixeff)

#Recode the variables names: rename bumblebees in social bees
library(dplyr)
coeff.mods.ReprSucc.GxL$fixeff <- recode_factor(coeff.mods.ReprSucc.GxL$fixeff, 
                                    'PlantS.z' = "Garden", 
                                    'PC1' = "Landscape",
                                    'PlantS.z:PC1'="GardenxLandscape",
                                    '(Intercept)' = "Intercept")

coeff.mods.ReprSucc.GxL.table <- cbind(coeff.mods.ReprSucc.GxL[,c(4:8)], as.data.frame(apply(coeff.mods.ReprSucc.GxL[,1:3], 2, function(x)round(x,2))))

coeff.mods.ReprSucc.GxL.table <- coeff.mods.ReprSucc.GxL.table[,c(7,6,8,2,3,5,4,1)] 

coeff.mods.ReprSucc.GxL.table

write.table(coeff.mods.ReprSucc.GxL, 
            file = "results/beta_coefficients_reproductive_success_garden_landscape.txt", 
            append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
