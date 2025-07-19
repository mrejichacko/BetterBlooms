#####################################################################
#####################################################################
#####################################################################
###
### 3h. Modelling plant reproductive success:  Bayesian inference ###
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
#####################################################################
#####################################################################
#####################################################################
###
### Prep

#####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
### RUN 3a-g!!! It's not enough to load the environment, you need all the models!! ####
#####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###

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
coeff.mod.1_Carrot$Effect <- "Abundance"
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
coeff.mod.1_Radish$Effect <- "Abundance"
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
coeff.mod.1_Radish_seedset$Effect <- "Abundance"
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
coeff.mod.1_Sainfoin$Effect <- "Abundance"
colnames(coeff.mod.1_Sainfoin)[1] <- "lower"
colnames(coeff.mod.1_Sainfoin)[2] <- "mean"
colnames(coeff.mod.1_Sainfoin)[3] <- "upper"
coeff.mod.1_Sainfoin

#Sainfoin: Fruit set: combined model (all bees together)
bsim.mod.1b_Sainfoin <- sim(mod.1b_Sainfoin, n.sim=nsim)
coeff.mod.1b_Sainfoin <- as.data.frame(t(as.data.frame(apply(bsim.mod.1b_Sainfoin@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.1b_Sainfoin$fixeff <- rownames(coeff.mod.1b_Sainfoin)
row.names(coeff.mod.1b_Sainfoin)<-NULL
coeff.mod.1b_Sainfoin$Model <- "mod.1b_Sainfoin"
coeff.mod.1b_Sainfoin$Phytometer <- "Sainfoin"
coeff.mod.1b_Sainfoin$Response <- "Fruit_set"
coeff.mod.1b_Sainfoin$Effect <- "Abundance"
colnames(coeff.mod.1b_Sainfoin)[1] <- "lower"
colnames(coeff.mod.1b_Sainfoin)[2] <- "mean"
colnames(coeff.mod.1b_Sainfoin)[3] <- "upper"
coeff.mod.1b_Sainfoin

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
coeff.mod.1_Comfrey$Effect <- "Abundance"
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
coeff.mod.1_Comfrey_seedset$Effect <- "Abundance"
colnames(coeff.mod.1_Comfrey_seedset)[1] <- "lower"
colnames(coeff.mod.1_Comfrey_seedset)[2] <- "mean"
colnames(coeff.mod.1_Comfrey_seedset)[3] <- "upper"
coeff.mod.1_Comfrey_seedset

##################################################

#Carrot: Seed set
bsim.mod.2_Carrot <- sim(mod.2_Carrot, n.sim=nsim)
coeff.mod.2_Carrot <- as.data.frame(t(as.data.frame(apply(bsim.mod.2_Carrot@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.2_Carrot$fixeff <- rownames(coeff.mod.2_Carrot)
row.names(coeff.mod.2_Carrot)<-NULL
coeff.mod.2_Carrot$Model <- "mod.2_Carrot"
coeff.mod.2_Carrot$Phytometer <- "Carrot"
coeff.mod.2_Carrot$Response <- "Seed_set"
coeff.mod.2_Carrot$Effect <- "SpeciesRichness"
colnames(coeff.mod.2_Carrot)[1] <- "lower"
colnames(coeff.mod.2_Carrot)[2] <- "mean"
colnames(coeff.mod.2_Carrot)[3] <- "upper"
coeff.mod.2_Carrot

#Radish: fruit set
bsim.mod.2_Radish <- sim(mod.2_Radish, n.sim=nsim) #simulate from the posterior distribution
coeff.mod.2_Radish <- as.data.frame(t(as.data.frame(apply(bsim.mod.2_Radish@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.2_Radish$fixeff <- rownames(coeff.mod.2_Radish)
row.names(coeff.mod.2_Radish)<-NULL
coeff.mod.2_Radish$Model <- "mod.2_Radish"
coeff.mod.2_Radish$Phytometer <- "Radish"
coeff.mod.2_Radish$Response <- "Fruit_set"
coeff.mod.2_Radish$Effect <- "SpeciesRichness"
colnames(coeff.mod.2_Radish)[1] <- "lower"
colnames(coeff.mod.2_Radish)[2] <- "mean"
colnames(coeff.mod.2_Radish)[3] <- "upper"
coeff.mod.2_Radish

#Radish: seed set
bsim.mod.2_Radish_seedset <- sim(mod.2_Radish_seedset, n.sim=nsim)
coeff.mod.2_Radish_seedset <- as.data.frame(t(as.data.frame(apply(bsim.mod.2_Radish_seedset@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.2_Radish_seedset$fixeff <- rownames(coeff.mod.2_Radish_seedset)
row.names(coeff.mod.2_Radish_seedset)<-NULL
coeff.mod.2_Radish_seedset$Model <- "mod.2_Radish_seedset"
coeff.mod.2_Radish_seedset$Phytometer <- "Radish"
coeff.mod.2_Radish_seedset$Response <- "Seed_set"
coeff.mod.2_Radish_seedset$Effect <- "SpeciesRichness"
colnames(coeff.mod.2_Radish_seedset)[1] <- "lower"
colnames(coeff.mod.2_Radish_seedset)[2] <- "mean"
colnames(coeff.mod.2_Radish_seedset)[3] <- "upper"
coeff.mod.2_Radish_seedset

#Sainfoin: Fruit set: full model
bsim.mod.2_Sainfoin <- sim(mod.2_Sainfoin, n.sim=nsim)
coeff.mod.2_Sainfoin <- as.data.frame(t(as.data.frame(apply(bsim.mod.2_Sainfoin@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.2_Sainfoin$fixeff <- rownames(coeff.mod.2_Sainfoin)
row.names(coeff.mod.2_Sainfoin)<-NULL
coeff.mod.2_Sainfoin$Model <- "mod.2_Sainfoin"
coeff.mod.2_Sainfoin$Phytometer <- "Sainfoin"
coeff.mod.2_Sainfoin$Response <- "Fruit_set"
coeff.mod.2_Sainfoin$Effect <- "SpeciesRichness"
colnames(coeff.mod.2_Sainfoin)[1] <- "lower"
colnames(coeff.mod.2_Sainfoin)[2] <- "mean"
colnames(coeff.mod.2_Sainfoin)[3] <- "upper"
coeff.mod.2_Sainfoin

#Sainfoin: Fruit set: combined model (all bees together)
bsim.mod.2b_Sainfoin <- sim(mod.2b_Sainfoin, n.sim=nsim)
coeff.mod.2b_Sainfoin <- as.data.frame(t(as.data.frame(apply(bsim.mod.2b_Sainfoin@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.2b_Sainfoin$fixeff <- rownames(coeff.mod.2b_Sainfoin)
row.names(coeff.mod.2b_Sainfoin)<-NULL
coeff.mod.2b_Sainfoin$Model <- "mod.2b_Sainfoin"
coeff.mod.2b_Sainfoin$Phytometer <- "Sainfoin"
coeff.mod.2b_Sainfoin$Response <- "Fruit_set"
coeff.mod.2b_Sainfoin$Effect <- "SpeciesRichness"
colnames(coeff.mod.2b_Sainfoin)[1] <- "lower"
colnames(coeff.mod.2b_Sainfoin)[2] <- "mean"
colnames(coeff.mod.2b_Sainfoin)[3] <- "upper"
coeff.mod.2b_Sainfoin

#Comfrey: fruit set
bsim.mod.2_Comfrey <- sim(mod.2_Comfrey, n.sim=nsim)
coeff.mod.2_Comfrey <- as.data.frame(t(as.data.frame(apply(bsim.mod.2_Comfrey@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.2_Comfrey
coeff.mod.2_Comfrey$fixeff <- rownames(coeff.mod.2_Comfrey)
row.names(coeff.mod.2_Comfrey)<-NULL
coeff.mod.2_Comfrey
coeff.mod.2_Comfrey$Model <- "mod.2_Comfrey"
coeff.mod.2_Comfrey$Phytometer <- "Comfrey"
coeff.mod.2_Comfrey$Response <- "Fruit_set"
coeff.mod.2_Comfrey$Effect <- "SpeciesRichness"
colnames(coeff.mod.2_Comfrey)[1] <- "lower"
colnames(coeff.mod.2_Comfrey)[2] <- "mean"
colnames(coeff.mod.2_Comfrey)[3] <- "upper"
coeff.mod.2_Comfrey

#Comfrey: seed set
bsim.mod.2_Comfrey_seedset <- sim(mod.2_Comfrey_seedset, n.sim=nsim)
coeff.mod.2_Comfrey_seedset <- as.data.frame(t(as.data.frame(apply(bsim.mod.2_Comfrey_seedset@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.2_Comfrey_seedset
coeff.mod.2_Comfrey_seedset$fixeff <- rownames(coeff.mod.2_Comfrey_seedset)
row.names(coeff.mod.2_Comfrey_seedset)<-NULL
coeff.mod.2_Comfrey_seedset
coeff.mod.2_Comfrey_seedset$Model <- "mod.2_Comfrey_seedset"
coeff.mod.2_Comfrey_seedset$Phytometer <- "Comfrey"
coeff.mod.2_Comfrey_seedset$Response <- "Seed_set"
coeff.mod.2_Comfrey_seedset$Effect <- "SpeciesRichness"
colnames(coeff.mod.2_Comfrey_seedset)[1] <- "lower"
colnames(coeff.mod.2_Comfrey_seedset)[2] <- "mean"
colnames(coeff.mod.2_Comfrey_seedset)[3] <- "upper"
coeff.mod.2_Comfrey_seedset

##############################################

#Combine them
coeff.mods.ReprSucc <- rbind(coeff.mod.1_Carrot, coeff.mod.1_Radish, coeff.mod.1_Radish_seedset, coeff.mod.1_Sainfoin, coeff.mod.1b_Sainfoin, coeff.mod.1_Comfrey, coeff.mod.1_Comfrey_seedset,
                               coeff.mod.2_Carrot, coeff.mod.2_Radish, coeff.mod.2_Radish_seedset, coeff.mod.2_Sainfoin, coeff.mod.2b_Sainfoin)
head(coeff.mods.ReprSucc)
tail(coeff.mods.ReprSucc)

coeff.mods.ReprSucc$fixeff <- as.factor(coeff.mods.ReprSucc$fixeff)

#Recode the variables names: rename bumblebees in social bees
library(dplyr)
coeff.mods.ReprSucc$fixeff <- recode_factor(coeff.mods.ReprSucc$fixeff, 
                                    'A_Apis_Carrot.dayly.z' = "Honeybees", 
                                    'A_socialBees_Carrot.dayly.z' = "social_Bees",
                                    'A_solitaryBees_Carrot.dayly.z' = "solitary_Bees",
                                    'A_otherAculeata_Carrot.dayly.z' = "Wasps",
                                    'A_Syrphidae_Carrot.dayly.z' = "Hoverflies",
                                    'A_Coleoptera_Carrot.dayly.z' = "Beetles",
                                    'A_Apis_Radish.dayly.z' = "Honeybees", 
                                    'A_socialBees_Radish.dayly.z' = "social_Bees",
                                    'A_solitaryBees_Radish.dayly.z' = "solitary_Bees",
                                    'A_otherAculeata_Radish.dayly.z' = "Wasps",
                                    'A_Syrphidae_Radish.dayly.z' = "Hoverflies",
                                    'A_Coleoptera_Radish.dayly.z' = "Beetles",
                                    'A_Coleoptera_Radish.dayly.t.z' = "Beetles",
                                    'A_Apis_Sainfoin.dayly.z' = "Honeybees",
                                    'A_solitaryBees_Sainfoin.dayly.z' ="solitary_Bees",
                                    'A_allBees_Sainfoin.dayly.z' = "all_Bees",
                                    'A_Bombus_Sainfoin.dayly.z' = "social_Bees",
                                    'A_Bombus_Comfrey.dayly.1.z' = "social_Bees", 
                                    'S_socialBees_Carrot.dayly.z' = "social_Bees",
                                    'S_solitaryBees_Carrot.dayly.z' = "solitary_Bees",
                                    'S_otherAculeata_Carrot.dayly.z' = "Wasps",
                                    'S_Syrphidae_Carrot.dayly.z' = "Hoverflies",
                                    'S_Coleoptera_Carrot.dayly.z' = "Beetles", 
                                    'S_socialBees_Radish.dayly.z' = "social_Bees",
                                    'S_solitaryBees_Radish.dayly.z' = "solitary_Bees",
                                    'S_Syrphidae_Radish.dayly.z' = "Hoverflies",
                                    'S_solitaryBees_Sainfoin.dayly.z' ="solitary_Bees",
                                    'S_allBees_Sainfoin.dayly.z' = "all_Bees",
                                    'S_Bombus_Sainfoin.dayly.z' = "social_Bees",
                                    '(Intercept)' = "Intercept")

coeff.mods.ReprSucc

coeff.mods.ReprSucc.table <- cbind(coeff.mods.ReprSucc[,c(4:8)], as.data.frame(apply(coeff.mods.ReprSucc[,1:3], 2, function(x)round(x,2))))

coeff.mods.ReprSucc.table <- coeff.mods.ReprSucc.table[,c(2,3,5,4,1,6,7,8)] 

coeff.mods.ReprSucc.table

write.table(coeff.mods.ReprSucc.table, file = "results/beta_coefficients_reproductive_success_all.txt", append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

