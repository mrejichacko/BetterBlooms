#####################################################################
#####################################################################
#####################################################################
###
### 3e. Pollinator experiment: 
### Modelling plant reproductive success: Sainfoin: Fruit set
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
#####################################################################
#####################################################################
#####################################################################
###
### Prep


#################### RUN 3a or load from environment!!! #############

rm(list=ls())
gc()

load("environments/3a_environment.RData")

########################
# Combine the datasets #
########################

df7d <- merge(datExpl, df6d)

str(df7d)
names(df7d)

#Kick out garden Nr. 39
df7d <- droplevels(df7d[which(df7d$Id != 39),])
dim(df7d)

#Transform factor variables
df7d$Id.fac <- as.factor(df7d$Id)
df7d$Plant_Id.fac <- as.factor(df7d$plant_id)

###############################
## Prepare the Sainfoin data ##
###############################

## We have Plants with 0 inflorescences, exclude them
## Exclude gardens with too few plants (19 & 28; only one plant scorable)

df8 <- subset(df7d, n_inflorescences_assessed != 0 & Id.fac != 19 & Id.fac != 28 & Id.fac != 52) #  
dim(df8)

###################
## Run the model ##
###################

#################
### Abundance ###
#################

#Abundance model: a binomial glmer (failure vrs. success); with dayly rates!

#mod.1_Sainfoin <- glmer(cbind(n_flowers_with_fruits, n_flowers_without_fruits) ~ A_Apis_Sainfoin.dayly.z +  A_Bombus_Sainfoin.dayly.z + A_solitaryBees_Sainfoin.dayly.z + (1|Id.fac/Plant_Id.fac), data=df8, family="binomial")
mod.1_Sainfoin <- do.call(glmer, list(
  formula = as.formula("cbind(n_flowers_with_fruits, n_flowers_without_fruits) ~ A_Apis_Sainfoin.dayly.z + A_Bombus_Sainfoin.dayly.z + A_solitaryBees_Sainfoin.dayly.z + (1 | Id.fac / Plant_Id.fac)"),
  data = df8,
  family = binomial()
))

summary(mod.1_Sainfoin) #For the coefficient plots

#this might look silly but mumin had issues with spaces/line breaks
mod.1b_Sainfoin <- do.call(glmer, list(
  formula = as.formula("cbind(n_flowers_with_fruits, n_flowers_without_fruits) ~ A_allBees_Sainfoin.dayly.z + (1 | Id.fac / Plant_Id.fac)"),
  data = df8,
  family = binomial()
))

summary(mod.1b_Sainfoin) #All together significant: for the effect plotsS

r.squaredGLMM(mod.1b_Sainfoin)

#### Model 1 ######################################################

library(car)
vif(mod.1_Sainfoin) #No extreme values: >5 would be problematic

#Assessing model assumptions: (check page 144)
par(mfrow=c(2,3))

#fitted vs. residuals
scatter.smooth(fitted(mod.1_Sainfoin), resid(mod.1_Sainfoin)) #this looks quite good (no strong positive correlation = no strong shrinkage)
abline(h=0, lty=2, col="red") #Some values do not fit the data well (very large residuals) -> for large seed sets the fitted values are not very reliable. 
title("Tukey-Anscombe Plot")

#qq-plot of residuals
qqnorm(resid(mod.1_Sainfoin)) #qq-plot of residuals
qqline(resid(mod.1_Sainfoin))

#qq-plot  of random effects: Id
qqnorm(ranef(mod.1_Sainfoin)$Id.fac[,1])
qqline(ranef(mod.1_Sainfoin)$Id.fac[,1])

#qq-plot of the random effects: Plant
qqnorm(ranef(mod.1_Sainfoin)$Plant_Id.fac[,1]) 
qqline(ranef(mod.1_Sainfoin)$Plant_Id.fac[,1])

#fitted vs. observed
plot(fitted(mod.1_Sainfoin), jitter(df8$n_flowers_with_fruits/(df8$n_flowers_with_fruits + df8$n_flowers_without_fruits),0.05))
abline(0,1) 

dev.off()  

#Heteroscedasticity (= nonhomogeneity of the residual variance)
par(mfrow=c(1,3))
scatter.smooth(df8$A_Apis_Sainfoin.dayly.z, resid(mod.1_Sainfoin)); abline(0,0, lty=2, col="red")
scatter.smooth(df8$A_Bombus_Sainfoin.dayly.z, resid(mod.1_Sainfoin)); abline(0,0, lty=2, col="red")
scatter.smooth(df8$A_solitaryBees_Sainfoin.dayly.z, resid(mod.1_Sainfoin)); abline(0,0, lty=2, col="red")
dev.off() #There are some outlyers

#Check if random effects are 0:
mean(ranef(mod.1_Sainfoin)$Id.fac[,1])
mean(ranef(mod.1_Sainfoin)$Plant_Id.fac[,1])

#Check for overdispersion:
dispersion_glmer(mod.1_Sainfoin) #without observation-level RF we have overdispersion.

#### Model 1b ######################################################

#Assessing model assumptions: (check page 144)
par(mfrow=c(2,3))

#fitted vs. residuals
scatter.smooth(fitted(mod.1b_Sainfoin), resid(mod.1b_Sainfoin)) #this looks quite good (no strong positive correlation = no strong shrinkage)
abline(h=0, lty=2, col="red") #Some values do not fit the data well (very large residuals) -> for large seed sets the fitted values are not very reliable. 
title("Tukey-Anscombe Plot")

#qq-plot of residuals
qqnorm(resid(mod.1b_Sainfoin)) #qq-plot of residuals
qqline(resid(mod.1b_Sainfoin))

#qq-plot  of random effects: Id
qqnorm(ranef(mod.1b_Sainfoin)$Id.fac[,1])
qqline(ranef(mod.1b_Sainfoin)$Id.fac[,1])

#qq-plot of the random effects: Plant
qqnorm(ranef(mod.1b_Sainfoin)$Plant_Id.fac[,1]) 
qqline(ranef(mod.1b_Sainfoin)$Plant_Id.fac[,1])

#fitted vs. observed
plot(fitted(mod.1b_Sainfoin), jitter(df8$n_flowers_with_fruits/(df8$n_flowers_with_fruits + df8$n_flowers_without_fruits),0.05))
abline(0,1) 

dev.off()  

#Heteroscedasticity (= nonhomogeneity of the residual variance)
par(mfrow=c(1,3))
scatter.smooth(df8$A_Apis_Sainfoin.dayly.z, resid(mod.1b_Sainfoin)); abline(0,0, lty=2, col="red")
scatter.smooth(df8$A_Bombus_Sainfoin.dayly.z, resid(mod.1b_Sainfoin)); abline(0,0, lty=2, col="red")
scatter.smooth(df8$A_solitaryBees_Sainfoin.dayly.z, resid(mod.1b_Sainfoin)); abline(0,0, lty=2, col="red")
dev.off() #There are some outlyers

#Check if random effects are 0:
mean(ranef(mod.1b_Sainfoin)$Id.fac[,1])
mean(ranef(mod.1b_Sainfoin)$Plant_Id.fac[,1])

#Check for overdispersion:
dispersion_glmer(mod.1b_Sainfoin) #without observation-level RF we have overdispersion.

#################################################################################################################

########################
### Species richness ###
########################

#Abundance model: a binomial glmer (failure vrs. success); with dayly rates!

mod.2_Sainfoin <- do.call(glmer, list(
  formula = as.formula("cbind(n_flowers_with_fruits, n_flowers_without_fruits) ~ S_Bombus_Sainfoin.dayly.z + S_solitaryBees_Sainfoin.dayly.z + (1 | Id.fac / Plant_Id.fac)"),
  data = df8,
  family = binomial()
))

summary(mod.2_Sainfoin) #For the coefficient plots

mod.2b_Sainfoin <- do.call(glmer, list(
  formula = as.formula("cbind(n_flowers_with_fruits, n_flowers_without_fruits) ~ S_allBees_Sainfoin.dayly.z + (1 | Id.fac / Plant_Id.fac)"),
  data = df8,
  family = binomial()
))

summary(mod.2b_Sainfoin) #All together significant: for the effect plotsS

r.squaredGLMM(mod.2_Sainfoin)
r.squaredGLMM(mod.2b_Sainfoin)

#### Model 2 ######################################################

library(car)
vif(mod.2_Sainfoin) #No extreme values: >5 would be problematic

#Assessing model assumptions: (check page 144)
par(mfrow=c(2,3))

#fitted vs. residuals
scatter.smooth(fitted(mod.2_Sainfoin), resid(mod.2_Sainfoin)) #this looks quite good (no strong positive correlation = no strong shrinkage)
abline(h=0, lty=2, col="red") #Some values do not fit the data well (very large residuals) -> for large seed sets the fitted values are not very reliable. 
title("Tukey-Anscombe Plot")

#qq-plot of residuals
qqnorm(resid(mod.2_Sainfoin)) #qq-plot of residuals
qqline(resid(mod.2_Sainfoin))

#qq-plot  of random effects: Id
qqnorm(ranef(mod.2_Sainfoin)$Id.fac[,1])
qqline(ranef(mod.2_Sainfoin)$Id.fac[,1])

#qq-plot of the random effects: Plant
qqnorm(ranef(mod.2_Sainfoin)$Plant_Id.fac[,1]) 
qqline(ranef(mod.2_Sainfoin)$Plant_Id.fac[,1])

#fitted vs. observed
plot(fitted(mod.2_Sainfoin), jitter(df8$n_flowers_with_fruits/(df8$n_flowers_with_fruits + df8$n_flowers_without_fruits),0.05))
abline(0,1) 

dev.off()  

#Heteroscedasticity (= nonhomogeneity of the residual variance)
par(mfrow=c(1,3))
scatter.smooth(df8$S_Apis_Sainfoin.dayly.z, resid(mod.2_Sainfoin)); abline(0,0, lty=2, col="red")
scatter.smooth(df8$S_Bombus_Sainfoin.dayly.z, resid(mod.2_Sainfoin)); abline(0,0, lty=2, col="red")
scatter.smooth(df8$S_solitaryBees_Sainfoin.dayly.z, resid(mod.2_Sainfoin)); abline(0,0, lty=2, col="red")
dev.off() #There are some outlyers

#Check if random effects are 0:
mean(ranef(mod.2_Sainfoin)$Id.fac[,1])
mean(ranef(mod.2_Sainfoin)$Plant_Id.fac[,1])

#Check for overdispersion:
dispersion_glmer(mod.2_Sainfoin) #without observation-level RF we have overdispersion.

#### Model 1b ######################################################

#Assessing model assumptions: (check page 144)
par(mfrow=c(2,3))

#fitted vs. residuals
scatter.smooth(fitted(mod.2b_Sainfoin), resid(mod.2b_Sainfoin)) #this looks quite good (no strong positive correlation = no strong shrinkage)
abline(h=0, lty=2, col="red") #Some values do not fit the data well (very large residuals) -> for large seed sets the fitted values are not very reliable. 
title("Tukey-Anscombe Plot")

#qq-plot of residuals
qqnorm(resid(mod.2b_Sainfoin)) #qq-plot of residuals
qqline(resid(mod.2b_Sainfoin))

#qq-plot  of random effects: Id
qqnorm(ranef(mod.2b_Sainfoin)$Id.fac[,1])
qqline(ranef(mod.2b_Sainfoin)$Id.fac[,1])

#qq-plot of the random effects: Plant
qqnorm(ranef(mod.2b_Sainfoin)$Plant_Id.fac[,1]) 
qqline(ranef(mod.2b_Sainfoin)$Plant_Id.fac[,1])

#fitted vs. observed
plot(fitted(mod.2b_Sainfoin), jitter(df8$n_flowers_with_fruits/(df8$n_flowers_with_fruits + df8$n_flowers_without_fruits),0.05))
abline(0,1) 

dev.off()  

#Heteroscedasticity (= nonhomogeneity of the residual variance)
par(mfrow=c(1,3))
scatter.smooth(df8$S_allBees_Sainfoin.dayly.z, resid(mod.2b_Sainfoin)); abline(0,0, lty=2, col="red")
dev.off() #

#Check if random effects are 0:
mean(ranef(mod.2b_Sainfoin)$Id.fac[,1])
mean(ranef(mod.2b_Sainfoin)$Plant_Id.fac[,1])

#Check for overdispersion:
dispersion_glmer(mod.2b_Sainfoin) #without observation-level RF we have overdispersion.

#################################################################################################################
#################################################################################################################
