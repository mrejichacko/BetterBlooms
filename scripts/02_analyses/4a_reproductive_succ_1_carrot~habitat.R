#####################################################################
#####################################################################
#####################################################################
###
### 4a. Pollinator experiment: 
### Modelling plant reproductive success: Daucus ~ Urban
### Code by: David Frey and Merin Reji Chacko
### Last edited: 17.07.2025
#####################################################################
#####################################################################
#####################################################################
###
### Prep
library(lme4)
library(blmeco)
library(MuMIn)

#################### RUN 3a or load from environment!!! #############

rm(list=ls())
gc()

load("environments/3a_environment.RData")

########################
# Combine the datasets #
########################

df7a <- merge(datExpl, df6a)

str(df7a)
names(df7a)

#Kick out garden Nr. 39
df7a <- droplevels(df7a[which(df7a$Id != 39),])
dim(df7a)

#Transform factor variables
df7a$Id.fac <- as.factor(df7a$Id)
df7a$Plant_Id.fac <- as.factor(df7a$plant_id)
df7a$Umbell_Id.fac <- as.factor(df7a$umbel_id)

###################
## Run the model ##
###################

#GLMM-version (poisson)
mod.1_Carrot <- glmer(n_seeds ~ PlantS.z + Urban_500.z + PlantS.z:Urban_500.z + (1|Id.fac/Plant_Id.fac/Umbell_Id.fac), data=df7a, family="poisson")
summary(mod.1_Carrot)

r.squaredGLMM(mod.1_Carrot)

#Calculate the variance inflation factor (due to correlated predictors) since all the abundances are a bit correlated:
#http://www.sthda.com/english/articles/39-regression-model-diagnostics/160-multicollinearity-essentials-and-vif-in-r/
library(car)
vif(mod.1_Carrot) #No extreme values: >5 would be problematic

#Assessing model assumptions:check page 149 for poisson models
par(mfrow=c(2,3))

#fitted vs. residuals
scatter.smooth(fitted(mod.1_Carrot), resid(mod.1_Carrot)) #this looks quite good (no strong positive correlation = no strong shrinkage)
abline(h=0, lty=2, col="red") #Some values do not fit the data well (very large residuals) -> for large seed sets the fitted values are not very reliable. 
title("Tukey-Anscombe Plot")

#qq-plot of residuals
qqnorm(resid(mod.1_Carrot))
qqline(resid(mod.1_Carrot)) #Most values fit!

#qq-plot  of random effects: Id
qqnorm(ranef(mod.1_Carrot)$Id.fac[,1])
qqline(ranef(mod.1_Carrot)$Id.fac[,1]) #there is no serious deviation in the distribution of rf from the normal distribution

#qq-plot  of random effects: Plant_Id
qqnorm(ranef(mod.1_Carrot)$Plant_Id.fac[,1])
qqline(ranef(mod.1_Carrot)$Plant_Id.fac[,1]) ##there is no serious deviation in the distribution of rf from the normal distribution

#qq-plot  of random effects: Umbell_Id
qqnorm(ranef(mod.1_Carrot)$Umbell_Id.fac[,1])
qqline(ranef(mod.1_Carrot)$Umbell_Id.fac[,1]) ##there is no serious deviation in the distribution of rf from the normal distribution

#fitted vs. observed
scatter.smooth(fitted(mod.1_Carrot), df7a$N_seeds)#data vs. fitted -> it is clear where the big residuals come from. 

dev.off()

#Heteroscedasticity (= nonhomogeneity of the residual variance)
par(mfrow=c(1,2))
scatter.smooth(df7a$PlantS.z, resid(mod.1_Carrot)); abline(0,0, lty=2, col="red")
scatter.smooth(df7a$Urban_500.z, resid(mod.1_Carrot)); abline(0,0, lty=2, col="red")
dev.off()

#Chec for overdispersion
dispersion_glmer(mod.1_Carrot)

#Check if random effects are 0:
mean(ranef(mod.1_Carrot)$Id.fac[,1])
mean(ranef(mod.1_Carrot)$Plant_Id.fac[,1])
mean(ranef(mod.1_Carrot)$Umbell_Id.fac[,1])

#################################################################################################################
#################################################################################################################