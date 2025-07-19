#####################################################################
###
### 4f. Pollinator experiment: 
### Modelling plant reproductive success: comfrey seed set ~ Urban
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

df7f <- merge(datExpl, df6f)

str(df7f)
names(df7f)

#Kick out garden Nr. 39
df7f <- droplevels(df7f[which(df7f$Id != 39),])
dim(df7f)

#Transform factor variables
df7f$Id.fac <- as.factor(df7f$Id)
df7f$Plant_Id.fac <- as.factor(df7f$plant_id)
df7f$Branch_Id.fac <- as.factor(df7f$branch_id)
df7f$Flower_Id.fac <- as.factor(df7f$flower_id)

#Make orthogonal polynomials:
poly.A_Bombus_Comfrey.dayly <- poly(df7f$A_Bombus_Comfrey.dayly,2)
df7f$A_Bombus_Comfrey.dayly.1.z <- scale(poly.A_Bombus_Comfrey.dayly[,1])
df7f$A_Bombus_Comfrey.dayly.2.z <- scale(poly.A_Bombus_Comfrey.dayly[,2])

###################
## Run the model ##
###################

#Abundance model: a binomial glmer (failure vrs. success); with dayly rates!

#Make a failure variable to use a binomial model: N_unfertilized_ovules
#The maximum number of seeds per flower is four; if N seeds = 4 -> 100% success = 0 failure
# -> note we have a fixed number of trials (4)

mod.1_Comfrey_seedset <- glmer(cbind(n_seeds, n_unfertilized_ovules) ~ PlantS.z * Urban_500.z  + (1|Id.fac/Plant_Id.fac), data=df7f, family="binomial", control = glmerControl(optimizer = "bobyqa")) #A_Bombus_Comfrey.dayly.2.z 
summary(mod.1_Comfrey_seedset) #

r.squaredGLMM(mod.1_Comfrey_seedset)

#mod.1b_Comfrey_seedset <- lmer(N_seeds/4 ~ A_Bombus_Comfrey.dayly.1.z + A_Bombus_Comfrey.dayly.2.z + (1|Id.fac/Plant_Id.fac), data=df7f)
#summary(mod.1b_Comfrey_seedset)

#mod.1c_Comfrey_seedset <- lmer(N_seeds/N_flowers ~ A_Bombus_Comfrey.dayly.1.z + A_Bombus_Comfrey.dayly.2.z + (1|Id.fac/Plant_Id.fac), data=df7f)
#summary(mod.1c_Comfrey_seedset)

#Assessing model assumptions: (check page 144)
par(mfrow=c(2,3))

plot(fitted(mod.1_Comfrey_seedset), resid(mod.1_Comfrey_seedset))
abline(h=0)

qqnorm(resid(mod.1_Comfrey_seedset)) #qq-plot of residuals
qqline(resid(mod.1_Comfrey_seedset))

qqnorm(ranef(mod.1_Comfrey_seedset)$Id.fac[,1]) #qq-plot of the random effects: Id
qqline(ranef(mod.1_Comfrey_seedset)$Id.fac[,1])

qqnorm(ranef(mod.1_Comfrey_seedset)$Plant_Id.fac[,1]) #qq-plot of the random effects: Plant
qqline(ranef(mod.1_Comfrey_seedset)$Plant_Id.fac[,1])

plot(fitted(mod.1_Comfrey_seedset), jitter(df7f$N_seeds/4,0.05))
abline(0,1)    

dev.off() 

#Heteroscedasticity (= nonhomogeneity of the residual variance)
scatter.smooth(df7f$A_Bombus_Comfrey.dayly.1.z, resid(mod.1_Comfrey_seedset)); abline(0,0, lty=2, col="red")
#scatter.smooth(df7f$A_Bombus_Comfrey.dayly.2.z, resid(mod.1_Comfrey_seedset)); abline(0,0, lty=2, col="red")

#Check if random effects are 0:
mean(ranef(mod.1_Comfrey_seedset)$Id.fac[,1])
mean(ranef(mod.1_Comfrey_seedset)$Plant_Id.fac[,1])

#Check for overdispersion:
dispersion_glmer(mod.1_Comfrey_seedset) # underdispersed...

#######################################################################################################################################
#######################################################################################################################################