#####################################################################
###
### 4e. Pollinator experiment: 
### Modelling plant reproductive success: comfrey fruit set ~ Urban
### Code by: David Frey and Merin Reji Chacko
### Last edited: 17.07.2025
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

df7e <- merge(datExpl, df6e)

str(df7e)
names(df7e)

#Kick out garden Nr. 39
df7e <- droplevels(df7e[which(df7e$Id != 39),])
dim(df7e)

#Transform factor variables
df7e$Id.fac <- as.factor(df7e$Id)
df7e$Plant_Id.fac <- as.factor(df7e$plant_id)
df7e$Branch_Id.fac <- as.factor(df7e$branch_id)

#Make orthogonal polynomials:
poly.A_Bombus_Comfrey.dayly <- poly(df7e$A_Bombus_Comfrey.dayly,2)
df7e$A_Bombus_Comfrey.dayly.1.z <- scale(poly.A_Bombus_Comfrey.dayly[,1])
df7e$A_Bombus_Comfrey.dayly.2.z <- scale(poly.A_Bombus_Comfrey.dayly[,2])

###################
## Run the model ##
###################

#Abundance model: a binomial glmer (failure vrs. success); with dayly rates!

mod.1_Comfrey <- glmer(cbind(n_flowers_with_seeds, n_flowers_without_seeds) ~ PlantS.z + Urban_500.z + PlantS.z:Urban_500.z + (1|Id.fac/Plant_Id.fac/Branch_Id.fac), data=df7e, family="binomial") #+ A_Bombus_Comfrey.dayly.2.z

r.squaredGLMM(mod.1_Comfrey)

#mod.1b_Comfrey <- glmer(cbind(n_flowers_with_seeds, n_flowers_without_seeds) ~ PlantS.z * Urban_500.z + (1|Id.fac/Plant_Id.fac), data=df7e, family="binomial")
#summary(mod.1b_Comfrey)

#Assessing model assumptions: (check page 144)
par(mfrow=c(2,3))

plot(fitted(mod.1_Comfrey), resid(mod.1_Comfrey))#fitted vs. observed values
abline(h=0)

qqnorm(resid(mod.1_Comfrey)) #qq-plot of residuals
qqline(resid(mod.1_Comfrey))

qqnorm(ranef(mod.1_Comfrey)$Id.fac[,1]) #qq-plot of the random effects: Id
qqline(ranef(mod.1_Comfrey)$Id.fac[,1])

qqnorm(ranef(mod.1_Comfrey)$Plant_Id.fac[,1]) #qq-plot of the random effects: Plant
qqline(ranef(mod.1_Comfrey)$Plant_Id.fac[,1])

qqnorm(ranef(mod.1_Comfrey)$Branch_Id.fac[,1]) #qq-plot of the random effects: Branch
qqline(ranef(mod.1_Comfrey)$Branch_Id.fac[,1])

plot(fitted(mod.1_Comfrey), jitter(df7e$n_flowers_with_seeds/(df7e$n_flowers_with_seeds + df7e$n_flowers_without_seeds),0.05))
abline(0,1)    

dev.off() 

#Heteroscedasticity (= nonhomogeneity of the residual variance)
scatter.smooth(df7e$A_Bombus_Comfrey.dayly.1.z, resid(mod.1_Comfrey)); abline(0,0, lty=2, col="red")
scatter.smooth(df7e$A_Bombus_Comfrey.dayly.2.z, resid(mod.1_Comfrey)); abline(0,0, lty=2, col="red")

#Check if random effects are 0:
mean(ranef(mod.1_Comfrey)$Id.fac[,1])
mean(ranef(mod.1_Comfrey)$Plant_Id.fac[,1])
mean(ranef(mod.1_Comfrey)$Branch_Id.fac[,1])

#Check for overdispersion:
dispersion_glmer(mod.1_Comfrey) #overdispersed without Branch_ID

#################################################################################################################
#################################################################################################################
