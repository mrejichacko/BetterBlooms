#####################################################################
#####################################################################
#####################################################################
###
### 4c. Pollinator experiment: 
### Modelling plant reproductive success: Radish seedset ~ Urban
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.01.2026
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

df7c <- merge(datExpl, df6c)

str(df7c)
names(df7c)

#Kick out garden Nr. 39
df7c <- droplevels(df7c[which(df7c$Id != 39),])
dim(df7c)

#Transform factor variables
df7c$Id.fac <- as.factor(df7c$Id)
df7c$Plant_Id.fac <- as.factor(df7c$plant_id)
df7c$Fruit_Id.fac <- as.factor(df7c$fruit_id)

###################
## Run the model ##
###################

## Here there are no plants to exclude. 

#Abundance model: a poisson GLMM; with dayly rates!

mod.1_Radish_seedset <- glmer(n_seeds ~ PlantS.z * Urban_500.z +  (1|Id.fac/Plant_Id.fac), data=df7c, family="poisson")
summary(mod.1_Radish_seedset) 

r.squaredGLMM(mod.1_Radish_seedset)

library(car)
vif(mod.1_Radish_seedset) #No extreme values: >5 would be problematic

#Assessing model assumptions: (check page 144)
par(mfrow=c(2,3))

#fitted vs. residuals
scatter.smooth(fitted(mod.1_Radish_seedset), resid(mod.1_Radish_seedset)) #
abline(h=0, lty=2, col="red") # 
title("Tukey-Anscombe Plot")

#qq-plot of residuals
qqnorm(resid(mod.1_Radish_seedset))
qqline(resid(mod.1_Radish_seedset)) #Most values fit!

#qq-plot  of random effects: Id
qqnorm(ranef(mod.1_Radish_seedset)$Id.fac[,1])
qqline(ranef(mod.1_Radish_seedset)$Id.fac[,1]) #

#qq-plot  of random effects: Plant_Id
qqnorm(ranef(mod.1_Radish_seedset)$Plant_Id.fac[,1])
qqline(ranef(mod.1_Radish_seedset)$Plant_Id.fac[,1]) 

#qq-plot  of random effects: Fruit_Id #-> no additional variance explained!
#qqnorm(ranef(mod.1_Radish_seedset)$Fruit_Id.fac[,1])
#qqline(ranef(mod.1_Radish_seedset)$Fruit_Id.fac[,1]) 

#fitted vs. observed
scatter.smooth(fitted(mod.1_Radish_seedset), df7c$N_seeds)#data vs. fitted

dev.off()  #We have a couple of outlyers (the plants with <100 seeds)

#Heteroscedasticity (= nonhomogeneity of the residual variance)
par(mfrow=c(2,3))
scatter.smooth(df7c$A_Apis_Radish.dayly.z, resid(mod.1_Radish_seedset)); abline(0,0, lty=2, col="red")
scatter.smooth(df7c$A_socialBees_Radish.dayly.z, resid(mod.1_Radish_seedset)); abline(0,0, lty=2, col="red")
scatter.smooth(df7c$A_solitaryBees_Radish.dayly.z, resid(mod.1_Radish_seedset)); abline(0,0, lty=2, col="red")
scatter.smooth(df7c$A_otherAculeata_Radish.dayly.z, resid(mod.1_Radish_seedset)); abline(0,0, lty=2, col="red")
scatter.smooth(df7c$A_Syrphidae_Radish.dayly.z, resid(mod.1_Radish_seedset)); abline(0,0, lty=2, col="red")
scatter.smooth(df7c$A_Coleoptera_Radish.dayly.z, resid(mod.1_Radish_seedset)); abline(0,0, lty=2, col="red") 
dev.off()

#Check if random effects are 0:
mean(ranef(mod.1_Radish_seedset)$Id.fac[,1])
mean(ranef(mod.1_Radish_seedset)$Plant_Id.fac[,1])

#Check for overdispersion:
dispersion_glmer(mod.1_Radish_seedset) #no overdispersion, a bit underdispersion

# check for spatial autocorrelation

library(DHARMa)

sim_radish <- simulateResiduals(mod.1_Radish_seedset)

# Aggregate to unique spatial locations (garden = Id.fac)
sim_radish_garden <- recalculateResiduals(sim_radish, group = df7c$Id.fac)

# (unique x/y per garden)
sp_test_radish <- testSpatialAutocorrelation(
  sim_radish_garden,
  x = df7c$X_KOORDINATE[!duplicated(df7c$Id.fac)],
  y = df7c$Y_KOORDINATE[!duplicated(df7c$Id.fac)]
)

sp_test_radish
sp_test_radish$p.value
sp_test_radish$statistic

table_s4 <- data.frame(Predictor = "Habitat",
                       Response = "Radish seed set", 
                       Observed = unname(sp_test_radish$statistic[1]),
                       Expected = unname(sp_test_radish$statistic[2]),
                       SD = unname(sp_test_radish$statistic[3]),
                       P_value = sp_test_radish$p.value
)

table_s4

write.table(
  table_s4,
  "results/Table_S4_MoransI_DHARMa_seedset_abundance_models.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,  # IMPORTANT: do not rewrite the header
  append = TRUE
)


#################################################################################################################
#################################################################################################################
