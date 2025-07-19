#####################################################################
#####################################################################
#####################################################################
###
### 2. Pollinator experiment: 
### Trait values ~ Garden + Urban: Graphical representation + Models
### Code by: David Frey and Merin Reji Chacko
### Last edited: 19.07.2025
#####################################################################
#####################################################################
#####################################################################
###
### Prep

rm(list=ls())
gc()

#Load the pollinator data
df1 <-read.csv("raw_data/RejiChacko_etal_2025_EnviDat/06_trait_data/individual_traits.csv", header = TRUE, sep = ",")
df1$capture_period <- as.Date(df1$capture_date, format = "%d.%m.%Y")
df2 <-readxl::read_xlsx("raw_data/RejiChacko_etal_2025_EnviDat/05_field_data/raw_sampling_data.xlsx")
df2$capture_period <- as.Date(df2$capture_period)
df3 <- merge(df1, df2[,c("Id","capture_period","capture_window","wind_speed","cloudiness")], by=c("Id","capture_period","capture_window"))
df4<-read.table("raw_data/plant_floristic_data.txt", header = TRUE, sep = ";")
df5<-read.table("raw_data/explanatory_variables.txt", header = TRUE, sep = ";")

#Combine the datasets: predictors
df6 <- Reduce(function(x, y) merge(x, y, by="Id", sort=TRUE), list(df3,df4[,c(1,5)],df5[,c(1,4,5,13:16)]))

str(df6)
head(df6)
names(df6)

#Rename the plant variable
colnames(df6)[colnames(df6)=="SR_all_insect_pollinated_May_August"]<-"PlantS"

#Kick out garden Nr. 39 and bad samples, and nectar robbers!
df6 <- droplevels(subset(df6, Id != 39 & comment != "not a hoverfly" & comment != "sample lost" & nectar_robber != 1)) 
dim(df6)

###########################
### Scale the variables ###
###########################

df6$Id.fac <- as.factor(df6$Id)

##Landscape, all scales:

poly.Urban_500 <- poly(df6$Urban_500,2)
df6$Urban_500.1.z <- scale(poly.Urban_500[,1])
df6$Urban_500.2.z <- scale(poly.Urban_500[,2])

poly.Urban_250 <- poly(df6$Urban_250,2)
df6$Urban_250.1.z <- scale(poly.Urban_250[,1])
df6$Urban_250.2.z <- scale(poly.Urban_250[,2])

poly.Urban_100 <- poly(df6$Urban_100,2)
df6$Urban_100.1.z <- scale(poly.Urban_100[,1])
df6$Urban_100.2.z <- scale(poly.Urban_100[,2])

poly.Urban_50 <- poly(df6$Urban_50,2)
df6$Urban_50.1.z <- scale(poly.Urban_50[,1])
df6$Urban_50.2.z <- scale(poly.Urban_50[,2])

##Garden:

#Plants species richness May-August
poly.PlantS <- poly(df6$PlantS, 2)
df6$PlantS.1.z <- scale(poly.PlantS[,1])
df6$PlantS.2.z <- scale(poly.PlantS[,2])

##Variables related to the capturing:

df6$capture_window.z <-scale(df6$capture_window)
df6$capture_window.ord.fac <-as.ordered(df6$capture_window)

table(df6$capture_window.ord.fac)

#Check the date:
summary(df6$capture_period)
str(df6$capture_period)

df6$capture_periodMMDD <- as.Date(df6$capture_period, tryFormats = c("%d.%m.%Y"))
#df6$capture_periodMMDD.ord.fac <- as.ordered(df6$capture_periodMMDD)

#cloudiness: rescale in proportions:

df6$cloudiness[df6$cloudiness==1]<-0
df6$cloudiness[df6$cloudiness==2]<-0.125
df6$cloudiness[df6$cloudiness==3]<-0.25
df6$cloudiness[df6$cloudiness==4]<-0.375

df6$cloudiness.z <- scale(df6$cloudiness)
table(df6$cloudiness)

#Beaufort-Scale: an ordered factor
df6$wind_speed.z <- scale(df6$wind_speed)
df6$wind_speed.ord.fac <- as.ordered(df6$wind_speed)

summary(df6$wind_speed.ord.fac)

#Create the weighted tongue data:
df6$PRB.wgt <- df6$proboscis_length/df6$intertegular_distance

names(df6)

#Rename ITD
colnames(df6)[colnames(df6)=="intertegular_distance"]<-"ITD"
colnames(df6)[colnames(df6)=="forewing_length"]<-"FWL"
colnames(df6)[colnames(df6)=="proboscis_length"]<-"PRL"
colnames(df6)[colnames(df6)=="labellum_prementum_ratio"]<-"LPR"

#Make a sex variable

df6$sex <- df6$male
df6$sex[df6$sex==1] <- "male"
df6$sex[df6$sex==0] <- "female"

df6$sex <- as.factor(df6$sex)
summary(df6$sex)


#Check the variables
hist(df6$ITD) #-> its bimodal!
hist(log(df6$ITD))

#ITD is bimodally distributed
hist(df6[which(df6$ITD < 1.648721),]$ITD)
hist(df6[which(df6$ITD > 1.648721),]$ITD)

checklist <- read.csv("raw_data/RejiChacko_etal_2025_EnviDat/04_taxonomic_data/taxa_checklist.csv")
df6 <- merge(df6, checklist, by = "taxon")
rm(checklist)

hist(df6[which(df6$pollinator_group == "Anthophila" & df6$ITD > 1.648721),]$ITD)
hist(df6[which(df6$pollinator_group == "Anthophila" & df6$ITD < 1.648721),]$ITD)

hist(df6[which(df6$ITD < 2),]$ITD)
hist(df6[which(df6$ITD < 2),]$ITD)

hist(df6[which(df6$pollinator_group == "Anthophila" & df6$ITD >= 2),]$ITD)
hist(df6[which(df6$pollinator_group == "Anthophila" & df6$ITD < 2),]$ITD)

summary(df6$ITD)
sd(df6$ITD)

hist(df6$PRB.wgt)
hist(log(df6$PRB.wgt))

names(df6)

###################################
### Happy modelling: Abundance~ ###
###################################

library(lme4)
library(blmeco)
library(plyr)
library(AICcmodavg)
library(piecewiseSEM)
#Useful discussion for R2: https://stats.stackexchange.com/questions/277371/are-r2-for-glmm-useful-for-modelers-but-not-necessarily-for-readers
library(MuMIn)

#Check for NA's:

summary(df6$capture_window.ord.fac) #no
summary(df6$wind_speed.ord.fac) #no
summary(df6$cloudiness.z)#no
summary(df6$sex) #two NA's
summary(df6$phytometer_plant) #no

summary(df6$PRB.wgt)#~700 NA's

df7 <- droplevels(df6[!is.na(df6$sex),])
df8 <- droplevels(df7[which(df7$capture_window !=10),])

########################
### The trait models ###
########################

df8$PRB.wgt <- df8$PRL/df8$ITD
df8 <- subset(df8, !is.na(PRB.wgt))
#log(Proboscis length) and log(LPR): capture_window is random effect
mod.1 <- lmer(log(PRB.wgt) ~ PlantS.1.z * Urban_500.1.z + PlantS.1.z * phytometer_plant + sex + (1|Id.fac) + (1|capture_periodMMDD) + (1|capture_window.ord.fac)+ (1|family/genus/species), data=df8, REML = F, subset = pollinator_group =="Anthophila")
mod.2 <- lmer(log(PRB.wgt) ~ PlantS.1.z * Urban_500.1.z + PlantS.1.z * phytometer_plant + sex + (1|Id.fac) + (1|capture_periodMMDD) + (1|capture_window.ord.fac)+ (1|family/genus/species), data=df8[which(df8$genus!="Apis"),], REML = F, subset = pollinator_group =="Anthophila")
mod.3 <- lmer(log(PRB.wgt) ~ PlantS.1.z * Urban_500.1.z + PlantS.1.z * phytometer_plant + sex + (1|Id.fac) + (1|capture_periodMMDD) + (1|capture_window.ord.fac)+ (1|family/genus/species), data=df8[which(df8$genus!="Apis"),], REML = F, subset = pollinator_group_sociality =="social_Bees")
mod.4 <- lmer(log(PRB.wgt+0.00001) ~ PlantS.1.z * Urban_500.1.z + PlantS.1.z * phytometer_plant + sex + (1|Id.fac) + (1|capture_periodMMDD) + (1|capture_window.ord.fac)+ (1|family/genus/species), data=df8, REML = F, subset = pollinator_group_sociality=="solitary_Bees")
mod.5 <- lmer(log(LPR) ~ PlantS.1.z * Urban_500.1.z + PlantS.1.z * phytometer_plant + sex + (1|Id.fac) + (1|capture_periodMMDD) + (1|capture_window.ord.fac)+ (1|genus/species), data=df8, REML = F, subset = family =="Syrphidae")

#log(ITD): capture_window is fixed effect
mod.6 <- lmer(log(ITD) ~ PlantS.1.z * Urban_500.1.z + PlantS.1.z * phytometer_plant + capture_window.ord.fac + wind_speed.ord.fac + cloudiness.z + sex + (1|Id.fac) + (1|capture_periodMMDD) + (1|family/genus/species), data=df8, REML = F, subset = pollinator_group =="Anthophila")
mod.7 <- lmer(log(ITD)~ PlantS.1.z * Urban_500.1.z + PlantS.1.z * phytometer_plant + capture_window.ord.fac + wind_speed.ord.fac + cloudiness.z + sex + (1|Id.fac) + (1|capture_periodMMDD)+ (1|family/genus/species), data=df8[which(df8$genus!="Apis"),], REML = F, subset = pollinator_group =="Anthophila")
mod.8 <- lmer(log(ITD) ~ PlantS.1.z * Urban_500.1.z + PlantS.1.z * phytometer_plant + capture_window.ord.fac + wind_speed.ord.fac + cloudiness.z + sex + (1|Id.fac) + (1|capture_periodMMDD)+ (1|family/genus/species), data=df8[which(df8$genus!="Apis"),], REML = F, subset = pollinator_group_sociality =="social_Bees")
mod.9 <- lmer(log(ITD) ~ PlantS.1.z * Urban_500.1.z + PlantS.1.z * phytometer_plant + capture_window.ord.fac + wind_speed.ord.fac + cloudiness.z + sex + (1|Id.fac) + (1|capture_periodMMDD)+ (1|family/genus/species), data=df8, REML = F, subset = pollinator_group_sociality=="solitary_Bees")
mod.10 <- lmer(log(FWL) ~ PlantS.1.z * Urban_500.1.z + PlantS.1.z * phytometer_plant + capture_window.ord.fac + wind_speed.ord.fac + cloudiness.z + sex + (1|Id.fac) + (1|capture_periodMMDD) + (1|genus/species), data=df8, REML = F, subset = family =="Syrphidae")

modlist.1 <- list(mod.1=mod.1,
                  mod.2=mod.2,
                  mod.3=mod.3,
                  mod.4=mod.4,
                  mod.5=mod.5,
                  mod.6=mod.6,
                  mod.7=mod.7,
                  mod.8=mod.8,
                  mod.9=mod.9,
                  mod.10=mod.10)#

library(plyr)

modlist.1_Fit <- ldply(modlist.1, function(mod) {data.frame(R2m = round(r.squaredGLMM(mod)[1],3),
                                                            R2c = round(r.squaredGLMM(mod)[2],3),
                                                            meanRF_Id.fac = round(mean(ranef(mod)$Id.fac[,1]),5),
                                                            meanRF_FangPer = round(mean(ranef(mod)$capture_periodMMDD[,1]),5))})
modlist.1_Fit

#How to deal with phylogenetic corrections, see for instance:
#https://bbolker.github.io/mixedmodels-misc/notes/phylog.html

summary(mod.1)
rsquared(mod.1)
summary(mod.2)
summary(mod.3)
summary(mod.4)
summary(mod.5)

summary(mod.6)
summary(mod.7)
summary(mod.8)
summary(mod.9)
summary(mod.10)

#########################################################
##################### Tongue models #####################
#########################################################
#dev.off()

#setwd("C:/Users/David/Desktop/Pollinator Experiment/Analysis/Plots/Plots - TraitFiltering - ResidualAnalysis")

clab = 1
cmain = 1
caxis = 0.5

#Mod.1
#tiff(filename = "data/PollinatorExperiment2025/Figures/20181230_mod.1_TraitFilter_resid.tif",width = 20, height = 20, units = "cm", res = 250, compression = c("none"))
par(mfrow=c(4,4))
plot(fitted(mod.1), resid(mod.1), main="(A) Residuals vs. fitted", xlab = "Fitted values", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.1), main="Q-Q plot residuals", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(resid(mod.1))
qqnorm(ranef(mod.1)$Id.fac[,1], main="Q-Q plot random \neffect: Garden", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.1)$Id.fac[,1]) #
qqnorm(ranef(mod.1)$capture_periodMMDD[,1], main="Q-Q plot random \neffect: Date ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.1)$capture_periodMMDD[,1])
qqnorm(ranef(mod.1)$capture_window.ord.fac[,1], main="Q-Q plot random \neffect: Time window ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.1)$capture_window.ord.fac[,1])
qqnorm(ranef(mod.1)$family[,1], main="Q-Q plot random \neffect: family ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.1)$family[,1])
qqnorm(ranef(mod.1)$genus[,1], main="Q-Q plot random \neffect: Genus ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.1)$genus[,1])
qqnorm(ranef(mod.1)$species[,1], main="Q-Q plot random \neffect: Species ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.1)$species[,1])
plot(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$PRB.wgt)),]$PlantS.1.z, resid(mod.1), main="Residuals vs. \nPredictor 1", xlab="Garden plant \nspecies richness S", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$PRB.wgt)),]$Urban_500.1.z, resid(mod.1), main="Residuals vs. \nPredictor 2", xlab= "Urban intensity", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$PRB.wgt)),]$phytometer_plant, resid(mod.1), main="Residuals vs. \nPredictor 3", xlab= "phytometer_plant", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$PRB.wgt)),]$sex, resid(mod.1), main="Residuals vs. \nPredictor 4", xlab= "sex", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.1), log(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$PRB.wgt)),]$PRB.wgt), main="Observed vs. \nfitted", xlab = "Fitted values", ylab = "Observed values", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted


#Mod.2
#tiff(filename = "data/PollinatorExperiment2025/Figures/20181230_mod.2_TraitFilter_resid.tif",width = 20, height = 20, units = "cm", res = 250, compression = c("none"))
par(mfrow=c(4,4))
plot(fitted(mod.2), resid(mod.2), main="(B) Residuals vs. fitted", xlab = "Fitted values", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.2), main="Q-Q plot residuals", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(resid(mod.2))
qqnorm(ranef(mod.2)$Id.fac[,1], main="Q-Q plot random \neffect: Garden", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.2)$Id.fac[,1]) #
qqnorm(ranef(mod.2)$capture_periodMMDD[,1], main="Q-Q plot random \neffect: Date ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.2)$capture_periodMMDD[,1]) #
qqnorm(ranef(mod.2)$capture_window.ord.fac[,1], main="Q-Q plot random \neffect: Time window ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.2)$capture_window.ord.fac[,1])
qqnorm(ranef(mod.2)$family[,1], main="Q-Q plot random \neffect: family ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.2)$family[,1])
qqnorm(ranef(mod.2)$genus[,1], main="Q-Q plot random \neffect: Genus ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.2)$genus[,1])
qqnorm(ranef(mod.2)$species[,1], main="Q-Q plot random \neffect: Species ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.2)$species[,1])
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & !is.na(df8$PRB.wgt)),]$PlantS.1.z, resid(mod.2), main="Residuals vs. \nPredictor 1", xlab="Garden plant species S", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & !is.na(df8$PRB.wgt)),]$Urban_500.1.z, resid(mod.2), main="Residuals vs. \nPredictor 2", xlab= "Urban intensity", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & !is.na(df8$PRB.wgt)),]$phytometer_plant, resid(mod.2), main="Residuals vs.\nPredictor 2", xlab= "phytometer_plant species", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & !is.na(df8$PRB.wgt)),]$sex, resid(mod.2), main="Residuals vs. \nPredictor 7", xlab= "sex", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.2), log(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & !is.na(df8$PRB.wgt)),]$PRB.wgt), main="Observed vs. fitted", xlab = "Fitted values", ylab = "Observed values", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted


#Mod.3
#tiff(filename = "data/PollinatorExperiment2025/Figures/20181230_mod.3_TraitFilter_resid.tif",width = 20, height = 20, units = "cm", res = 250, compression = c("none"))
par(mfrow=c(4,4))
plot(fitted(mod.3), resid(mod.3), main="(C) Residuals vs. fitted", xlab = "Fitted values", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.3), main="Q-Q plot residuals", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(resid(mod.3))
qqnorm(ranef(mod.3)$Id.fac[,1], main="Q-Q plot random \neffect: Garden", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.3)$Id.fac[,1]) #
qqnorm(ranef(mod.3)$capture_periodMMDD[,1], main="Q-Q plot random \neffect: Date ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.3)$capture_periodMMDD[,1]) #
qqnorm(ranef(mod.3)$capture_window.ord.fac[,1], main="Q-Q plot random \neffect: Time window ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.3)$capture_window.ord.fac[,1])
qqnorm(ranef(mod.3)$family[,1], main="Q-Q plot random \neffect: family ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.3)$family[,1])
qqnorm(ranef(mod.3)$genus[,1], main="Q-Q plot random \neffect: Genus ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.3)$genus[,1])
qqnorm(ranef(mod.3)$species[,1], main="Q-Q plot random \neffect: Species ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.3)$species[,1])
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & df8$pollinator_group_sociality == "social_Bees" & !is.na(df8$PRB.wgt)),]$PlantS.1.z, resid(mod.3), main="Residuals vs. \nPredictor 1", xlab="Garden plant species S", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & df8$pollinator_group_sociality == "social_Bees" & !is.na(df8$PRB.wgt)),]$Urban_500.1.z, resid(mod.3), main="Residuals vs. \nPredictor 2", xlab= "Urban intensity", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & df8$pollinator_group_sociality == "social_Bees" & !is.na(df8$PRB.wgt)),]$phytometer_plant, resid(mod.3), main="Residuals vs. \nPredictor 2", xlab= "phytometer_plant species", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & df8$pollinator_group_sociality == "social_Bees" & !is.na(df8$PRB.wgt)),]$sex, resid(mod.3), main="Residuals vs. \nPredictor 7", xlab= "sex", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.3), log(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & df8$pollinator_group_sociality == "social_Bees" & !is.na(df8$PRB.wgt)),]$PRB.wgt), main="Observed vs. fitted", xlab = "Fitted values", ylab = "Observed values", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted


#Mod.4
#tiff(filename = "data/PollinatorExperiment2025/Figures/20181230_mod.4_TraitFilter_resid.tif",width = 20, height = 20, units = "cm", res = 250, compression = c("none"))
par(mfrow=c(4,4))
plot(fitted(mod.4), resid(mod.4), main="(D) Residuals vs. fitted", xlab = "Fitted values", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.4), main="Q-Q plot residuals", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(resid(mod.4))
qqnorm(ranef(mod.4)$Id.fac[,1], main="Q-Q plot random \neffect: Garden", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.4)$Id.fac[,1]) #
qqnorm(ranef(mod.4)$capture_periodMMDD[,1], main="Q-Q plot random \neffect: Date ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.4)$capture_periodMMDD[,1]) #
qqnorm(ranef(mod.4)$capture_window.ord.fac[,1], main="Q-Q plot random \neffect: Time window ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.4)$capture_window.ord.fac[,1])
qqnorm(ranef(mod.4)$family[,1], main="Q-Q plot random \neffect: family ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.4)$family[,1])
qqnorm(ranef(mod.4)$genus[,1], main="Q-Q plot random \neffect: Genus ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.4)$genus[,1])
qqnorm(ranef(mod.4)$species[,1], main="Q-Q plot random \neffect: Species ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.4)$species[,1])
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$PRB.wgt)),]$PlantS.1.z, resid(mod.4), main="Residuals vs. \nPredictor 1", xlab="Garden plant species S", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$PRB.wgt)),]$Urban_500.1.z, resid(mod.4), main="Residuals vs. \nPredictor 2", xlab= "Urban intensity", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$PRB.wgt)),]$phytometer_plant, resid(mod.4), main="Residuals vs. \nPredictor 2", xlab= "phytometer_plant species", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$PRB.wgt)),]$sex, resid(mod.4), main="Residuals vs. \nPredictor 7", xlab= "sex", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.4), log(df8[which(df8$pollinator_group == "Anthophila" & df8$genus != "Apis" & df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$PRB.wgt)),]$PRB.wgt), main="Observed vs. fitted", xlab = "Fitted values", ylab = "Observed values", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

#Mod.5 
#tiff(filename = "data/PollinatorExperiment2025/Figures/20181230_mod.5_TraitFilter_resid.tif",width = 20, height = 20, units = "cm", res = 250, compression = c("none"))
par(mfrow=c(4,4))
plot(fitted(mod.5), resid(mod.5), main="(E) Residuals vs. fitted", xlab = "Fitted values", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.5), main="Q-Q plot residuals", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(resid(mod.5))
qqnorm(ranef(mod.5)$Id.fac[,1], main="Q-Q plot random \neffect: Garden", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.5)$Id.fac[,1]) #
qqnorm(ranef(mod.5)$capture_periodMMDD[,1], main="Q-Q plot random \neffect: Date ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.5)$capture_periodMMDD[,1])
qqnorm(ranef(mod.5)$capture_window.ord.fac[,1], main="Q-Q plot random \neffect: Time window ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.5)$capture_window.ord.fac[,1])
qqnorm(ranef(mod.5)$genus[,1], main="Q-Q plot random \neffect: Genus ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.5)$genus[,1])
qqnorm(ranef(mod.5)$species[,1], main="Q-Q plot random \neffect: Species ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.5)$species[,1])
plot(df8[which(df8$family == "Syrphidae" & !is.na(df8$LPR)),]$PlantS.1.z, resid(mod.5), main="Residuals vs. \nPredictor 1", xlab="Garden plant species S", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$family == "Syrphidae" & !is.na(df8$LPR)),]$Urban_500.1.z, resid(mod.5), main="Residuals vs. \nPredictor 2", xlab= "Urban intensity", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$family == "Syrphidae" & !is.na(df8$LPR)),]$phytometer_plant, resid(mod.5), main="Residuals vs. \nPredictor 3", xlab= "phytometer_plant species", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.5), log(df8[which(df8$family == "Syrphidae" & !is.na(df8$LPR)),]$LPR), main="Observed vs. \nfitted", xlab = "Fitted values", ylab = "Observed values", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

#########################################################
##################### ITD models ########################
#########################################################

#Mod.6 :
#tiff(filename = "data/PollinatorExperiment2025/Figures/20181230_mod.6_TraitFilter_resid.tif",width = 20, height = 20, units = "cm", res = 250, compression = c("none"))
par(mfrow=c(4,4))
plot(fitted(mod.6), resid(mod.6), main="(A) Residuals vs. fitted", xlab = "Fitted values", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.6), main="Q-Q plot residuals", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(resid(mod.6))
qqnorm(ranef(mod.6)$Id.fac[,1], main="Q-Q plot random \neffect: Garden", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.6)$Id.fac[,1]) #
qqnorm(ranef(mod.6)$capture_periodMMDD[,1], main="Q-Q plot random \neffect: Date ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.6)$capture_periodMMDD[,1])
qqnorm(ranef(mod.6)$family[,1], main="Q-Q plot random \neffect: family", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.6)$family[,1])
qqnorm(ranef(mod.6)$genus[,1], main="Q-Q plot random \neffect: Genus", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.6)$genus[,1])
qqnorm(ranef(mod.6)$species[,1], main="Q-Q plot random \neffect: Species", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.6)$species[,1])
plot(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$ITD)),]$PlantS.1.z, resid(mod.6), main="Residuals vs. \nPredictor 1", xlab="Garden plant species S", ylab="Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$ITD)),]$Urban_500.1.z, resid(mod.6), main="Residuals vs. \nPredictor 2", xlab= "Urban intensity", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$ITD)),]$phytometer_plant, resid(mod.6), main="Residuals vs. \nPredictor 3", xlab= "phytometer_plant species", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$ITD)),]$capture_window.ord.fac, resid(mod.6), main="Residuals vs. \nPredictor 4", xlab= "Time window", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$ITD)),]$wind_speed.ord.fac, resid(mod.6), main="Residuals vs. \nPredictor 5", xlab= "wind_speed", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$ITD)),]$cloudiness.z, resid(mod.6), main="Residuals vs. \nPredictor 6", xlab= "cloudiness", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$ITD)),]$sex, resid(mod.6), main="Residuals vs. \nPredictor 7", xlab= "sex", ylab= "Residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.6), log(df8[which(df8$pollinator_group == "Anthophila" & !is.na(df8$ITD)),]$ITD), main="Observed vs. fitted", xlab = "Fitted values", ylab = "Observed values", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

#mod.7
#tiff(filename = "data/PollinatorExperiment2025/Figures/20181230_mod.7_TraitFilter_resid.tif",width = 20, height = 20, units = "cm", res = 250, compression = c("none"))
par(mfrow=c(4,4))
plot(fitted(mod.7), resid(mod.7), main="(B) Residuals vs. fitted", xlab = "Fitted values", ylab="Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.7), main="Q-Q plot residuals", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(resid(mod.7))
qqnorm(ranef(mod.7)$Id.fac[,1], main="Q-Q plot random \neffect: Garden", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.7)$Id.fac[,1]) #
qqnorm(ranef(mod.7)$capture_periodMMDD[,1], main="Q-Q plot random \neffect: Date ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.7)$capture_periodMMDD[,1]) #
qqnorm(ranef(mod.7)$family[,1], main="Q-Q plot random \neffect: family ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.7)$family[,1]) #
qqnorm(ranef(mod.7)$genus[,1], main="Q-Q plot random \neffect: Genus ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.7)$genus[,1]) #
qqnorm(ranef(mod.7)$species[,1], main="Q-Q plot random \neffect: Species ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.7)$species[,1]) #
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus!="Apis" & !is.na(df8$ITD)),]$PlantS.1.z, resid(mod.7), main="Residuals vs. \nPredictor 1", xlab="Garden plant species S", ylab="Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus!="Apis" & !is.na(df8$ITD)),]$Urban_500.1.z, resid(mod.7), main="Residuals vs. \nPredictor 2", xlab= "Urban intensity", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus!="Apis" & !is.na(df8$ITD)),]$phytometer_plant, resid(mod.7), main="Residuals vs. \nPredictor 3", xlab= "phytometer_plant species", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus!="Apis" & !is.na(df8$ITD)),]$capture_window.ord.fac, resid(mod.7), main="Residuals vs. \nPredictor 4", xlab= "Time window", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus!="Apis" & !is.na(df8$ITD)),]$wind_speed.ord.fac, resid(mod.7), main="Residuals vs. \nPredictor 5", xlab= "wind_speed", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus!="Apis" & !is.na(df8$ITD)),]$cloudiness.z, resid(mod.7), main="Residuals vs. \nPredictor 6", xlab= "cloudiness", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group == "Anthophila" & df8$genus!="Apis" & !is.na(df8$ITD)),]$sex, resid(mod.7), main="Residuals vs. \nPredictor 7", xlab= "sex", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.7), log(df8[which(df8$pollinator_group == "Anthophila" & df8$genus!="Apis" & !is.na(df8$ITD)),]$ITD), main="Observed vs. fitted", xlab = "Fitted values", ylab = "Observed values", cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

#mod.8
#tiff(filename = "20181230_mod.8_TraitFilter_resid.tif",width = 20, height = 20, units = "cm", res = 250, compression = c("none"))
par(mfrow=c(4,4))
plot(fitted(mod.8), resid(mod.8), main="(C) Residuals vs. fitted", xlab = "Fitted values", ylab="Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.8), main="Q-Q plot residuals", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(resid(mod.8))
qqnorm(ranef(mod.8)$Id.fac[,1], main="Q-Q plot random \neffect: Garden", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.8)$Id.fac[,1]) #
qqnorm(ranef(mod.8)$capture_periodMMDD[,1], main="Q-Q plot random \neffect: Date ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.8)$capture_periodMMDD[,1])
qqnorm(ranef(mod.8)$family[,1], main="Q-Q plot random \neffect: family ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.8)$family[,1]) #
qqnorm(ranef(mod.8)$genus[,1], main="Q-Q plot random \neffect: Genus ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.8)$genus[,1]) #
qqnorm(ranef(mod.8)$species[,1], main="Q-Q plot random \neffect: Species ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.8)$species[,1]) #
plot(df8[which(df8$pollinator_group_sociality == "social_Bees" & df8$genus!="Apis" & !is.na(df8$ITD)),]$PlantS.1.z, resid(mod.8), main="Residuals vs. \nPredictor 1", xlab="Garden plant species S", ylab="Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "social_Bees" & df8$genus!="Apis" & !is.na(df8$ITD)),]$Urban_500.1.z, resid(mod.8), main="Residuals vs. \nPredictor 2", xlab= "Urban intensity", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "social_Bees" & df8$genus!="Apis" & !is.na(df8$ITD)),]$phytometer_plant, resid(mod.8), main="Residuals vs. \nPredictor 3", xlab= "phytometer_plant species", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "social_Bees" & df8$genus!="Apis" & !is.na(df8$ITD)),]$capture_window.ord.fac, resid(mod.8), main="Residuals vs. \nPredictor 4", xlab= "Time window", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "social_Bees" & df8$genus!="Apis" & !is.na(df8$ITD)),]$wind_speed.ord.fac, resid(mod.8), main="Residuals vs. \nPredictor 5", xlab= "wind_speed", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "social_Bees" & df8$genus!="Apis" & !is.na(df8$ITD)),]$cloudiness.z, resid(mod.8), main="Residuals vs. \nPredictor 6", xlab= "cloudiness", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "social_Bees" & df8$genus!="Apis" & !is.na(df8$ITD)),]$sex, resid(mod.8), main="Residuals vs. \nPredictor 7", xlab= "sex", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.8), log(df8[which(df8$pollinator_group_sociality == "social_Bees" & df8$genus!="Apis" & !is.na(df8$ITD)),]$ITD), main="Observed vs. fitted", xlab = "Fitted values", ylab = "Observed values", cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()


#mod.9
#tiff(filename = "20181230_mod.9_TraitFilter_resid.tif",width = 20, height = 20, units = "cm", res = 250, compression = c("none"))
par(mfrow=c(4,4))
plot(fitted(mod.9), resid(mod.9), main="(D) Residuals vs. fitted", xlab = "Fitted values", ylab="Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.9), main="Q-Q plot residuals", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(resid(mod.9))
qqnorm(ranef(mod.9)$Id.fac[,1], main="Q-Q plot random \neffect: Garden", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.9)$Id.fac[,1]) #
qqnorm(ranef(mod.9)$capture_periodMMDD[,1], main="Q-Q plot random \neffect: Date ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.9)$capture_periodMMDD[,1]) #
qqnorm(ranef(mod.9)$family[,1], main="Q-Q plot random \neffect: family ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.9)$family[,1]) #
qqnorm(ranef(mod.9)$genus[,1], main="Q-Q plot random \neffect: Genus ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.9)$genus[,1]) #
qqnorm(ranef(mod.9)$species[,1], main="Q-Q plot random \neffect: Species ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.9)$species[,1]) #
plot(df8[which(df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$ITD)),]$PlantS.1.z, resid(mod.9), main="Residuals vs. \nPredictor 1", xlab="Garden plant species S", ylab="Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$ITD)),]$Urban_500.1.z, resid(mod.9), main="Residuals vs. \nPredictor 2", xlab= "Urban intensity", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$ITD)),]$phytometer_plant, resid(mod.9), main="Residuals vs. \nPredictor 3", xlab= "phytometer_plant species", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$ITD)),]$capture_window.ord.fac, resid(mod.9), main="Residuals vs. \nPredictor 4", xlab= "Time window", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$ITD)),]$wind_speed.ord.fac, resid(mod.9), main="Residuals vs. \nPredictor 5", xlab= "wind_speed", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$ITD)),]$cloudiness.z, resid(mod.9), main="Residuals vs. \nPredictor 6", xlab= "cloudiness", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$ITD)),]$sex, resid(mod.9), main="(Residuals vs. \nPredictor 7", xlab= "sex", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.9), log(df8[which(df8$pollinator_group_sociality == "solitary_Bees" & !is.na(df8$ITD)),]$ITD), main="Observed vs. fitted", xlab = "Fitted values", ylab = "Observed values", cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

#mod.10
#tiff(filename = "20181230_mod.10_TraitFilter_resid.tif",width = 20, height = 20, units = "cm", res = 250, compression = c("none"))
par(mfrow=c(4,4))
plot(fitted(mod.10), resid(mod.10), main="(E) Residuals vs. fitted", xlab = "Fitted values", ylab="Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.10), main="Q-Q plot residuals", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(resid(mod.10))
qqnorm(ranef(mod.10)$Id.fac[,1], main="Q-Q plot random \neffect: Garden", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.10)$Id.fac[,1]) #
qqnorm(ranef(mod.10)$capture_periodMMDD[,1], main="Q-Q plot random \neffect: Date ", xlab="Theoretical quantiles", ylab="Sample quantiles", cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.10)$capture_periodMMDD[,1]) #
plot(df8[which(df8$family == "Syrphidae" & !is.na(df8$FWL)),]$PlantS.1.z, resid(mod.10), main="Residuals vs. \nPredictor 1", xlab="Garden plant species S", ylab="Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$family == "Syrphidae" & !is.na(df8$FWL)),]$Urban_500.1.z, resid(mod.10), main="Residuals vs. \nPredictor 2", xlab= "Urban intensity", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$family == "Syrphidae" & !is.na(df8$FWL)),]$phytometer_plant, resid(mod.10), main="Residuals vs. \nPredictor 3", xlab= "phytometer_plant species", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$family == "Syrphidae" & !is.na(df8$FWL)),]$capture_window.ord.fac, resid(mod.10), main="Residuals vs. \nPredictor 4", xlab= "Time window", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$family == "Syrphidae" & !is.na(df8$FWL)),]$wind_speed.ord.fac, resid(mod.10), main="Residuals vs. \nPredictor 5", xlab= "wind_speed", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$family == "Syrphidae" & !is.na(df8$FWL)),]$cloudiness.z, resid(mod.10), main="Residuals vs. \nPredictor 6", xlab= "cloudiness", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df8[which(df8$family == "Syrphidae" &!is.na(df8$FWL)),]$sex, resid(mod.10), main="Residuals vs. \nPredictor 7", xlab= "sex", ylab= "Residuals", cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.10), log(df8[which(df8$family == "Syrphidae" & !is.na(df8$FWL)),]$FWL), main="Observed vs. fitted", xlab = "Fitted values", ylab = "Observed values", cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

##########################################################
#Draw Baysian conclusions: 1. Make the Coefficient plots #
##########################################################

library(arm)
nsim <- 10000

#Run them separately for the plots (you need the bsim.mod-objects!)
bsim.mod.1 <- sim(mod.1, n.sim=nsim)
bsim.mod.2 <- sim(mod.2, n.sim=nsim)
bsim.mod.3 <- sim(mod.3, n.sim=nsim)
bsim.mod.4 <- sim(mod.4, n.sim=nsim)
bsim.mod.5 <- sim(mod.5, n.sim=nsim)
bsim.mod.6 <- sim(mod.6, n.sim=nsim)
bsim.mod.7 <- sim(mod.7, n.sim=nsim)
bsim.mod.8 <- sim(mod.8, n.sim=nsim)
bsim.mod.9 <- sim(mod.9, n.sim=nsim)
bsim.mod.10 <- sim(mod.10, n.sim=nsim)

#############################################
############### Fig4 Tongue length ##########
#############################################

#mod.1
coeff.mod.1 <- as.data.frame(t(as.data.frame(apply(bsim.mod.1@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.1$fixeff <- rownames(coeff.mod.1)
row.names(coeff.mod.1)<-NULL
colnames(coeff.mod.1)[1] <- "lower"
colnames(coeff.mod.1)[2] <- "mean"
colnames(coeff.mod.1)[3] <- "upper"
coeff.mod.1$Model <- "mod.1"
coeff.mod.1$phytometer_plant <- "all"
coeff.mod.1$Response <- "PRB.wgt"
coeff.mod.1$Apis <- "yes"
coeff.mod.1$Scale <- 500
coeff.mod.1$SizeClass <- "all"
coeff.mod.1$FG <- "allBees"
coeff.mod.1

#mod.2
coeff.mod.2 <- as.data.frame(t(as.data.frame(apply(bsim.mod.2@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.2$fixeff <- rownames(coeff.mod.2)
row.names(coeff.mod.2)<-NULL
colnames(coeff.mod.2)[1] <- "lower"
colnames(coeff.mod.2)[2] <- "mean"
colnames(coeff.mod.2)[3] <- "upper"
coeff.mod.2$Model <- "mod.2"
coeff.mod.2$phytometer_plant <- "all"
coeff.mod.2$Response <- "PRB.wgt"
coeff.mod.2$Apis <- "no"
coeff.mod.2$Scale <- 500
coeff.mod.2$SizeClass <- "all"
coeff.mod.2$FG <- "allBees_noApis"
coeff.mod.2

#mod.3
coeff.mod.3 <- as.data.frame(t(as.data.frame(apply(bsim.mod.3@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.3$fixeff <- rownames(coeff.mod.3)
row.names(coeff.mod.3)<-NULL
colnames(coeff.mod.3)[1] <- "lower"
colnames(coeff.mod.3)[2] <- "mean"
colnames(coeff.mod.3)[3] <- "upper"
coeff.mod.3$Model <- "mod.3"
coeff.mod.3$phytometer_plant <- "all"
coeff.mod.3$Response <- "PRB.wgt"
coeff.mod.3$Apis <- "no"
coeff.mod.3$Scale <- 500
coeff.mod.3$SizeClass <- "all"
coeff.mod.3$FG <- "Social_bees"
coeff.mod.3

#mod.4
coeff.mod.4 <- as.data.frame(t(as.data.frame(apply(bsim.mod.4@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.4$fixeff <- rownames(coeff.mod.4)
row.names(coeff.mod.4)<-NULL
colnames(coeff.mod.4)[1] <- "lower"
colnames(coeff.mod.4)[2] <- "mean"
colnames(coeff.mod.4)[3] <- "upper"
coeff.mod.4$Model <- "mod.4"
coeff.mod.4$phytometer_plant <- "all"
coeff.mod.4$Response <- "PRB.wgt"
coeff.mod.4$Apis <- "no"
coeff.mod.4$Scale <- 500
coeff.mod.4$SizeClass <- "all"
coeff.mod.4$FG <- "Solitary_bees"
coeff.mod.4

#Add empty vectors for pics
coeff.mod.4[nrow(coeff.mod.4)+2,] <- NA
coeff.mod.4[10,]$fixeff<-"phytometer_plantComfrey"
coeff.mod.4[11,]$fixeff<-"PlantS.1.z:phytometer_plantComfrey"
coeff.mod.4

#mod.5
coeff.mod.5 <- as.data.frame(t(as.data.frame(apply(bsim.mod.5@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.5$fixeff <- rownames(coeff.mod.5)
row.names(coeff.mod.5)<-NULL
colnames(coeff.mod.5)[1] <- "lower"
colnames(coeff.mod.5)[2] <- "mean"
colnames(coeff.mod.5)[3] <- "upper"
coeff.mod.5$Model <- "mod.5"
coeff.mod.5$phytometer_plant <- "all"
coeff.mod.5$Response <- "LPR"
coeff.mod.5$Apis <- "no"
coeff.mod.5$Scale <- 500
coeff.mod.5$SizeClass <- "all"
coeff.mod.5$FG <- "Hoverflies"
coeff.mod.5

#Add empty vectors for plots
coeff.mod.5[nrow(coeff.mod.5)+4,] <- NA
coeff.mod.5[8,]$fixeff<-"phytometer_plantSainfoin"
coeff.mod.5[9,]$fixeff<-"phytometer_plantComfrey"
coeff.mod.5[10,]$fixeff<-"PlantS.1.z:phytometer_plantSainfoin"
coeff.mod.5[11,]$fixeff<-"PlantS.1.z:phytometer_plantComfrey"

#############################################
################### ITD #####################
#############################################

#mod.6
coeff.mod.6 <- as.data.frame(t(as.data.frame(apply(bsim.mod.6@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.6$fixeff <- rownames(coeff.mod.6)
row.names(coeff.mod.6)<-NULL
colnames(coeff.mod.6)[1] <- "lower"
colnames(coeff.mod.6)[2] <- "mean"
colnames(coeff.mod.6)[3] <- "upper"
coeff.mod.6$Model <- "mod.6"
coeff.mod.6$phytometer_plant <- "all"
coeff.mod.6$Response <- "ITD"
coeff.mod.6$Apis <- "yes"
coeff.mod.6$Scale <- 500
coeff.mod.6$SizeClass <- "all"
coeff.mod.6$FG <- "allBees"
coeff.mod.6

dim(coeff.mod.6)

#mod.7
coeff.mod.7 <- as.data.frame(t(as.data.frame(apply(bsim.mod.7@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.7$fixeff <- rownames(coeff.mod.7)
row.names(coeff.mod.7)<-NULL
colnames(coeff.mod.7)[1] <- "lower"
colnames(coeff.mod.7)[2] <- "mean"
colnames(coeff.mod.7)[3] <- "upper"
coeff.mod.7$Model <- "mod.7"
coeff.mod.7$phytometer_plant <- "all"
coeff.mod.7$Response <- "ITD"
coeff.mod.7$Apis <- "no"
coeff.mod.7$Scale <- 500
coeff.mod.7$SizeClass <- "all"
coeff.mod.7$FG <- "allBees_noApis"
coeff.mod.7

dim(coeff.mod.7)

#mod.8
coeff.mod.8 <- as.data.frame(t(as.data.frame(apply(bsim.mod.8@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.8$fixeff <- rownames(coeff.mod.8)
row.names(coeff.mod.8)<-NULL
colnames(coeff.mod.8)[1] <- "lower"
colnames(coeff.mod.8)[2] <- "mean"
colnames(coeff.mod.8)[3] <- "upper"
coeff.mod.8$Model <- "mod.8"
coeff.mod.8$phytometer_plant <- "all"
coeff.mod.8$Response <- "ITD"
coeff.mod.8$Apis <- "no"
coeff.mod.8$Scale <- 500
coeff.mod.8$SizeClass <- "all"
coeff.mod.8$FG <- "Social_bees"
coeff.mod.8

dim(coeff.mod.8)

#mod.9
coeff.mod.9 <- as.data.frame(t(as.data.frame(apply(bsim.mod.9@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.9$fixeff <- rownames(coeff.mod.9)
row.names(coeff.mod.9)<-NULL
colnames(coeff.mod.9)[1] <- "lower"
colnames(coeff.mod.9)[2] <- "mean"
colnames(coeff.mod.9)[3] <- "upper"
coeff.mod.9$Model <- "mod.9"
coeff.mod.9$phytometer_plant <- "all"
coeff.mod.9$Response <- "ITD"
coeff.mod.9$Apis <- "no"
coeff.mod.9$Scale <- 500
coeff.mod.9$SizeClass <- "all"
coeff.mod.9$FG <- "Solitary_bees"
coeff.mod.9

dim(coeff.mod.9)

#Add empty vectors for plots
coeff.mod.9[nrow(coeff.mod.9)+2,] <- NA
coeff.mod.9[22,]$fixeff<-"phytometer_plantComfrey"
coeff.mod.9[23,]$fixeff<-"PlantS.1.z:phytometer_plantComfrey"

#mod.10
coeff.mod.10 <- as.data.frame(t(as.data.frame(apply(bsim.mod.10@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975)))))
coeff.mod.10$fixeff <- rownames(coeff.mod.10)
row.names(coeff.mod.10)<-NULL
colnames(coeff.mod.10)[1] <- "lower"
colnames(coeff.mod.10)[2] <- "mean"
colnames(coeff.mod.10)[3] <- "upper"
coeff.mod.10$Model <- "mod.10"
coeff.mod.10$phytometer_plant <- "all"
coeff.mod.10$Response <- "FWL"
coeff.mod.10$Apis <- "no"
coeff.mod.10$Scale <- 500
coeff.mod.10$SizeClass <- "no"
coeff.mod.10$FG <- "Hoverflies"
coeff.mod.10

dim(coeff.mod.10)

#Add empty vectors for plots
coeff.mod.10

coeff.mod.10[nrow(coeff.mod.10)+4,] <- NA
coeff.mod.10[20,]$fixeff<-"phytometer_plantComfrey"
coeff.mod.10[21,]$fixeff<-"phytometer_plantSainfoin"
coeff.mod.10[22,]$fixeff<-"PlantS.1.z:phytometer_plantSainfoin"
coeff.mod.10[23,]$fixeff<-"PlantS.1.z:phytometer_plantComfrey"
coeff.mod.10$Model <- "mod.10"
coeff.mod.10$phytometer_plant <- "all"
coeff.mod.10$Response <- "FWL"
coeff.mod.10$Apis <- "no"
coeff.mod.10$Scale <- 500
coeff.mod.10$SizeClass <- "no"
coeff.mod.10$FG <- "Hoverflies"

##############################
### Make coefficient plots ###
##############################

fixeff_labels <- c(
  "sexmale" = "sex: male",
  "cloudiness.z" = "% cloudiness",
  "wind_speed.ord.fac.C" = "WS: 3 Bft.",
  "wind_speed.ord.fac.Q" = "WS: 2 Bft.",
  "wind_speed.ord.fac.L" = "WS: 1 Bft.",
  "capture_window.ord.fac^8" = "17:00-18:00",
  "capture_window.ord.fac^7" = "16:00-17:00",
  "capture_window.ord.fac^6" = "15:00-16:00",
  "capture_window.ord.fac^5" = "14:00-15:00",
  "capture_window.ord.fac^4" = "13:00-14:00",
  "capture_window.ord.fac.C" = "12:00-13:00",
  "capture_window.ord.fac.Q" = "11:00-12:00",
  "capture_window.ord.fac.L" = "10:00-11:00",
  "PlantS.1.z:phytometer_plantComfrey" = "Local × comfrey",
  "PlantS.1.z:phytometer_plantSainfoin" = "Local × sainfoin",
  "PlantS.1.z:phytometer_plantRadish" = "Local × radish",
  "phytometer_plantComfrey" = "Comfrey",
  "phytometer_plantSainfoin" = "Sainfoin",
  "phytometer_plantRadish" = "Radish",
  "PlantS.1.z:Urban_500.1.z" = "Local × landscape",
  "Urban_500.1.z" = "Landscape-scale\nhabitat loss",
  "PlantS.1.z" = "Local plant\nspecies richness",
  "All Bees" = "with honeybees",
  "All Bees (no honeybees)" = "without honeybees"
)

########################
## 1. Tongue lengths! ##
########################

coeff.mod.1.ord <- coeff.mod.1[c(7,9,11,10,4,6,5,8,3,2),]
coeff.mod.2.ord <- coeff.mod.2[c(7,9,11,10,4,6,5,8,3,2),]
coeff.mod.3.ord <- coeff.mod.3[c(7,9,11,10,4,6,5,8,3,2),]
coeff.mod.4.ord <- coeff.mod.4[c(6,11,9,8,10,5,4,7,3,2),]
coeff.mod.5.ord <- coeff.mod.5[c(5,11,10,7,9,8,4,6,3,2),]

coeff.mod.1.ord
coeff.mod.2.ord
coeff.mod.3.ord
coeff.mod.4.ord
coeff.mod.5.ord


#plot all bees
library(dplyr)
library(ggplot2)

# Combine data into a single dataframe
plot_data <- bind_rows(
  coeff.mod.2.ord %>% mutate(group = "All Bees"),
  coeff.mod.1.ord %>% mutate(group = "All Bees (no honeybees)")
)

# Ensure fixeff is a factor for correct ordering
plot_data <- plot_data %>%
  mutate(fixeff = factor(fixeff, levels = rev((unique(fixeff)))))  # Reverse order for readability

fig2a <- 
  ggplot(plot_data, aes(x = mean, y = fixeff, shape = group)) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "#E01A4F", linewidth=.8) +  # Black vertical reference line
    #geom_point(size = 4, color = "black", stroke = 0.8) +
    geom_point(aes(color = (lower * upper > 0)), size = 4, stroke = 0.8) +  # Keep all points black
    # Define colors manually: black for significant, grey for uncertain
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
    coord_flip()+
    # geom_errorbarh(aes(xmin = lower, xmax = upper), 
    #                height = 0, linewidth = 0.8, color = "black") +  # Black error bars
    # # Error bars: grey if non-significant, black otherwise
    geom_errorbarh(aes(xmin = lower, xmax = upper, color = (lower * upper > 0)), 
                 height = 0, linewidth = .8) +  # Thicker error bars
    scale_x_continuous(limits = c(min(plot_data$lower), max(plot_data$upper))) +
    scale_shape_manual(values = c("All Bees" = 16, "All Bees (no honeybees)" = 4)) +  # Different shapes
    labs(x = "", y = "", shape = "Group",
         title = "(a) all bees") +  # Shape legend
    theme_classic(base_size = 14) +
    scale_y_discrete(labels = NULL) +
    theme(
      
      axis.text.y = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.title.x = element_text(hjust = 0, color = "black"),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(color = "black"),
      legend.position = "none",  # Move legend inside (x, y) in relative coordinates
      legend.justification = c(1,0),  # Align bottom-right corner of legend at this position
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),  
      legend.key = element_blank()  # Optional: removes legend key box background
    )
fig2a
#Plot social bees
coeff.mod.3.ord <- coeff.mod.3.ord %>%
  mutate(fixeff = factor(fixeff, levels = rev(unique(fixeff)))) 

fig2b <- ggplot(coeff.mod.3.ord, aes(x = mean, y = fixeff)) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "#E01A4F", linewidth=.8) +  # Black vertical reference line
  geom_point(aes(color = (lower * upper > 0)), size = 4, stroke = 0.8) +  # Keep all points black
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  coord_flip()+
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = (lower * upper > 0)), 
                 height = 0, linewidth = .8) +  # Thicker error bars
  scale_x_continuous(limits = c(min(coeff.mod.3.ord$lower), max(coeff.mod.3.ord$upper))) +
  labs(x = "", 
       y = "",
       title = "(b) wild social bees") +
  scale_y_discrete(labels = NULL) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",  # Keep legend for shape
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank(),
  )
fig2b

#Plot solitary bees

coeff.mod.4.ord <- coeff.mod.4.ord %>%
  mutate(fixeff = factor(fixeff, levels = rev(unique(fixeff))))  # Reverse order for readability

fig2c <- 
  ggplot(coeff.mod.4.ord, aes(x = mean, y = fixeff)) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "#E01A4F", linewidth=.8) +  # Black vertical reference line
  geom_point(aes(color = (lower * upper > 0)), size = 4, stroke = 0.8) +  # Keep all points black
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = (lower * upper > 0)), 
                 height = 0, linewidth = .8) +  # Thicker error bars
  coord_flip()+
    scale_x_continuous(limits = c(min(coeff.mod.4.ord$lower, na.rm = T), 
                                  max(coeff.mod.4.ord$upper, na.rm = T))) +
    labs(x = "", 
         y = "",
         title = "(c) solitary bees") +
    theme_classic(base_size = 14) +
    scale_y_discrete(labels = NULL) +
    theme(
      legend.position = "none",  # Keep legend for shape
      axis.text.y = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.title.x = element_text(hjust = 0, color = "black"),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.y = element_blank()
    )

fig2c

#Plot hoverflies

coeff.mod.5.ord <- coeff.mod.5.ord %>%
  mutate(fixeff = factor(fixeff, levels = rev(unique(fixeff))))  # Reverse order for readability

fig2d <- 
  ggplot(coeff.mod.5.ord, aes(x = mean, y = fixeff)) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "#E01A4F", linewidth=.8) +  # Black vertical reference line
    geom_point(aes(color = (lower * upper > 0)), size = 4, stroke = 0.8) +  # Keep all points black
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
    geom_errorbarh(aes(xmin = lower, xmax = upper, color = (lower * upper > 0)), 
                 height = 0, linewidth = .8) +  # Thicker error bars
  coord_flip()+
    scale_x_continuous(limits = c(min(coeff.mod.5.ord$lower, na.rm = T), 
                                  max(coeff.mod.5.ord$upper, na.rm = T))) +
    labs(x = "", 
         y = "",
         title = "(d) hoverflies") +
    scale_y_discrete(labels = fixeff_labels) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",  # Keep legend for shape
      axis.text.y = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = .5,hjust = 1),
      axis.title.x = element_text(hjust = 0, color = "black"),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.y= element_blank()
    )

fig2d

library(patchwork)
fig2a /fig2b/fig2c /fig2d
# save pdf 8 x 14

####################
## 2. Body sizes! ##
####################

coeff.mod.6.ord <- coeff.mod.6[c(19,18,17,16,15,14,13,12,11,10,9,8,7,21,23,22,4,6,5,20,3,2),]
coeff.mod.7.ord <- coeff.mod.7[c(19,18,17,16,15,14,13,12,11,10,9,8,7,21,23,22,4,6,5,20,3,2),]
coeff.mod.8.ord <- coeff.mod.8[c(19,18,17,16,15,14,13,12,11,10,9,8,7,21,23,22,4,6,5,20,3,2),]
coeff.mod.9.ord <- coeff.mod.9[c(18,17,16,15,14,13,12,11,10,9,8,7,6,23,21,20,22,5,4,19,3,2),]
coeff.mod.10.ord <- coeff.mod.10[c(17,16,15,14,13,12,11,10,9,8,7,6,5,23,22,19,20,21,4, 18,3,2),]

dim(coeff.mod.6.ord)
dim(coeff.mod.7.ord)
dim(coeff.mod.8.ord)
dim(coeff.mod.9.ord)
dim(coeff.mod.10.ord)

coeff.mod.6.ord

#Define the position of the factors on the y axis
c <- cbind(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43)
d <- c-0.75
d

#Plot all bees (+/- honeybees)

# Combine datasets into a single dataframe
plot_data <- bind_rows(
  coeff.mod.7.ord %>% mutate(group = "All Bees"),
  coeff.mod.6.ord %>% mutate(group = "All Bees (no honeybees)")
)

# Ensure fixeff is a factor for correct ordering
plot_data <- plot_data %>%
  mutate(fixeff = factor(fixeff, levels = rev(unique(fixeff))))  # Reverse order for readability

fig2e <- 
  ggplot(plot_data, aes(x = mean, y = fixeff, shape = group)) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "#E01A4F", linewidth = 0.8) +  # Reference line
    geom_point(aes(color = (lower * upper > 0)), size = 4, stroke = 0.8) +  # Keep all points black
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgray")) +
  coord_flip()+
    geom_errorbarh(aes(xmin = lower, xmax = upper, color = (lower * upper > 0)), 
                   height = 0, linewidth = .8) +  # Thicker error bars
    scale_x_continuous(limits = c(min(plot_data$lower), max(plot_data$upper))) +
    scale_shape_manual(values = c("All Bees" = 16, "All Bees (no honeybees)" = 4)) +  # Different shapes
    
    labs(x = "", 
         y = "", title = "(a) all bees")  +  # Shape legend
    
    theme_classic(base_size = 14) +
    scale_y_discrete(labels = NULL) +
    theme(
      legend.position = "none",  # Keep legend for shape
      axis.text.y = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.title.x = element_text(hjust = 0, color = "black"),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(color = "black")
    )

fig2e 

# social bees
coeff.mod.8.ord <- coeff.mod.8.ord %>%
  mutate(fixeff = factor(fixeff, levels = rev(unique(fixeff))))  # Reverse order for readability

fig2f <- 
  ggplot(coeff.mod.8.ord, aes(x = mean, y = fixeff)) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "#E01A4F", linewidth=.8) +  # Vertical reference line
  
  geom_point(aes(color = (lower * upper > 0)), size = 4, stroke = 0.8) +  # Keep all points black
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  coord_flip()+
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = (lower * upper > 0)), 
                 height = 0, linewidth = .8) +  # Thicker error bars
  
  scale_x_continuous(limits = c(min(coeff.mod.8.ord$lower), max(coeff.mod.8.ord$upper))) +
  
  labs(x = "", y = "",
       title = "(b) wild social bees") +
  theme_classic(base_size = 14) +
  scale_y_discrete(labels = NULL) +
  theme(
    legend.position = "none",  # Remove legend
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank()
  )

fig2e / fig2f
#Plot solitary bees

# social bees
coeff.mod.9.ord <- coeff.mod.9.ord %>%
  mutate(fixeff = factor(fixeff, levels = rev(unique(fixeff))))  # Reverse order for readability

fig2g <- ggplot(coeff.mod.9.ord, aes(x = mean, y = fixeff)) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "#E01A4F", linewidth=.8) +  # Black vertical reference line
  geom_point(aes(color = (lower * upper > 0)), size = 4, stroke = 0.8) +  # Keep all points black
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  coord_flip()+
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = (lower * upper > 0)), 
                 height = 0, linewidth = .8) +  # Thicker error bars
  scale_x_continuous(limits = c(min(coeff.mod.9.ord$lower), max(coeff.mod.9.ord$upper))) +
  labs(x = "", y = "", title = "(c) solitary bees") +
  theme_classic(base_size = 14) +
  scale_y_discrete(labels = NULL) +
  theme(
    legend.position = "none",  # Keep legend for shape
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank()
  )

fig2g

#Plot hoverflies
coeff.mod.10.ord <- coeff.mod.10.ord %>%
  mutate(fixeff = factor(fixeff, levels = rev(unique(fixeff))))  # Reverse order for readability

fig2h <- ggplot(coeff.mod.10.ord, aes(x = mean, y = fixeff)) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "#E01A4F", linewidth=.8) +  # Black vertical reference line
  geom_point(aes(color = (lower * upper > 0)), size = 4, stroke = 0.8) +  # Keep all points black
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  coord_flip()+
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = (lower * upper > 0)), 
                 height = 0, linewidth = .8) +  # Thicker error bars
  scale_x_continuous(limits = c(min(coeff.mod.10.ord$lower), max(coeff.mod.10.ord$upper))) +
  labs(x = "", y = "", title = "(d) hoverflies") +
  theme_classic(base_size = 14) +
  scale_y_discrete(labels = fixeff_labels) +
  theme(
    legend.position = "none",  # Keep legend for shape
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = .5,hjust = 1),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank()
  )

fig2e / fig2f / fig2g / fig2h

# Save the figure as a single PDF

###############################################
#### Merge all coefficients ###################
###############################################

coeff.mods.TraitFiltering <- rbind(coeff.mod.1,
                                   coeff.mod.2,
                                   coeff.mod.3,
                                   coeff.mod.4,
                                   coeff.mod.5,
                                   coeff.mod.6,
                                   coeff.mod.7,
                                   coeff.mod.8,
                                   coeff.mod.9,
                                   coeff.mod.10)
head(coeff.mods.TraitFiltering)
tail(coeff.mods.TraitFiltering)
str(coeff.mods.TraitFiltering)

#Recode the variables:
table(coeff.mods.TraitFiltering$fixeff)

coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="(Intercept)"]<-"Intercept"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="PlantS.1.z"]<-"GardenPlant_S"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="Urban_500.1.z"]<-"UrbanIntensity_500"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="capture_window.ord.fac.L"]<-"TimeWindow:10-11"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="capture_window.ord.fac.Q"]<-"TimeWindow:11-12"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="capture_window.ord.fac.C"]<-"TimeWindow:12-13"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="capture_window.ord.fac^4"]<-"TimeWindow:13-14"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="capture_window.ord.fac^5"]<-"TimeWindow:14-15"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="capture_window.ord.fac^6"]<-"TimeWindow:15-16"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="capture_window.ord.fac^7"]<-"TimeWindow:16-17"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="capture_window.ord.fac^8"]<-"TimeWindow:17-18"
#coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="capture_window.ord.fac^9"]<-"TimeWindow:18-19"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="wind_speed.ord.fac.L"]<-"wind_speed:1"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="wind_speed.ord.fac.Q"]<-"wind_speed:2"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="wind_speed.ord.fac.C"]<-"wind_speed:3"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="cloudiness.z"]<-"cloudiness"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="sexmale"]<-"sex:Male"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="PlantS.1.z:Urban_500.1.z"]<-"GardenPlantxUrbanIntensity"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="phytometer_plantComfrey"]<-"phytometer_plant:Comfrey"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="phytometer_plantRadish"]<-"phytometer_plant:Radish"
coeff.mods.TraitFiltering$fixeff[coeff.mods.TraitFiltering$fixeff=="phytometer_plantSainfoin"]<-"phytometer_plant:Sainfoin"

coeff.mods.TraitFiltering$fixeff <- as.factor(coeff.mods.TraitFiltering$fixeff)

#Get rid of new rows for the plots
dim(coeff.mods.TraitFiltering)
coeff.mods.TraitFiltering<-na.omit(coeff.mods.TraitFiltering)
dim(coeff.mods.TraitFiltering)

#Round them and take only the relevant columns
coeff.mods.TraitFiltering

coeff.mods.TraitFiltering.table <- cbind(coeff.mods.TraitFiltering[,c(5,4)], as.data.frame(apply(coeff.mods.TraitFiltering[,1:3], 2, function(x)round(x,3))))
coeff.mods.TraitFiltering.table

coeff.mods.TraitFiltering.table <- coeff.mods.TraitFiltering.table[,c(1,2,4,3,5)]
coeff.mods.TraitFiltering.table

dim(coeff.mods.TraitFiltering.table)
str(coeff.mods.TraitFiltering.table)

rownames(coeff.mods.TraitFiltering.table) <- c(1:158) 

coeff.mods.TraitFiltering.table_tongue <- coeff.mods.TraitFiltering.table[1:49,]
coeff.mods.TraitFiltering.table_ITD <- coeff.mods.TraitFiltering.table[50:158,]

#Save them

#Tongue
write.table(coeff.mods.TraitFiltering.table_tongue, file = "results/Beta_coefficients_TraitFiltering_tongue.txt", 
            append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

#ITD
write.table(coeff.mods.TraitFiltering.table_ITD, file = "results/Beta_coefficients_TraitFiltering_ITD.txt", 
            append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

###################################################
#Make a coefficents table out of it and merge them#
###################################################
table(coeff.mods.TraitFiltering$fixeff)

#Divide first the table accoridng to their predictor variable:
dat.x0 <- subset(coeff.mods.TraitFiltering, fixeff == "Intercept")
dat.x1 <- subset(coeff.mods.TraitFiltering, fixeff == "GardenPlant_S")
dat.x2 <- subset(coeff.mods.TraitFiltering, fixeff == "UrbanIntensity_500")
dat.x3 <- subset(coeff.mods.TraitFiltering, fixeff == "GardenPlantxUrbanIntensity")
dat.x4 <- subset(coeff.mods.TraitFiltering, fixeff == "Phytometer:Radish")
dat.x5 <- subset(coeff.mods.TraitFiltering, fixeff == "Phytometer:Sainfoin")
dat.x6 <- subset(coeff.mods.TraitFiltering, fixeff == "Phytometer:Comfrey")
dat.x7 <- subset(coeff.mods.TraitFiltering, fixeff == "TimeWindow:10-11")
dat.x8 <- subset(coeff.mods.TraitFiltering, fixeff == "TimeWindow:11-12")
dat.x9 <- subset(coeff.mods.TraitFiltering, fixeff == "TimeWindow:12-13")
dat.x10 <- subset(coeff.mods.TraitFiltering, fixeff == "TimeWindow:13-14")
dat.x11 <- subset(coeff.mods.TraitFiltering, fixeff == "TimeWindow:14-15")
dat.x12 <- subset(coeff.mods.TraitFiltering, fixeff == "TimeWindow:15-16")
dat.x13 <- subset(coeff.mods.TraitFiltering, fixeff == "TimeWindow:16-17")
dat.x14 <- subset(coeff.mods.TraitFiltering, fixeff == "TimeWindow:17-18")
dat.x15 <- subset(coeff.mods.TraitFiltering, fixeff == "TimeWindow:18-19")
#dat.x16 <- subset(coeff.mods.TraitFiltering, fixeff == "TimeWindow:18-19")
dat.x16 <- subset(coeff.mods.TraitFiltering, fixeff == "Windspeed:1")
dat.x17 <- subset(coeff.mods.TraitFiltering, fixeff == "Windspeed:2")
dat.x18 <- subset(coeff.mods.TraitFiltering, fixeff == "Windspeed:3")
#dat.x19 <- subset(coeff.mods.TraitFiltering, fixeff == "Windspeed")
dat.x19 <- subset(coeff.mods.TraitFiltering, fixeff == "Cloudiness")
dat.x20 <- subset(coeff.mods.TraitFiltering, fixeff == "Sex:Male")

#Rename the Mean-CI columns:
colnames(dat.x0)[c(1,2,3)]<-c("2.5%_CI_x0","Mean_x0","97.5%_CI_x0")
colnames(dat.x1)[c(1,2,3)]<-c("2.5%_CI_x1","Mean_x1","97.5%_CI_x1")
colnames(dat.x2)[c(1,2,3)]<-c("2.5%_CI_x2","Mean_x2","97.5%_CI_x2")
colnames(dat.x3)[c(1,2,3)]<-c("2.5%_CI_x3","Mean_x3","97.5%_CI_x3")
colnames(dat.x4)[c(1,2,3)]<-c("2.5%_CI_x4","Mean_x4","97.5%_CI_x4")
colnames(dat.x5)[c(1,2,3)]<-c("2.5%_CI_x5","Mean_x5","97.5%_CI_x5")
colnames(dat.x6)[c(1,2,3)]<-c("2.5%_CI_x6","Mean_x6","97.5%_CI_x6")
colnames(dat.x7)[c(1,2,3)]<-c("2.5%_CI_x7","Mean_x7","97.5%_CI_x7")
colnames(dat.x8)[c(1,2,3)]<-c("2.5%_CI_x8","Mean_x8","97.5%_CI_x8")
colnames(dat.x9)[c(1,2,3)]<-c("2.5%_CI_x9","Mean_x9","97.5%_CI_x9")
colnames(dat.x10)[c(1,2,3)]<-c("2.5%_CI_x10","Mean_x10","97.5%_CI_x10")
colnames(dat.x11)[c(1,2,3)]<-c("2.5%_CI_x11","Mean_x11","97.5%_CI_x11")
colnames(dat.x12)[c(1,2,3)]<-c("2.5%_CI_x12","Mean_x12","97.5%_CI_x12")
colnames(dat.x13)[c(1,2,3)]<-c("2.5%_CI_x13","Mean_x13","97.5%_CI_x13")
colnames(dat.x14)[c(1,2,3)]<-c("2.5%_CI_x14","Mean_x14","97.5%_CI_x14")
colnames(dat.x15)[c(1,2,3)]<-c("2.5%_CI_x15","Mean_x15","97.5%_CI_x15")
colnames(dat.x16)[c(1,2,3)]<-c("2.5%_CI_x16","Mean_x16","97.5%_CI_x16")
colnames(dat.x17)[c(1,2,3)]<-c("2.5%_CI_x17","Mean_x17","97.5%_CI_x17")
colnames(dat.x18)[c(1,2,3)]<-c("2.5%_CI_x18","Mean_x18","97.5%_CI_x18")
colnames(dat.x19)[c(1,2,3)]<-c("2.5%_CI_x19","Mean_x19","97.5%_CI_x19")
colnames(dat.x20)[c(1,2,3)]<-c("2.5%_CI_x20","Mean_x20","97.5%_CI_x20")

dat.comb <- Reduce(function(x, y) merge(x, y, by="Model", all = T, sort=T), list(dat.x0[,c(1,2,3,5)], 
                                                                                 dat.x1[,c(1,2,3,5)], 
                                                                                 dat.x2[,c(1,2,3,5)],
                                                                                 dat.x3[,c(1,2,3,5)],
                                                                                 dat.x4[,c(1,2,3,5)],
                                                                                 dat.x5[,c(1,2,3,5)],
                                                                                 dat.x6[,c(1,2,3,5)],
                                                                                 dat.x7[,c(1,2,3,5)],
                                                                                 dat.x8[,c(1,2,3,5)],
                                                                                 dat.x9[,c(1,2,3,5)],
                                                                                 dat.x10[,c(1,2,3,5)],
                                                                                 dat.x11[,c(1,2,3,5)],
                                                                                 dat.x12[,c(1,2,3,5)],
                                                                                 dat.x13[,c(1,2,3,5)],
                                                                                 dat.x14[,c(1,2,3,5)],
                                                                                 dat.x15[,c(1,2,3,5)],
                                                                                 dat.x16[,c(1,2,3,5)],
                                                                                 dat.x17[,c(1,2,3,5)],
                                                                                 dat.x18[,c(1,2,3,5)],
                                                                                 dat.x19[,c(1,2,3,5)],
                                                                                 dat.x20[,c(1,2,3,5)]
))

dim(dat.comb)
head(dat.comb)
tail(dat.comb)
dat.comb

dat.comb <- dat.comb[c(1,3:10,2),c("Model",
                                   "Mean_x0","2.5%_CI_x0","97.5%_CI_x0",
                                   "Mean_x1","2.5%_CI_x1","97.5%_CI_x1",
                                   "Mean_x2","2.5%_CI_x2","97.5%_CI_x2",
                                   "Mean_x3","2.5%_CI_x3","97.5%_CI_x3",
                                   "Mean_x4","2.5%_CI_x4","97.5%_CI_x4",
                                   "Mean_x5","2.5%_CI_x5","97.5%_CI_x5",
                                   "Mean_x6","2.5%_CI_x6","97.5%_CI_x6",
                                   "Mean_x7","2.5%_CI_x7","97.5%_CI_x7",
                                   "Mean_x8","2.5%_CI_x8","97.5%_CI_x8",
                                   "Mean_x9","2.5%_CI_x9","97.5%_CI_x9",
                                   "Mean_x10","2.5%_CI_x10","97.5%_CI_x10",
                                   "Mean_x11","2.5%_CI_x11","97.5%_CI_x11",
                                   "Mean_x12","2.5%_CI_x12","97.5%_CI_x12",
                                   "Mean_x13","2.5%_CI_x13","97.5%_CI_x13",
                                   "Mean_x14","2.5%_CI_x14","97.5%_CI_x14",
                                   "Mean_x15","2.5%_CI_x15","97.5%_CI_x15",
                                   "Mean_x16","2.5%_CI_x16","97.5%_CI_x16",
                                   "Mean_x17","2.5%_CI_x17","97.5%_CI_x17",
                                   "Mean_x18","2.5%_CI_x18","97.5%_CI_x18",
                                   "Mean_x19","2.5%_CI_x19","97.5%_CI_x19",
                                   "Mean_x20","2.5%_CI_x20","97.5%_CI_x20")]
dat.comb

#Round the values
dim(dat.comb)

dat.comb <- cbind(dat.comb[,1], as.data.frame(apply(dat.comb[,2:64], 2, function(x)round(x,2))))
colnames(dat.comb)[1]<-"Model"

dat.comb

#Save the values
write.table(dat.comb, file = "results/Beta_coefficients_TraitFiltering_all_Table4Paper.txt", 
            append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
