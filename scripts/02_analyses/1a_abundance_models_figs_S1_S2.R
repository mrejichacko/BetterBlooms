#####################################################################
#####################################################################
#####################################################################
###
### 1a. Modelling abundance
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
### 
#####################################################################
#####################################################################
#####################################################################
###
### Prep

rm(list=ls())

floristic<-read.table("raw_data/plant_floristic_data.txt", header = TRUE, sep = ";")
explanatory<-read.table("raw_data/explanatory_variables.txt", header = TRUE, sep = ";")
abundance<-read.table("cleaned_data/pollinator_abundance_aggregated_garden_phytometer.txt", header = TRUE, sep = ";")
richness<-read.table("cleaned_data/pollinator_SR_aggregated_garden_phytometer.txt", header = TRUE, sep = ";")
effort<-read.table("cleaned_data/sampling_effort_gardens_plants_aggregated.txt", header = TRUE, sep = ";")

str(floristic)
str(explanatory)
str(abundance)
str(richness)
str(effort)

dim(explanatory)

#Combine the datasets
df <- Reduce(function(x, y) merge(x, y, by="Id", sort=TRUE), 
             list(floristic[,c(1,5)],
                  explanatory[,c(1,4,5,13:16)],
                  abundance,
                  richness,
                  effort))

summary(df)
dim(df)
str(df)
names(df)

#Rename the plant variable
colnames(df)[2]<-"PlantS"
colnames(df)[2]

#Kick out garden Nr. 39, the sampling effort is too low
df <- droplevels(df[which(df$Id != 39),])
dim(df)

rm(abundance, effort, explanatory, floristic, richness)

###########################
### Scale the variables ###
###########################

df$Id.fac <- as.factor(df$Id)

##Landscape, all scales:

poly.Urban_500 <- poly(df$Urban_500,2)
df$Urban_500.1.z <- scale(poly.Urban_500[,1])
df$Urban_500.2.z <- scale(poly.Urban_500[,2])

poly.Urban_250 <- poly(df$Urban_250,2)
df$Urban_250.1.z <- scale(poly.Urban_250[,1])
df$Urban_250.2.z <- scale(poly.Urban_250[,2])

poly.Urban_100 <- poly(df$Urban_100,2)
df$Urban_100.1.z <- scale(poly.Urban_100[,1])
df$Urban_100.2.z <- scale(poly.Urban_100[,2])

poly.Urban_50 <- poly(df$Urban_50,2)
df$Urban_50.1.z <- scale(poly.Urban_50[,1])
df$Urban_50.2.z <- scale(poly.Urban_50[,2])

##Garden:

#Plants species richness May-August
poly.PlantS <- poly(df$PlantS, 2)
df$PlantS.1.z <- scale(poly.PlantS[,1])
df$PlantS.2.z <- scale(poly.PlantS[,2])

##Sampling:

#Check the sampling variables and transform them in daily rates:
df$Total_sampling_effort_d <- df$total_sampling_effort_min/60/9 #first minutes to h, then  per fieldwork day = 9h (one fieldwork day lasted from 9:00 to 18:00)
hist(df$Total_sampling_effort_d) 

df$n_flowers_carrot.z<-scale(df$n_flowers_carrot)
df$n_flowers_radish.z<-scale(df$n_flowers_radish)
df$n_flowers_sainfoin.z<-scale(df$n_flowers_sainfoin)
df$n_flowers_comfrey.z<-scale(df$n_flowers_comfrey)

df$n_flowers.z <-scale((df$n_flowers_carrot+df$n_flowers_radish+df$n_flowers_sainfoin+df$n_flowers_comfrey))

###################################
### Happy modelling: Abundance~ ###
###################################

library(lme4)
library(blmeco)
library(plyr)
library(AICcmodavg)
library(MuMIn)
?r.squaredGLMM

####################################################
## First step: Make a list of the selected models ##
####################################################

modlist.1 <- list(mod.1=glmer(A_allBees~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                  mod.2=glmer(A_allBees_noApis ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                  mod.3=glmer(A_Apis ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                  mod.4=glmer(A_Bombus ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                  mod.5=glmer(A_socialBees ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                  mod.6=glmer(A_solitaryBees ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                  mod.7=glmer(A_otherAculeata ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                  mod.8=glmer(A_Syrphidae ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                  mod.9=glmer(A_Coleoptera ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                  mod.10=glmer(A_allPollinators ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)))

#run them also outside the list:
mod.1<-glmer(A_allBees ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.2<-glmer(A_allBees_noApis ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.3<-glmer(A_Apis ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +(1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.4<-glmer(A_Bombus ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.5<-glmer(A_socialBees ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +(1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.6<-glmer(A_solitaryBees ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.7<-glmer(A_otherAculeata ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +   (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.8<-glmer(A_Syrphidae ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.9<-glmer(A_Coleoptera ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.10<-glmer(A_allPollinators ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))

summary(mod.1)

#Build table with model names, function, AIC and AICc
#Increase WAIC simulations to 10000! 
modlist.1_Fit <- ldply(modlist.1, 
                       function(mod) {
                         data.frame(#WAIC = round(WAIC(mod,nsim=1000)$WAIC2,2), 
                           AICc = round(AICc(mod),2),
                           #R2m = round(r.squaredGLMM(mod)[3],2), # function from piecewiseSEM
                           #R2c = round(r.squaredGLMM(mod)[6],2), # function from piecewiseSEM
                           DispPar = dispersion_glmer(mod), 
                           meanRF = round(mean(ranef(mod)$Id.fac[,1]),5), 
                           shapiroRANEF = round(shapiro.test(ranef(mod)$Id.fac[,1])$p.value,3), 
                           shapiroRESID = round(shapiro.test(resid(mod))$p.value,3))}) #model = deparse(formula(mod)) -> to plot the model formula

colnames(modlist.1_Fit)[1]<-"Model"
modlist.1_Fit 

###################################################################
## Second step: check model assumptions with the selected models ##
###################################################################

clab = 0.8
cmain = 0.7
caxis = 0.8

par(mfrow=c(10,6),
    oma = c(1,4.8,4.8,0) + 0.1,
    mar = c(1,0,1,2.9) + 0.1,
    cex.lab = 2.5)

#Model 1 
#par(mfrow=c(1,6))
plot(fitted(mod.1), resid(mod.1), main="Residuals vs. fitted", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.1), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main="Q-Q plot residuals") #
qqline(resid(mod.1))
qqnorm(ranef(mod.1)$Id.fac[,1], main="Q-Q plot RF: Garden", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.1)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.1), main="Residuals vs. Predictor 1", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.1), main="Residuals vs. Predictor 2", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.1), df$A_allBees, main="Observed vs. fitted", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

#Model 2 
#par(mfrow=c(1,6))
plot(fitted(mod.2), resid(mod.2), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.2), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.2))
qqnorm(ranef(mod.2)$Id.fac[,1], cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(ranef(mod.2)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.2), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.2), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.2), df$A_allBees_noApis, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

#Model 3 
#par(mfrow=c(1,6))
plot(fitted(mod.3), resid(mod.3), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.3), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.3))
qqnorm(ranef(mod.3)$Id.fac[,1], cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(ranef(mod.3)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.3), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.3), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.3), df$A_Apis, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

#Model 4 
#par(mfrow=c(1,6))
plot(fitted(mod.4), resid(mod.4), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.4), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.4))
qqnorm(ranef(mod.4)$Id.fac[,1], cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(ranef(mod.4)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.4), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.4), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.4), df$A_Bombus, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

#Model 5 
#par(mfrow=c(1,6))
plot(fitted(mod.5), resid(mod.5), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.5),  cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.5))
qqnorm(ranef(mod.5)$Id.fac[,1],cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(ranef(mod.5)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.5), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.5), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.5), df$A_socialBees, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

#Model 6 
#par(mfrow=c(1,6))
plot(fitted(mod.6), resid(mod.6), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.6), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.6))
qqnorm(ranef(mod.6)$Id.fac[,1], cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(ranef(mod.6)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.6), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.6), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.6), df$A_solitaryBees, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

#Model 7 
plot(fitted(mod.7), resid(mod.7), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.7), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.7))
qqnorm(ranef(mod.7)$Id.fac[,1],cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(ranef(mod.7)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.7), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.7),cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.7), df$A_otherAculeata, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

#Model 8 
plot(fitted(mod.8), resid(mod.8), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.8), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.8))
qqnorm(ranef(mod.8)$Id.fac[,1],  cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(ranef(mod.8)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.8),  cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.8), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.8), df$A_Syrphidae, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

#Model 9 
plot(fitted(mod.9), resid(mod.9), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.9), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.9))
qqnorm(ranef(mod.9)$Id.fac[,1], cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(ranef(mod.9)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.9), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.9), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.9), df$A_Coleoptera, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

#Model 10 
plot(fitted(mod.10), resid(mod.10),xlab = "Fitted values", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.10), cex.lab=clab,  xlab="Theoretical quantiles", cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.10))
qqnorm(ranef(mod.10)$Id.fac[,1], cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(ranef(mod.10)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.10), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.10), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.10), df$A_allPollinators, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

#dev.off()

##########################################################################
#Spatial correlation: make one big multipanel plot including all figures:#
##########################################################################

library(sp)

#Model 1:
spdata_mod.1 <- data.frame(resid=resid(mod.1), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata_mod.1) <- c("x","y")
bubble(spdata_mod.1, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 2:
spdata_mod.2 <- data.frame(resid=resid(mod.2), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata_mod.2) <- c("x","y")
bubble(spdata_mod.2, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 3:
spdata.mod.3 <- data.frame(resid=resid(mod.3), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata.mod.3) <- c("x","y")
bubble(spdata.mod.3, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 4:
spdata.mod.4 <- data.frame(resid=resid(mod.4), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata.mod.4) <- c("x","y")
bubble(spdata.mod.4, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 5:
spdata.mod.5 <- data.frame(resid=resid(mod.5), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata.mod.5) <- c("x","y")
bubble(spdata.mod.5, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 6:
spdata.mod.6 <- data.frame(resid=resid(mod.6), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata.mod.6) <- c("x","y")
bubble(spdata.mod.6, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 7:
spdata.mod.7 <- data.frame(resid=resid(mod.7), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata.mod.7) <- c("x","y")
bubble(spdata.mod.7, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 8: !!
spdata.mod.8 <- data.frame(resid=resid(mod.8), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata.mod.8) <- c("x","y")
bubble(spdata.mod.8, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 9:
spdata.mod.9 <- data.frame(resid=resid(mod.9), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata.mod.9) <- c("x","y")
bubble(spdata.mod.9, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 10:
spdata.mod.10 <- data.frame(resid=resid(mod.10), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata.mod.10) <- c("x","y")
bubble(spdata.mod.10, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")


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

#Collect them also in a list:
modlist.500_Coeffs <- ldply(modlist.1, function(mod) {coeff.mod = (as.data.frame(t(as.data.frame(apply(sim(mod, n.sim=nsim)@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975))))))})

modlist.500_Coeffs

rownames(modlist.500_Coeffs)<-NULL
colnames(modlist.500_Coeffs)[1]<-"Model"
colnames(modlist.500_Coeffs)[2] <- "lower"
colnames(modlist.500_Coeffs)[3] <- "mean"
colnames(modlist.500_Coeffs)[4] <- "upper"
modlist.500_Coeffs$Phytometer <- "all"
modlist.500_Coeffs$Response <- "Abundance"
modlist.500_Coeffs$Scale <- 500
modlist.500_Coeffs$Pollinator_group <- c(rep("allBees",4),
                                             rep("allBees_noApis",4),
                                             rep("Apis",4),
                                             rep("Bombus",4),
                                             rep("socialBees",4),
                                             rep("solitaryBees",4),
                                             rep("otherAculeata",4),
                                             rep("Syrphidae",4),
                                             rep("Coleoptera",4),
                                             rep("allPollinators",4))  
  
modlist.500_Coeffs$fixeff <- rep(c("Intercept","PlantS.1.z","Urban_500.1.z","PlantsxUrban"),10)

modlist.500_Coeffs

#-> Continue with script 10.b)
##################################################################################################################
##################################################################################################################

# For supp-mats

##################################################
###Run the models with the other spatial scales###
##################################################

modlist.allScales <- list(
  mod.1_500=glmer(A_allBees~ PlantS.1.z * Urban_500.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.2_500=glmer(A_allBees_noApis ~ PlantS.1.z * Urban_500.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.3_500=glmer(A_Apis ~ PlantS.1.z * Urban_500.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.4_500=glmer(A_Bombus ~ PlantS.1.z * Urban_500.1.z  +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.5_500=glmer(A_socialBees ~ PlantS.1.z * Urban_500.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.6_500=glmer(A_solitaryBees ~ PlantS.1.z * Urban_500.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.7_500=glmer(A_otherAculeata ~ PlantS.1.z * Urban_500.1.z  +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.8_500=glmer(A_Syrphidae ~ PlantS.1.z * Urban_500.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.9_500=glmer(A_Coleoptera ~ PlantS.1.z * Urban_500.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.10_500=glmer(A_allPollinators ~ PlantS.1.z * Urban_500.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.1_250=glmer(A_allBees~ PlantS.1.z * Urban_250.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.2_250=glmer(A_allBees_noApis ~ PlantS.1.z * Urban_250.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.3_250=glmer(A_Apis ~ PlantS.1.z * Urban_250.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.4_250=glmer(A_Bombus ~ PlantS.1.z * Urban_250.1.z  +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.5_250=glmer(A_socialBees ~ PlantS.1.z * Urban_250.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.6_250=glmer(A_solitaryBees ~ PlantS.1.z * Urban_250.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.7_250=glmer(A_otherAculeata ~ PlantS.1.z * Urban_250.1.z  +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.8_250=glmer(A_Syrphidae ~ PlantS.1.z * Urban_250.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.9_250=glmer(A_Coleoptera ~ PlantS.1.z * Urban_250.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.10_250=glmer(A_allPollinators ~ PlantS.1.z * Urban_250.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.1_100=glmer(A_allBees~ PlantS.1.z * Urban_100.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.2_100=glmer(A_allBees_noApis ~ PlantS.1.z * Urban_100.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.3_100=glmer(A_Apis ~ PlantS.1.z * Urban_100.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.4_100=glmer(A_Bombus ~ PlantS.1.z * Urban_100.1.z  +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.5_100=glmer(A_socialBees ~ PlantS.1.z * Urban_100.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.6_100=glmer(A_solitaryBees ~ PlantS.1.z * Urban_100.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.7_100=glmer(A_otherAculeata ~ PlantS.1.z * Urban_100.1.z  +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.8_100=glmer(A_Syrphidae ~ PlantS.1.z * Urban_100.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.9_100=glmer(A_Coleoptera ~ PlantS.1.z * Urban_100.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.10_100=glmer(A_allPollinators ~ PlantS.1.z * Urban_100.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.1_50=glmer(A_allBees~ PlantS.1.z * Urban_50.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.2_50=glmer(A_allBees_noApis ~ PlantS.1.z * Urban_50.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.3_50=glmer(A_Apis ~ PlantS.1.z * Urban_50.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.4_50=glmer(A_Bombus ~ PlantS.1.z * Urban_50.1.z  +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.5_50=glmer(A_socialBees ~ PlantS.1.z * Urban_50.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.6_50=glmer(A_solitaryBees ~ PlantS.1.z * Urban_50.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.7_50=glmer(A_otherAculeata ~ PlantS.1.z * Urban_50.1.z  +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.8_50=glmer(A_Syrphidae ~ PlantS.1.z * Urban_50.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.9_50=glmer(A_Coleoptera ~ PlantS.1.z * Urban_50.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
  mod.10_50=glmer(A_allPollinators ~ PlantS.1.z * Urban_50.1.z  + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
)

#Calculate R2s
modlist.allScales_Fit <- ldply(modlist.allScales, function(mod) {data.frame(#WAIC = round(WAIC(mod,nsim=10000)$WAIC2,2), 
                                                            AICc = round(AICc(mod),2), 
                                                            R2m = round(r.squaredGLMM(mod)[3],2), 
                                                            R2c = round(r.squaredGLMM(mod)[6],2),
                                                            DispPar = dispersion_glmer(mod), 
                                                            meanRF = round(mean(ranef(mod)$Id.fac[,1]),5), 
                                                            shapiroRANEF = round(shapiro.test(ranef(mod)$Id.fac[,1])$p.value,3), 
                                                            shapiroRESID = round(shapiro.test(resid(mod))$p.value,3))}) #model = deparse(formula(mod)) -> to plot the model formula
modlist.allScales_Fit
colnames(modlist.allScales_Fit)[1]<-"Model"
modlist.allScales_Fit$Model.bis <- c(rep(c("mod.1","mod.2","mod.3","mod.4","mod.5","mod.6","mod.7","mod.8","mod.9","mod.10"),4))
modlist.allScales_Fit$Scale <- c(rep(500,10),rep(250,10),rep(100,10),rep(50,10))
modlist.allScales_Fit

write.table(modlist.allScales_Fit, file = "results/abundance_models_model_fit_all_scales.txt", 
            append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


######################################################
#Draw Baysian conclusions for all models in one line!#
######################################################

##Takes time##!
modlist.allScales_Coeffs <- ldply(modlist.allScales, function(mod) {
  simres <- sim(mod, n.sim = nsim)@fixef
  qtiles <- apply(simres, 2, quantile, prob = c(0.025, 0.5, 0.975))
  out <- as.data.frame(t(qtiles))
  out$fixeff <- rownames(out)
  out
})

#modlist.allScales_Coeffs <- ldply(modlist.allScales, function(mod) {coeff.mod = (as.data.frame(t(as.data.frame(apply(sim(mod, n.sim=nsim)@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975))))))})

rownames(modlist.allScales_Coeffs)<-NULL
colnames(modlist.allScales_Coeffs)[1]<-"Model"
colnames(modlist.allScales_Coeffs)[2] <- "lower"
colnames(modlist.allScales_Coeffs)[3] <- "mean"
colnames(modlist.allScales_Coeffs)[4] <- "upper"

#modlist.allScales_Coeffs$fixeff<-rep(c(rep(c("Intercept","PlantS.1.z","Urban_500.1.z","PlantsxUrban"),4)),10) #Needs to be fixed!
modlist.allScales_Coeffs

modlist.allScales_Coeffs$Pollinator_group <- rep(c(rep("allBees",4),
                                                   rep("allBees_noApis",4),
                                                   rep("Apis",4),
                                                   rep("Bombus",4),
                                                   rep("socialBees",4),
                                                   rep("solitaryBees",4),
                                                   rep("otherAculeata",4),
                                                   rep("Syrphidae",4),
                                                   rep("Coleoptera",4),
                                                   rep("allPollinators",4)),4)

modlist.allScales_Coeffs$Model.bis <- rep(c(rep("mod.1",4),
                                                   rep("mod.2",4),
                                                   rep("mod.3",4),
                                                   rep("mod.4",4),
                                                   rep("mod.5",4),
                                                   rep("mod.6",4),
                                                   rep("mod.7",4),
                                                   rep("mod.8",4),
                                                   rep("mod.9",4),
                                                   rep("mod.10",4)),4)


modlist.allScales_Coeffs$Phytometer <- "all"
modlist.allScales_Coeffs$Response <- "Abundance"
modlist.allScales_Coeffs$Scale <- c(rep(500,40),rep(250,40),rep(100,40),rep(50,40))

modlist.allScales_Coeffs

write.table(modlist.allScales_Coeffs, file = "results/abundance_models_coefficents_all_scales.txt", 
            append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

###############################################################################################################
###############################################################################################################

##################################
## Plot the parameter estimates ##
##################################

df7 <-read.table("results/abundance_models_coefficents_all_scales.txt", header = TRUE, sep = ";")
head(df7)

#Plot the parameter estimates for each scale:

str(df7)
summary(df7$fixeff)

library(ggplot2)
library(dplyr)
library(stringr)

df_plot <- df7 %>%
  filter(fixeff %in% c("Urban_500.1.z", "Urban_250.1.z", "Urban_100.1.z", "Urban_50.1.z")) %>%
  mutate(
    Scale = as.factor(Scale),
    Pollinator_group = recode(Pollinator_group,
                              "Apis" = "Honeybees",
                              "Bombus" = "Bumblebees",
                              "socialBees" = "Social bees",
                              "solitaryBees" = "Solitary bees",
                              "Syrphidae" = "Hoverflies",
                              "Coleoptera" = "Beetles",
                              "otherAculeata" = "Wasps",
                              "allBees" = "All bees",
                              "allBees_noApis" = "All bees (no honeybees)",
                              "allPollinators" = "All pollinators"
    )
  ) %>%
  mutate(
    Pollinator_group = factor(Pollinator_group,
                              levels = c("All bees","All bees (no honeybees)",
                                         "Honeybees","Bumblebees","Social bees",
                                         "Solitary bees","Wasps","Hoverflies",
                                         "Beetles","All pollinators") ))
  

### Supplementary figure S1

ggplot(df_plot, aes(x = Scale, y = mean, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#E01A4F", linewidth = 0.8) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, linewidth = 0.8) +
  facet_wrap(~Pollinator_group, ncol = 5, scales  = "free_y") +
  labs(
    x = "Spatial scale (m)",
    y = "Mean parameter estimate and 95% CI",
    color = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "white", color = "black"),
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    panel.grid.major = element_blank(),   # remove major grid lines
    panel.grid.minor = element_blank()    # remove minor grid lines
  )

#12.5*6  inches landscape

# Supplementary Figure S2

library(MuMIn)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# Extract R2 from all models
r2_results <- lapply(modlist.allScales, function(mod) {
  r2 <- r.squaredGLMM(mod)
  tibble(
    R2m = r2[1],
    R2c = r2[2]
  )
}) %>%
  bind_rows(.id = "Model")

# Add metadata (pollinator group and scale) from model names
r2_results <- r2_results %>%
  mutate(
    Model_num = gsub("mod\\.(\\d+)_.*", "\\1", Model),
    Scale = as.numeric(gsub(".*_(\\d+)$", "\\1", Model)),
    Pollinator_group = recode(Model_num,
                              "1" = "All bees",
                              "2" = "All bees (no honeybees)",
                              "3" = "Honeybees",
                              "4" = "Bumblebees",
                              "5" = "Social bees",
                              "6" = "Solitary bees",
                              "7" = "Wasps",
                              "8" = "Hoverflies",
                              "9" = "Beetles",
                              "10" = "All pollinators"
    ),
    Pollinator_group = factor(Pollinator_group, levels = c(
      "All bees", "All bees (no honeybees)", "Honeybees", "Bumblebees",
      "Social bees", "Solitary bees", "Wasps", "Hoverflies", "Beetles", "All pollinators"
    ))
  )

r2_long <- r2_results %>%
  pivot_longer(cols = c(R2m, R2c), names_to = "Type", values_to = "R2") %>%
  mutate(
    Type = recode(Type,
                  "R2m" = "R²m",
                  "R2c" = "R²c"),
    Type = factor(Type, levels = c("R²m", "R²c"))
  ) %>%
  mutate(
    Scale = factor(Scale, levels = c(50, 100, 250, 500)),
    Pollinator_group = factor(Pollinator_group, levels = c(
      "All bees", "All bees (no honeybees)", "Honeybees", "Bumblebees",
      "Social bees", "Solitary bees", "Wasps", "Hoverflies",
      "Beetles", "All pollinators"
    ))
  )

r2_long$Scale <- factor(r2_long$Scale, levels = c(50, 100, 250, 500))

library(patchwork)

# Define group order and split into two
group_order <- c("All bees", "All bees (no honeybees)", "Honeybees", "Bumblebees", "Social bees",
                 "Solitary bees", "Wasps", "Hoverflies", "Beetles", "All pollinators")

r2_long$Pollinator_group <- factor(r2_long$Pollinator_group, levels = group_order)

# Subset into top 5 and bottom 5
df_top <- r2_long %>% filter(Pollinator_group %in% group_order[1:5])
df_bottom <- r2_long %>% filter(Pollinator_group %in% group_order[6:10])

# Create each plot
p_top <- ggplot(df_top, aes(x = Type, y = R2, fill = factor(Scale))) +
  geom_bar(stat = "identity", width = 0.6) +
  facet_grid(Scale ~ Pollinator_group, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(
    values = c("500" = "#22577a", 
               "250" = "#38a3a5", 
               "100" = "#57cc99", 
               "50" = "#80ed99")
  )+
  scale_y_continuous(limits = c(0, 0.7))+
  labs(
    x = "Spatial scale (m)",
    y = expression(paste("Explained variance (", R^2, ")")),
    fill = "Scale  (m)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    legend.position = "none",
    panel.grid.major = element_blank(),   # remove major grid lines
    panel.grid.minor = element_blank()    # remove minor grid lines
  )

p_bottom <- ggplot(df_bottom, aes(x = Type, y = R2, fill = factor(Scale))) +
  geom_bar(stat = "identity", width = 0.6) +
  facet_grid(Scale ~ Pollinator_group, scales = "free_y") +
  scale_y_continuous(limits = c(0, 0.7))+
  coord_flip() +
  scale_fill_manual(
    values = c("500" = "#22577a", 
               "250" = "#38a3a5", 
               "100" = "#57cc99", 
               "50" = "#80ed99")
  )+
  labs(
    x = "Spatial scale (m)",
    y = expression(paste("Explained variance (", R^2, ")")),
    fill = "Scale (m)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    legend.position = "bottom",
    panel.grid.major = element_blank(),   # remove major grid lines
    panel.grid.minor = element_blank()    # remove minor grid lines
  )

# Combine
p_top / p_bottom + plot_layout(heights = c(1, 1))

#  9 x 12 in  landscape
