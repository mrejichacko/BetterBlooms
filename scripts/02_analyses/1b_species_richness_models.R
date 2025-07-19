#####################################################################
#####################################################################
#####################################################################
###
### 1b. Modelling SR: Final Models
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
### 
#####################################################################
#####################################################################
#####################################################################


######################!!!!!! RUN 1a FIRST!!!!!######################

##########################################################################
## First step: Make a list of the selected models and run them together ##
##########################################################################

# (a) The GLMM models:
modlist.2a <- list(mod.17=glmer(S_otherAculeata ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                   mod.19=glmer(S_Coleoptera ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                   mod.20=glmer(S_allPollinators ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d)))

#run them also outside the list:
mod.17<-glmer(S_otherAculeata ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z + (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.19<-glmer(S_Coleoptera ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.20<-glmer(S_allPollinators ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z +  (1|Id.fac), data=df, family="poisson", offset=log(Total_sampling_effort_d))

modlist.2a_Fit <- ldply(modlist.2a, function(mod) {data.frame(#WAIC = round(WAIC(mod,nsim=1000)$WAIC2,2), 
                                                            AICc = round(AICc(mod),2), 
                                                            #R2m = round(r.squaredGLMM(mod)[3],2), 
                                                            #R2c = round(r.squaredGLMM(mod)[6],2),
                                                            DispPar = dispersion_glmer(mod), 
                                                            meanRF = round(mean(ranef(mod)$Id.fac[,1]),5), 
                                                            shapiroRANEF = round(shapiro.test(ranef(mod)$Id.fac[,1])$p.value,3), 
                                                            shapiroRESID = round(shapiro.test(resid(mod))$p.value,3))})
colnames(modlist.2a_Fit)[1]<-"Model"
modlist.2a_Fit 

#Check also the model output:
summary(mod.17)
summary(mod.19)
summary(mod.20)

# (b) The GLM models:
mod.11<-glm(S_allBees~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.12<-glm(S_allBees_noApis ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.14<-glm(S_Bombus ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.15<-glm(S_socialBees ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.16<-glm(S_solitaryBees ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d))
mod.18<-glm(S_Syrphidae ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d))

summary(mod.11)#ok/ slightly underdisp. 
summary(mod.12)#ok/ slightly underdisp.
summary(mod.14)#ok; all NS
summary(mod.15)#massively underdispersed; a quasi-poisson could be a solution (increases the p-value)
summary(mod.16)#perfect
summary(mod.18)#perfect!

##Make another list:
modlist.2b <- list(mod.11=glm(S_allBees~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                 mod.12=glm(S_allBees_noApis ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                 mod.14=glm(S_Bombus ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                 mod.15=glm(S_socialBees ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                 mod.16=glm(S_solitaryBees ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d)),
                 mod.18=glm(S_Syrphidae ~ PlantS.1.z + Urban_500.1.z  + PlantS.1.z:Urban_500.1.z, data=df, family="poisson", offset=log(Total_sampling_effort_d)))

modlist.2b_Fit <- ldply(modlist.2b, function(mod) {data.frame(#WAIC = round(WAIC(mod, nsim=1000)$WAIC2,2),
                                                            #R2 = round(rsquared(mod)[5],2), piecewiseSEM-package
                                                            shapiroRESID = round(shapiro.test(resid(mod))$p.value,3))})

colnames(modlist.2b_Fit)[1]<-"Model"
modlist.2b_Fit 

###################################################################
## Second step: check model assumptions with the selected models ##
###################################################################

#Check the residuals:

par(mfrow=c(6,5),
    oma = c(1,4.8,4.8,0) + 0.1,
    mar = c(1,0,1,2.9) + 0.1,
    cex.lab = 2.5)

clab = 0.8
cmain = 0.8
caxis = 0.8

#Model 11 
plot(fitted(mod.11), resid(mod.11), main="Residuals vs. fitted", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.11), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main="Q-Q plot residuals") #
qqline(resid(mod.11))
#qqnorm(ranef(mod.11)$Id.fac[,1], main="Q-Q plot RF: Garden", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
#qqline(ranef(mod.11)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.11), main="Residuals vs. Predictor 1", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.11), main="Residuals vs. Predictor 2", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.11), df$A_allBees, main="Observed vs. fitted", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

#Model 12 
plot(fitted(mod.12), resid(mod.12), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.12), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.12))
#qqnorm(ranef(mod.12)$Id.fac[,1], cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
#qqline(ranef(mod.12)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.12), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.12), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.12), df$A_allBees_noApis, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

#Model 14
plot(fitted(mod.14), resid(mod.14), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.14), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.14))
#qqnorm(ranef(mod.14)$Id.fac[,1], cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
#qqline(ranef(mod.14)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.14), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.14), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.14), df$A_Bombus, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted
#dev.off()

#Model 15
plot(fitted(mod.15), resid(mod.15), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.15),  cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.15))
#qqnorm(ranef(mod.15)$Id.fac[,1],cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
#qqline(ranef(mod.15)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.15), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.15), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.15), df$A_socialBees, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

#Model 16
plot(fitted(mod.16), resid(mod.16), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.16), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.16))
#qqnorm(ranef(mod.16)$Id.fac[,1], cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
#qqline(ranef(mod.16)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.16), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.16), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.16), df$A_solitaryBees, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

#Model 18 
plot(fitted(mod.18), resid(mod.18), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.18), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.18))
#qqnorm(ranef(mod.18)$Id.fac[,1],  cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
#qqline(ranef(mod.18)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.18),  cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.18), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.18), df$A_Syrphidae, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

##################
#The mixed models#
##################

#tiff(filename = "20190117_GLMMs_SpeRich_resid.tif",width = 22.5, height = 11.5, units = "cm", res = 250, compression = c("none"))

par(mfrow=c(3,6),
    oma = c(1,4.8,4.8,0) + 0.1,
    mar = c(1,0,1,2.9) + 0.1,
    cex.lab = 2.5)

clab = 0.8
cmain = 0.8
caxis = 0.8

#Model 17 
plot(fitted(mod.17), resid(mod.17), main="Residuals vs. fitted",  cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.17), main="Q-Q plot residuals", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(resid(mod.17))
qqnorm(ranef(mod.17)$Id.fac[,1],main="Q-Q plot RF: Garden", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #
qqline(ranef(mod.17)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.17), main="Residuals vs. Predictor 1", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.17), main="Residuals vs. Predictor 2", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.17), df$A_otherAculeata, main="Observed vs. fitted", cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

#Model 19 
plot(fitted(mod.19), resid(mod.19), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.19), cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.19))
qqnorm(ranef(mod.19)$Id.fac[,1], cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(ranef(mod.19)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.19), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.19), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.19), df$A_Coleoptera, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

#Model 10 
plot(fitted(mod.20), resid(mod.20),xlab = "Fitted values", cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")#fitted vs. residuals
qqnorm(resid(mod.20), cex.lab=clab,  xlab="Theoretical quantiles", cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(resid(mod.20))
qqnorm(ranef(mod.20)$Id.fac[,1], cex.lab=clab, cex.main=cmain, cex.axis=caxis, main=NULL) #
qqline(ranef(mod.20)$Id.fac[,1]) #
plot(df$PlantS.1.z, resid(mod.20), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
plot(df$Urban_500.1.z, resid(mod.20), cex.lab=clab, cex.main=cmain, cex.axis=caxis); abline(0,0, lty=2, col="red")
scatter.smooth(fitted(mod.20), df$A_allPollinators, cex.lab=clab, cex.main=cmain, cex.axis=caxis) #data vs. fitted

##########################################################################
#Spatial correlation: make one big multipanel plot including all figures:#
##########################################################################

library(sp)

#Model 11:
spdata <- data.frame(resid=resid(mod.11), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata) <- c("x","y")
bubble(spdata, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 12:
spdata <- data.frame(resid=resid(mod.12), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata) <- c("x","y")
bubble(spdata, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 14:
spdata <- data.frame(resid=resid(mod.14), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata) <- c("x","y")
bubble(spdata, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 15:
spdata <- data.frame(resid=resid(mod.15), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata) <- c("x","y")
bubble(spdata, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 16:
spdata <- data.frame(resid=resid(mod.16), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata) <- c("x","y")
bubble(spdata, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 17:
spdata <- data.frame(resid=resid(mod.17), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata) <- c("x","y")
bubble(spdata, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 18: there might be some structure. 
spdata <- data.frame(resid=resid(mod.18), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata) <- c("x","y")
bubble(spdata, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 19:. 
spdata <- data.frame(resid=resid(mod.19), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata) <- c("x","y")
bubble(spdata, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")

#Model 20: 
spdata <- data.frame(resid=resid(mod.20), x=df$X_KOORDINATE, y=df$Y_KOORDINATE) 
coordinates(spdata) <- c("x","y")
bubble(spdata, "resid", col=c("blue", "orange"), main="Residuals", xlab="X-coordinates", ylab="Y-coordinates")


#####################################################################################
#Third step: Draw Baysian conclusions: Extract the coefficients for tables and plots#
#####################################################################################

library(arm)
nsim <- 10000

#Run them separately for the plots (you need the bsim.mod-objects!)
bsim.mod.11 <- sim(mod.11, n.sim=nsim)
bsim.mod.12 <- sim(mod.12, n.sim=nsim)
bsim.mod.14 <- sim(mod.14, n.sim=nsim)
bsim.mod.15 <- sim(mod.15, n.sim=nsim)
bsim.mod.16 <- sim(mod.16, n.sim=nsim)
bsim.mod.17 <- sim(mod.17, n.sim=nsim)
bsim.mod.18 <- sim(mod.18, n.sim=nsim)
bsim.mod.19 <- sim(mod.19, n.sim=nsim)
bsim.mod.20 <- sim(mod.20, n.sim=nsim)

#The glmmms
modlist.500_Coeffs_S.a <- ldply(modlist.2a, function(mod) {coeff.mod = (as.data.frame(t(as.data.frame(apply(sim(mod, n.sim=nsim)@fixef, 2, quantile, prob=c(0.025, 0.5, 0.975))))))})

rownames(modlist.500_Coeffs_S.a)<-NULL
colnames(modlist.500_Coeffs_S.a)[1]<-"Model"
colnames(modlist.500_Coeffs_S.a)[2] <- "lower"
colnames(modlist.500_Coeffs_S.a)[3] <- "mean"
colnames(modlist.500_Coeffs_S.a)[4] <- "upper"
modlist.500_Coeffs_S.a$Phytometer <- "all"
modlist.500_Coeffs_S.a$Response <- "SpeciesRichness"
modlist.500_Coeffs_S.a$Scale <- 500
modlist.500_Coeffs_S.a$Pollinator_group <- c(rep("otherAculeata",4),
                                         rep("Coleoptera",4),
                                         rep("allPollinators",4))  
modlist.500_Coeffs_S.a$fixeff <- rep(c("Intercept","PlantS.1.z","Urban_500.1.z","PlantsxUrban"),3)

#The glms
modlist.500_Coeffs_S.b <- ldply(modlist.2b, function(mod) {coeff.mod = (as.data.frame(t(as.data.frame(apply(sim(mod, n.sim=nsim)@coef, 2, quantile, prob=c(0.025, 0.5, 0.975))))))})
modlist.500_Coeffs_S.b

rownames(modlist.500_Coeffs_S.b)<-NULL
colnames(modlist.500_Coeffs_S.b)[1]<-"Model"
colnames(modlist.500_Coeffs_S.b)[2] <- "lower"
colnames(modlist.500_Coeffs_S.b)[3] <- "mean"
colnames(modlist.500_Coeffs_S.b)[4] <- "upper"
modlist.500_Coeffs_S.b$Phytometer <- "all"
modlist.500_Coeffs_S.b$Response <- "SpeciesRichness"
modlist.500_Coeffs_S.b$Scale <- 500
modlist.500_Coeffs_S.b$Pollinator_group <- c(rep("allBees",4),
                                             rep("allBees_noApis",4),
                                             rep("Bombus",4),
                                             rep("socialBees",4),
                                             rep("solitaryBees",4),
                                             rep("Syrphidae",4))  


modlist.500_Coeffs_S.b$fixeff <- rep(c("Intercept","PlantS.1.z","Urban_500.1.z","PlantsxUrban"),6)

#Combine all three tables:
coeff_comb <- rbind(modlist.500_Coeffs, modlist.500_Coeffs_S.a, modlist.500_Coeffs_S.b)
coeff_comb

write.table(coeff_comb, file = "results/beta_coefficients_all_plants.txt", append = FALSE, quote = TRUE, sep = " ;",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

#->Continue with script 1c)

save.image(file = "environments/1a_1b_environment.RData")


######################################################################################################################################################
########################################################################################################################