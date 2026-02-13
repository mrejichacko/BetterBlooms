#####################################################################
#####################################################################
#####################################################################
###
### 4b. Pollinator experiment: 
### Modelling plant reproductive success: Radish fruit set ~ Urban
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
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

df7b <- merge(datExpl, df6b)

# Have a look at the fruit set data: for convenience, they were counted branch-wise, 
# but we could aggregate (sum) them at the level of entire plants

#table(df6$Id, df6$Plant_Nr) #All gardens have six plants with varying number of branches
#names(df6)
#dim(df6)

#Aggregate:
#dat1 <- aggregate(cbind(n_flowers_with_fruits, n_flowers_without_fruits) ~ Id + Plant_Id, data=df6, function(x)sum(x))
#head(dat1)
#dim(dat1)

str(df7b)
names(df7b)

#Kick out garden Nr. 39
df7b <- droplevels(df7b[which(df7b$Id != 39),])
dim(df7b)

#Transform factor variables
df7b$Id.fac <- as.factor(df7b$Id)
df7b$Plant_Id.fac <- as.factor(df7b$plant_id)
df7b$Branch_Id.fac <- as.factor(df7b$branch_id)

str(df7b)
summary(df7b$A_Coleoptera_Radish.dayly.t.z)

###################
## Run the model ##
###################

#Abundance model: a binomial glmer (failure vrs. success); with dayly rates!

mod.1_Radish <- glmer(cbind(n_flowers_with_fruits, n_flowers_without_fruits) ~ 
                        PlantS.z + Urban_500.z  + PlantS.z:Urban_500.z + (1|Id.fac/Plant_Id.fac), data=df7b, family="binomial") #
summary(mod.1_Radish) 

r.squaredGLMM(mod.1_Radish)

library(car)
vif(mod.1_Radish) #No extreme values: >5 would be problematic

#Assessing model assumptions: (check page 144)
par(mfrow=c(2,3))

#fitted vs. residuals
scatter.smooth(fitted(mod.1_Radish), resid(mod.1_Radish)) #this looks quite good (no strong positive correlation = no strong shrinkage)
abline(h=0, lty=2, col="red") #Some values do not fit the data well (very large residuals) -> for large seed sets the fitted values are not very reliable. 
title("Tukey-Anscombe Plot")

#qq-plot of residuals
qqnorm(resid(mod.1_Radish)) #qq-plot of residuals
qqline(resid(mod.1_Radish))

#qq-plot  of random effects: Id
qqnorm(ranef(mod.1_Radish)$Id.fac[,1])
qqline(ranef(mod.1_Radish)$Id.fac[,1])

#qq-plot of the random effects: Plant
qqnorm(ranef(mod.1_Radish)$Plant_Id.fac[,1]) 
qqline(ranef(mod.1_Radish)$Plant_Id.fac[,1])

#fitted vs. observed
plot(fitted(mod.1_Radish), jitter(df7b$n_flowers_with_fruits/(df7b$n_flowers_with_fruits + df7b$n_flowers_without_fruits),0.05))
abline(0,1) 

dev.off()  

#Heteroscedasticity (= nonhomogeneity of the residual variance)
par(mfrow=c(1,2))
scatter.smooth(df7b$PlantS.z, resid(mod.1_Radish)); abline(0,0, lty=2, col="red")
scatter.smooth(df7b$Urban_500.z, resid(mod.1_Radish)); abline(0,0, lty=2, col="red")
dev.off()

#Check if random effects are 0:
mean(ranef(mod.1_Radish)$Id.fac[,1])
mean(ranef(mod.1_Radish)$Plant_Id.fac[,1])

#Check for overdispersion:
dispersion_glmer(mod.1_Radish) #!

# check for spatial autocorrelation
library(DHARMa)

sim_radish <- simulateResiduals(mod.1_Radish)

# Aggregate to unique spatial locations (garden = Id.fac)
sim_radish_garden <- recalculateResiduals(sim_radish, group = df7b$Id.fac)

# (unique x/y per garden)
sp_test_radish <- testSpatialAutocorrelation(
  sim_radish_garden,
  x = df7b$X_KOORDINATE[!duplicated(df7b$Id.fac)],
  y = df7b$Y_KOORDINATE[!duplicated(df7b$Id.fac)]
)

sp_test_radish
sp_test_radish$p.value
sp_test_radish$statistic

table_s4 <- data.frame(Predictor = "Habitat",
                       Response = "Radish fruit set", 
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

# let's check the results of a spatial model

library(spdep)
library(adespatial)

# coordinates at garden level
coords <- as.matrix(
  df7b[!duplicated(df7b$Id.fac), c("X_KOORDINATE", "Y_KOORDINATE")]
)

# k-nearest neighbours (robust with n = 23)
knn <- knearneigh(coords, k = 4)
nb  <- knn2nb(knn)
lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Moran Eigenvector Map (positive autocorrelation)
mem <- scores.listw(lw, MEM.autocor = "positive")

# attach MEM1 to full dataset (by garden)
df7b$MEM1 <- mem[match(df7b$Id.fac,
                       df7b$Id.fac[!duplicated(df7b$Id.fac)]), 1]

mod.1_Radish_MEM <- glmer(
  cbind(n_flowers_with_fruits, n_flowers_without_fruits) ~
    PlantS.z + Urban_500.z  + PlantS.z:Urban_500.z +
    MEM1 +
    (1|Id.fac/Plant_Id.fac),
  data = df7b,
  family = binomial
)

summary(mod.1_Radish_MEM)
summary(mod.1_Radish)

## fixed effects comparison (Estimate + p only; drop intercept)
get_fixef_tab <- function(mod, model_label, response_label) {
  cf <- as.data.frame(summary(mod)$coefficients)
  cf$Term <- rownames(cf)
  rownames(cf) <- NULL
  names(cf)[1:4] <- c("Estimate", "SE", "z", "p_value")
  
  out <- data.frame(
    Response = response_label,
    Model    = model_label,
    Term     = cf$Term,
    Estimate = round(cf$Estimate, 3),
    p_value  = signif(cf$p_value, 3),
    stringsAsFactors = FALSE
  )
  
  out[out$Term != "(Intercept)", ]
}

tab_radish_coef <- rbind(
  get_fixef_tab(mod.1_Radish,     "Original", "Radish"),
  get_fixef_tab(mod.1_Radish_MEM, "Original + MEM1",     "Radish")
)

tab_radish_coef <- tab_radish_coef[order(tab_radish_coef$Term), ]

tab_radish_coef







# DHARMa residuals
sim_radish_mem <- simulateResiduals(mod.1_Radish_MEM)

# aggregate to garden level
sim_radish_mem_garden <- recalculateResiduals(
  sim_radish_mem,
  group = df7b$Id.fac
)

# unique coordinates
xy <- df7b[!duplicated(df7b$Id.fac), c("X_KOORDINATE", "Y_KOORDINATE")]

# Moran’s I
sp_test_radish_mem <- testSpatialAutocorrelation(
  sim_radish_mem_garden,
  x = xy$X_KOORDINATE,
  y = xy$Y_KOORDINATE
)

sp_test_radish_mem
sp_test_radish_mem$p.value

table_s4 <- data.frame(Predictor = "Habitat + MEM1",
                       Response = "Radish fruit set", 
                       Observed = unname(sp_test_radish_mem$statistic[1]),
                       Expected = unname(sp_test_radish_mem$statistic[2]),
                       SD = unname(sp_test_radish_mem$statistic[3]),
                       P_value = sp_test_radish_mem$p.value
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
