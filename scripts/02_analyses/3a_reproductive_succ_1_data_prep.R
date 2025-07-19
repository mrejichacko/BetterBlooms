#####################################################################
#####################################################################
#####################################################################
###
### 3a. Pollinator experiment: 
### Modelling plant reproductive success: Load all data tables and scale the variables
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
#####################################################################
#####################################################################
#####################################################################
###
### Prep

#https://nicercode.github.io/guides/repeating-things/

rm(list=ls())
gc()

#The predictors
df1<-read.table("raw_data/plant_floristic_data.txt", header = TRUE, sep = ";")
df2<-read.table("raw_data/explanatory_variables.txt", header = TRUE, sep = ";")
df3<-read.table("cleaned_data/pollinator_abundance_aggregated_garden_phytometer.txt", header = TRUE, sep = ";")
df4<-read.table("cleaned_data/pollinator_SR_aggregated_garden_phytometer.txt", header = TRUE, sep = ";")
df5<-read.table("cleaned_data/sampling_effort_gardens_plants_aggregated.txt", header = TRUE, sep = ";")

#The reproductive success data
df6a<-read.table("raw_data/RejiChacko_etal_2025_EnviDat/07_pollination_success/daucus_carota_seed_set.csv", header = TRUE, sep = ",")
df6b<-read.table("raw_data/RejiChacko_etal_2025_EnviDat/07_pollination_success/raphanus_sativus_fruit_set.csv", header = TRUE, sep = ",")
df6c<-read.table("raw_data/RejiChacko_etal_2025_EnviDat/07_pollination_success/raphanus_sativus_seed_set.csv", header = TRUE, sep = ",")
df6d<-read.table("raw_data/RejiChacko_etal_2025_EnviDat/07_pollination_success/onobrychis_viciifolia_fruit_set.csv", header = TRUE, sep = ",")
df6e<-read.table("raw_data/RejiChacko_etal_2025_EnviDat/07_pollination_success/symphytum_officinale_fruit_set.csv", header = TRUE, sep = ",")
df6f<-read.table("raw_data/RejiChacko_etal_2025_EnviDat/07_pollination_success/symphytum_officinale_seed_set.csv", header = TRUE, sep = ",")

#Combine the datasets: predictors
datExpl <- Reduce(function(x, y) merge(x, y, by="Id", sort=TRUE), list(df1[,c(1,5)],df2[,c(1,4,5,13:16)],df3,df4,df5))
str(datExpl)

#Prepare the explanatory variables:

#1.Calculate the sampling effort

#sampling effort in fieldwork days : 1 fieldwork day = 9h (-> one fieldwork day lasted from 9:00 to 18:00)
datExpl$Dayly_sampling_effort_Carrot <- datExpl$sampling_effort_min_carrot/60/9 
datExpl$Dayly_sampling_effort_Radish <- datExpl$sampling_effort_min_radish/60/9 
datExpl$Dayly_sampling_effort_Sainfoin <- datExpl$sampling_effort_min_sainfoin/60/9 
datExpl$Dayly_sampling_effort_Comfrey <- datExpl$sampling_effort_min_comfrey/60/9

#2. Calculate dayly capture rates for each phytometer:

#2.1 Carrot 

#Scaled:
datExpl$A_allBees_Carrot.dayly.z <-  scale(datExpl$A_allBees_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$A_allBees_noApis_Carrot.dayly.z <-  scale(datExpl$A_allBees_noApis_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$A_Apis_Carrot.dayly.z <- scale(datExpl$A_Apis_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$A_Bombus_Carrot.dayly.z <- scale(datExpl$A_Bombus_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$A_socialBees_Carrot.dayly.z <- scale(datExpl$A_socialBees_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$A_solitaryBees_Carrot.dayly.z <- scale(datExpl$A_solitaryBees_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$A_otherAculeata_Carrot.dayly.z <- scale(datExpl$A_otherAculeata_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$A_Syrphidae_Carrot.dayly.z <- scale(datExpl$A_Syrphidae_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$A_Coleoptera_Carrot.dayly.z <- scale(datExpl$A_Coleoptera_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$A_allPollinators_Carrot.dayly.z <- scale(datExpl$A_allPollinators_Carrot/datExpl$Dayly_sampling_effort_Carrot)

datExpl$S_allBees_Carrot.dayly.z <-  scale(datExpl$S_allBees_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$S_allBees_noApis_Carrot.dayly.z <-  scale(datExpl$S_allBees_noApis_Carrot/datExpl$Dayly_sampling_effort_Carrot)
#datExpl$S_Apis_Carrot.dayly.z <- scale(datExpl$S_Apis_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$S_Bombus_Carrot.dayly.z <- scale(datExpl$S_Bombus_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$S_socialBees_Carrot.dayly.z <- scale(datExpl$S_socialBees_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$S_solitaryBees_Carrot.dayly.z <- scale(datExpl$S_solitaryBees_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$S_otherAculeata_Carrot.dayly.z <- scale(datExpl$S_otherAculeata_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$S_Syrphidae_Carrot.dayly.z <- scale(datExpl$S_Syrphidae_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$S_Coleoptera_Carrot.dayly.z <- scale(datExpl$S_Coleoptera_Carrot/datExpl$Dayly_sampling_effort_Carrot)
datExpl$S_allPollinators_Carrot.dayly.z <- scale(datExpl$S_allPollinators_Carrot/datExpl$Dayly_sampling_effort_Carrot)

#Unscaled:
datExpl$A_allBees_Carrot.dayly <-  datExpl$A_allBees_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$A_allBees_noApis_Carrot.dayly <-  datExpl$A_allBees_noApis_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$A_Apis_Carrot.dayly <- datExpl$A_Apis_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$A_Bombus_Carrot.dayly <- datExpl$A_Bombus_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$A_socialBees_Carrot.dayly <- datExpl$A_socialBees_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$A_solitaryBees_Carrot.dayly <- datExpl$A_solitaryBees_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$A_otherAculeata_Carrot.dayly <- datExpl$A_otherAculeata_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$A_Syrphidae_Carrot.dayly <- datExpl$A_Syrphidae_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$A_Coleoptera_Carrot.dayly <- datExpl$A_Coleoptera_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$A_allPollinators_Carrot.dayly <- datExpl$A_allPollinators_Carrot/datExpl$Dayly_sampling_effort_Carrot

datExpl$S_allBees_Carrot.dayly <-  datExpl$S_allBees_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$S_allBees_noApis_Carrot.dayly <-  datExpl$S_allBees_noApis_Carrot/datExpl$Dayly_sampling_effort_Carrot
#datExpl$S_Apis_Carrot.dayly <- datExpl$S_Apis_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$S_Bombus_Carrot.dayly <- datExpl$S_Bombus_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$S_socialBees_Carrot.dayly <- datExpl$S_socialBees_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$S_solitaryBees_Carrot.dayly <- datExpl$S_solitaryBees_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$S_otherAculeata_Carrot.dayly <- datExpl$S_otherAculeata_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$S_Syrphidae_Carrot.dayly <- datExpl$S_Syrphidae_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$S_Coleoptera_Carrot.dayly <- datExpl$S_Coleoptera_Carrot/datExpl$Dayly_sampling_effort_Carrot
datExpl$S_allPollinators_Carrot.dayly <- datExpl$S_allPollinators_Carrot/datExpl$Dayly_sampling_effort_Carrot

#2.2 Radish

#Scaled:
datExpl$A_allBees_Radish.dayly.z <-  scale(datExpl$A_allBees_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$A_allBees_noApis_Radish.dayly.z <-  scale(datExpl$A_allBees_noApis_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$A_Apis_Radish.dayly.z <- scale(datExpl$A_Apis_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$A_Bombus_Radish.dayly.z <- scale(datExpl$A_Bombus_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$A_socialBees_Radish.dayly.z <- scale(datExpl$A_socialBees_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$A_solitaryBees_Radish.dayly.z <- scale(datExpl$A_solitaryBees_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$A_otherAculeata_Radish.dayly.z <- scale(datExpl$A_otherAculeata_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$A_Syrphidae_Radish.dayly.z <- scale(datExpl$A_Syrphidae_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$A_Coleoptera_Radish.dayly.z <- scale(datExpl$A_Coleoptera_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$A_allPollinators_Radish.dayly.z <- scale(datExpl$A_allPollinators_Radish/datExpl$Dayly_sampling_effort_Radish)

datExpl$A_Coleoptera_Radish.dayly.t.z <- scale(log(sqrt(datExpl$A_Coleoptera_Radish/datExpl$Dayly_sampling_effort_Radish+1)))

#Scaled:
datExpl$S_allBees_Radish.dayly.z <-  scale(datExpl$S_allBees_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$S_allBees_noApis_Radish.dayly.z <-  scale(datExpl$S_allBees_noApis_Radish/datExpl$Dayly_sampling_effort_Radish)
#datExpl$S_Apis_Radish.dayly.z <- scale(datExpl$S_Apis_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$S_Bombus_Radish.dayly.z <- scale(datExpl$S_Bombus_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$S_socialBees_Radish.dayly.z <- scale(datExpl$S_socialBees_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$S_solitaryBees_Radish.dayly.z <- scale(datExpl$S_solitaryBees_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$S_otherAculeata_Radish.dayly.z <- scale(datExpl$S_otherAculeata_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$S_Syrphidae_Radish.dayly.z <- scale(datExpl$S_Syrphidae_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$S_Coleoptera_Radish.dayly.z <- scale(datExpl$S_Coleoptera_Radish/datExpl$Dayly_sampling_effort_Radish)
datExpl$S_allPollinators_Radish.dayly.z <- scale(datExpl$S_allPollinators_Radish/datExpl$Dayly_sampling_effort_Radish)



#Unscaled:
datExpl$A_allBees_Radish.dayly <-  datExpl$A_allBees_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$A_allBees_noApis_Radish.dayly <-  datExpl$A_allBees_noApis_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$A_Apis_Radish.dayly <- datExpl$A_Apis_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$A_Bombus_Radish.dayly <- datExpl$A_Bombus_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$A_socialBees_Radish.dayly <- datExpl$A_socialBees_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$A_solitaryBees_Radish.dayly <- datExpl$A_solitaryBees_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$A_otherAculeata_Radish.dayly <- datExpl$A_otherAculeata_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$A_Syrphidae_Radish.dayly <- datExpl$A_Syrphidae_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$A_Coleoptera_Radish.dayly <- datExpl$A_Coleoptera_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$A_allPollinators_Radish.dayly <- datExpl$A_allPollinators_Radish/datExpl$Dayly_sampling_effort_Radish

datExpl$S_allBees_Radish.dayly <-  datExpl$S_allBees_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$S_allBees_noApis_Radish.dayly <-  datExpl$S_allBees_noApis_Radish/datExpl$Dayly_sampling_effort_Radish
#datExpl$S_Apis_Radish.dayly <- datExpl$S_Apis_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$S_Bombus_Radish.dayly <- datExpl$S_Bombus_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$S_socialBees_Radish.dayly <- datExpl$S_socialBees_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$S_solitaryBees_Radish.dayly <- datExpl$S_solitaryBees_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$S_otherAculeata_Radish.dayly <- datExpl$S_otherAculeata_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$S_Syrphidae_Radish.dayly <- datExpl$S_Syrphidae_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$S_Coleoptera_Radish.dayly <- datExpl$S_Coleoptera_Radish/datExpl$Dayly_sampling_effort_Radish
datExpl$S_allPollinators_Radish.dayly <- datExpl$S_allPollinators_Radish/datExpl$Dayly_sampling_effort_Radish

#2.3 Sainfoin

#Scaled:
datExpl$A_allBees_Sainfoin.dayly.z <-  scale(datExpl$A_allBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
datExpl$A_allBees_noApis_Sainfoin.dayly.z <-  scale(datExpl$A_allBees_noApis_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
datExpl$A_Apis_Sainfoin.dayly.z <- scale(datExpl$A_Apis_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
datExpl$A_Bombus_Sainfoin.dayly.z <- scale(datExpl$A_Bombus_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
datExpl$A_socialBees_Sainfoin.dayly.z <- scale(datExpl$A_socialBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
datExpl$A_solitaryBees_Sainfoin.dayly.z <- scale(datExpl$A_solitaryBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
datExpl$A_allPollinators_Sainfoin.dayly.z <- scale(datExpl$A_allPollinators_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)

datExpl$S_allBees_Sainfoin.dayly.z <-  scale(datExpl$S_allBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
datExpl$S_allBees_noApis_Sainfoin.dayly.z <-  scale(datExpl$S_allBees_noApis_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
#datExpl$S_Apis_Sainfoin.dayly.z <- scale(datExpl$S_Apis_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
datExpl$S_Bombus_Sainfoin.dayly.z <- scale(datExpl$S_Bombus_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
datExpl$S_socialBees_Sainfoin.dayly.z <- scale(datExpl$S_socialBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
datExpl$S_solitaryBees_Sainfoin.dayly.z <- scale(datExpl$S_solitaryBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)
datExpl$S_allPollinators_Sainfoin.dayly.z <- scale(datExpl$S_allPollinators_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin)

#Unscaled:
datExpl$A_allBees_Sainfoin.dayly <-  datExpl$A_allBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
datExpl$A_allBees_noApis_Sainfoin.dayly <-  datExpl$A_allBees_noApis_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
datExpl$A_Apis_Sainfoin.dayly <- datExpl$A_Apis_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
datExpl$A_Bombus_Sainfoin.dayly <- datExpl$A_Bombus_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
datExpl$A_socialBees_Sainfoin.dayly <- datExpl$A_socialBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
datExpl$A_solitaryBees_Sainfoin.dayly <- datExpl$A_solitaryBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
datExpl$A_allPollinators_Sainfoin.dayly <- datExpl$A_allPollinators_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin

datExpl$S_allBees_Sainfoin.dayly <-  datExpl$S_allBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
datExpl$S_allBees_noApis_Sainfoin.dayly <-  datExpl$S_allBees_noApis_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
#datExpl$S_Apis_Sainfoin.dayly <- datExpl$S_Apis_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
datExpl$S_Bombus_Sainfoin.dayly <- datExpl$S_Bombus_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
datExpl$S_socialBees_Sainfoin.dayly <- datExpl$S_socialBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
datExpl$S_solitaryBees_Sainfoin.dayly <- datExpl$S_solitaryBees_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin
datExpl$S_allPollinators_Sainfoin.dayly <- datExpl$S_allPollinators_Sainfoin/datExpl$Dayly_sampling_effort_Sainfoin

#2.4 Comfrey

#Scaled:
datExpl$A_Bombus_Comfrey.dayly.z <- scale(datExpl$A_Bombus_Comfrey/datExpl$Dayly_sampling_effort_Comfrey)

datExpl$S_Bombus_Comfrey.dayly.z <- scale(datExpl$S_Bombus_Comfrey/datExpl$Dayly_sampling_effort_Comfrey)


#Unscaled:
datExpl$A_Bombus_Comfrey.dayly <- datExpl$A_Bombus_Comfrey/datExpl$Dayly_sampling_effort_Comfrey

datExpl$S_Bombus_Comfrey.dayly.z <- scale(datExpl$S_Bombus_Comfrey/datExpl$Dayly_sampling_effort_Comfrey)

#3. Landscape
names(datExpl)
#datExpl$Urban_30.z <- scale(datExpl$Urban_30)
datExpl$Urban_50.z <- scale(datExpl$Urban_50)
datExpl$Urban_100.z <- scale(datExpl$Urban_100)
datExpl$Urban_250.z <- scale(datExpl$Urban_250)
datExpl$Urban_500.z <- scale(datExpl$Urban_500)

#4. Garden
datExpl$PlantS.z <- scale(datExpl$SR_all_insect_pollinated_May_August)

if (!dir.exists("environments")){
  dir.create("environments")
}else{
  print("dir exists")
}
save.image(file = "environments/3a_environment.RData")

# CONTINUE WITH ANALYSIS
#####################################################################################################################
#####################################################################################################################


