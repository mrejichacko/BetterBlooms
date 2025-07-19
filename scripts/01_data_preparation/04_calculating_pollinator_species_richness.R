#####################################################################
#####################################################################
#####################################################################
###
### 4.Calculate pollinator SR
### Code by: David Frey and Merin Reji Chacko
### contact: merin.rejichacko@gmail.com
### Last edited: 18/07/2025
### 
### 
#####################################################################
#####################################################################
#####################################################################
###
### Prep

rm(list = ls())

df1 <-read.csv("raw_data/RejiChacko_etal_2025_EnviDat/06_trait_data/individual_traits.csv", header = TRUE, sep = ",")
checklist <- read.csv("raw_data/RejiChacko_etal_2025_EnviDat/04_taxonomic_data/taxa_checklist.csv")

df1 <- merge(df1, checklist, by = "taxon")
rm(checklist)

####Prepare the dataset like in 07
df2 <- df1[df1$family!="",]
df2<- df2[df2$genus!="",]

###############################
#### All plants per garden ####
###############################

#How many taxa do we have?
unique(df2$taxon)#163!
sort(table(df2$taxon), decreasing = TRUE) #List of observations per species.  

#All bees (Anthophila)
df3 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df3)[2] <- "S_allBees"
df3[order(df3$S_allBees),] 

#All bees (Anthophila) without honeybees
df4 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & taxon != "Apis mellifera" & nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df4)[2] <- "S_allBees_noApis"
df4[order(df4$S_allBees_noApis),] 

#Honeybees
#df5 <- aggregate(taxon ~ Id, data=df2, subset = taxon == "Apis mellifera", FUN=function(x)length(unique(x)))
#colnames(df5)[2] <- "S_Apis"
#df5[order(df5$S_Apis),] 

#Bumblebees
df6 <- aggregate(taxon ~ Id, data=df2, subset = genus == "Bombus" & nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df6)[2] <- "S_Bombus"
df6[order(df6$S_Bombus),] 

#Social bees (without honeybees)
df7 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "social_Bees" & taxon != "Apis mellifera"& nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df7)[2] <- "S_socialBees"
df7[order(df7$S_socialBees),] 

#Solitary bees (without honeybees)
df8 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "solitary_Bees" & taxon != "Apis mellifera"& nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df8)[2] <- "S_solitaryBees"
df8[order(df8$S_solitaryBees),] 

#Other Aculeata
df9 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Wasps", FUN=function(x)length(unique(x)))
colnames(df9)[2] <- "S_otherAculeata"
df9[order(df9$S_otherAculeata),]

#Hoverflies
df10 <- aggregate(taxon ~ Id, data=df2, subset = family == "Syrphidae"& nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df10)[2] <- "S_Syrphidae"
df10[order(df10$S_Syrphidae),]

#Beetles
df11 <- aggregate(taxon ~ Id, data=df2, subset = order == "Coleoptera", FUN=function(x)length(unique(x)))
colnames(df11)[2] <- "S_Coleoptera"
df11[order(df11$S_Coleoptera),]

#All pollinators
df12 <- aggregate(taxon ~ Id, data=df2, subset = nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df12)[2] <- "S_allPollinators"
df12[order(df12$S_allPollinators),]

###################################
#### Per phytometer per garden ####
###################################

####################
#### Carrot ########
####################

#All bees (Anthophila)
df13 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & phytometer_plant == "Carrot", FUN=function(x)length(unique(x)))
colnames(df13)[2] <- "S_allBees_Carrot"
df13[order(df13$S_allBees_Carrot),] 

#All bees (Anthophila) without honeybees
df14 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & taxon != "Apis mellifera" & phytometer_plant == "Carrot", FUN=function(x)length(unique(x)))
colnames(df14)[2] <- "S_allBees_noApis_Carrot"
df14[order(df14$S_allBees_noApis_Carrot),] 

#Honeybees
#df15 <- aggregate(taxon ~ Id, data=df2, subset = taxon == "Apis mellifera" & phytometer_plant == "Carrot", FUN=function(x)length(unique(x)))
#colnames(df15)[2] <- "S_Apis_Carrot"
#df15[order(df15$S_Apis_Carrot),] 

#Bumblebees
df16 <- aggregate(taxon ~ Id, data=df2, subset = genus == "Bombus" & phytometer_plant == "Carrot", FUN=function(x)length(unique(x)))
colnames(df16)[2] <- "S_Bombus_Carrot"
df16[order(df16$S_Bombus_Carrot),] 

#Social bees (without honeybees)
df17 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "social_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Carrot", FUN=function(x)length(unique(x)))
colnames(df17)[2] <- "S_socialBees_Carrot"
df17[order(df17$S_socialBees_Carrot),] 

#Solitary bees (without honeybees)
df18 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "solitary_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Carrot", FUN=function(x)length(unique(x)))
colnames(df18)[2] <- "S_solitaryBees_Carrot"
df18[order(df18$S_solitaryBees_Carrot),] 

#Other Aculeata
df19 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Wasps" & phytometer_plant == "Carrot", FUN=function(x)length(unique(x)))
colnames(df19)[2] <- "S_otherAculeata_Carrot"
df19[order(df19$S_otherAculeata_Carrot),]

#Hoverflies
df20 <- aggregate(taxon ~ Id, data=df2, subset = family == "Syrphidae" & phytometer_plant == "Carrot", FUN=function(x)length(unique(x)))
colnames(df20)[2] <- "S_Syrphidae_Carrot"
df20[order(df20$S_Syrphidae_Carrot),]

#Beetles
df21 <- aggregate(taxon ~ Id, data=df2, subset = order == "Coleoptera" & phytometer_plant == "Carrot", FUN=function(x)length(unique(x)))
colnames(df21)[2] <- "S_Coleoptera_Carrot"
df21[order(df21$S_Coleoptera_Carrot),]

#All pollinators
df22 <- aggregate(taxon ~ Id, data=df2, subset = phytometer_plant == "Carrot", FUN=function(x)length(unique(x)))
colnames(df22)[2] <- "S_allPollinators_Carrot"
df22[order(df22$S_allPollinators_Carrot),]

######################
###### Radish ########
######################

#All bees (Anthophila)
df23 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & phytometer_plant == "Radish", FUN=function(x)length(unique(x)))
colnames(df23)[2] <- "S_allBees_Radish"
df23[order(df23$S_allBees_Radish),] 

#All bees (Anthophila) without honeybees
df24 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & taxon != "Apis mellifera" & phytometer_plant == "Radish", FUN=function(x)length(unique(x)))
colnames(df24)[2] <- "S_allBees_noApis_Radish"
df24[order(df24$S_allBees_noApis_Radish),] 

#Honeybees
#df25 <- aggregate(taxon ~ Id, data=df2, subset = taxon == "Apis mellifera" & phytometer_plant == "Radish", FUN=function(x)length(unique(x)))
#colnames(df25)[2] <- "S_Apis_Radish"
#df25[order(df25$S_Apis_Radish),] 

#Bumblebees
df26 <- aggregate(taxon ~ Id, data=df2, subset = genus == "Bombus" & phytometer_plant == "Radish", FUN=function(x)length(unique(x)))
colnames(df26)[2] <- "S_Bombus_Radish"
df26[order(df26$S_Bombus_Radish),] 

#Social bees (without honeybees)
df27 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "social_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Radish", FUN=function(x)length(unique(x)))
colnames(df27)[2] <- "S_socialBees_Radish"
df27[order(df27$S_socialBees_Radish),] 

#Solitary bees (without honeybees)
df28 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "solitary_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Radish", FUN=function(x)length(unique(x)))
colnames(df28)[2] <- "S_solitaryBees_Radish"
df28[order(df28$S_solitaryBees_Radish),] 

#Other Aculeata
df29 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Wasps" & phytometer_plant == "Radish", FUN=function(x)length(unique(x)))
colnames(df29)[2] <- "S_otherAculeata_Radish"
df29[order(df29$S_otherAculeata_Radish),]

#Hoverflies
df30 <- aggregate(taxon ~ Id, data=df2, subset = family == "Syrphidae" & phytometer_plant == "Radish", FUN=function(x)length(unique(x)))
colnames(df30)[2] <- "S_Syrphidae_Radish"
df30[order(df30$S_Syrphidae_Radish),]

#Beetles
df31 <- aggregate(taxon ~ Id, data=df2, subset = order == "Coleoptera" & phytometer_plant == "Radish", FUN=function(x)length(unique(x)))
colnames(df31)[2] <- "S_Coleoptera_Radish"
df31[order(df31$S_Coleoptera_Radish),]

#All pollinators
df32 <- aggregate(taxon ~ Id, data=df2, subset = phytometer_plant == "Radish", FUN=function(x)length(unique(x)))
colnames(df32)[2] <- "S_allPollinators_Radish"
df32[order(df32$S_allPollinators_Radish),]

########################
###### Sainfoin ########
########################

#All bees (Anthophila)
df33 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & phytometer_plant == "Sainfoin" & nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df33)[2] <- "S_allBees_Sainfoin"
df33[order(df33$S_allBees_Sainfoin),] 

#All bees (Anthophila) without honeybees
df34 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & taxon != "Apis mellifera" & phytometer_plant == "Sainfoin"& nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df34)[2] <- "S_allBees_noApis_Sainfoin"
df34[order(df34$S_allBees_noApis_Sainfoin),] 

#Honeybees
#df35 <- aggregate(taxon ~ Id, data=df2, subset = taxon == "Apis mellifera" & phytometer_plant == "Sainfoin"& nectar_robber != 1, FUN=function(x)length(unique(x)))
#colnames(df35)[2] <- "S_Apis_Sainfoin"
#df35[order(df35$S_Apis_Sainfoin),] 

#Bumblebees
df36 <- aggregate(taxon ~ Id, data=df2, subset = genus == "Bombus" & phytometer_plant == "Sainfoin"& nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df36)[2] <- "S_Bombus_Sainfoin"
df36[order(df36$S_Bombus_Sainfoin),] 

#Social bees (without honeybees)
df37 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "social_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Sainfoin"& nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df37)[2] <- "S_socialBees_Sainfoin"
df37[order(df37$S_socialBees_Sainfoin),] 

#Solitary bees (without honeybees)
df38 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "solitary_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Sainfoin"& nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df38)[2] <- "S_solitaryBees_Sainfoin"
df38[order(df38$S_solitaryBees_Sainfoin),] 

#Other Aculeata / no !
#df39 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Wasps" & phytometer_plant == "Sainfoin"& nectar_robber != 1, FUN=function(x)length(unique(x)))
#colnames(df39)[2] <- "S_otherAculeata_Sainfoin"
#df39[order(df39$S_otherAculeata_Sainfoin),]

#Hoverflies
#df40 <- aggregate(taxon ~ Id, data=df2, subset = family == "Syrphidae" & phytometer_plant == "Sainfoin"& nectar_robber != 1, FUN=function(x)length(unique(x)))
#colnames(df40)[2] <- "S_Syrphidae_Sainfoin"
#df40[order(df40$S_Syrphidae_Sainfoin),]

#Beetles / no !
#df41 <- aggregate(taxon ~ Id, data=df2, subset = order == "Coleoptera" & phytometer_plant == "Sainfoin"& nectar_robber != 1, FUN=function(x)length(unique(x)))
#colnames(df41)[2] <- "S_Coleoptera_Sainfoin"
#df41[order(df41$S_Coleoptera_Sainfoin),]

#All pollinators
df42 <- aggregate(taxon ~ Id, data=df2, subset = phytometer_plant == "Sainfoin"& nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df42)[2] <- "S_allPollinators_Sainfoin"
df42[order(df42$S_allPollinators_Sainfoin),]

########################
###### Comfrey ########
########################

names(df2)

#All bees (Anthophila)
df43 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & phytometer_plant == "Comfrey" & nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df43)[2] <- "S_allBees_Comfrey"
df43[order(df43$S_allBees_Comfrey),] 

#All bees (Anthophila) without honeybees
df44 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & taxon != "Apis mellifera" & phytometer_plant == "Comfrey" & nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df44)[2] <- "S_allBees_noApis_Comfrey"
df44[order(df44$S_allBees_noApis_Comfrey),] 

#Honeybees / no !
#df45 <- aggregate(taxon ~ Id, data=df2, subset = taxon == "Apis mellifera" & phytometer_plant == "Comfrey" & nectar_robber != 1, FUN=function(x)length(unique(x)))
#colnames(df45)[2] <- "S_Apis_Comfrey"
#df45[order(df45$S_Apis_Comfrey),] 

#Bumblebees
df46 <- aggregate(taxon ~ Id, data=df2, subset = genus == "Bombus" & phytometer_plant == "Comfrey" & nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df46)[2] <- "S_Bombus_Comfrey"
df46[order(df46$S_Bombus_Comfrey),] 

#Social bees (without honeybees)
df47 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "social_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Comfrey" & nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df47)[2] <- "S_socialBees_Comfrey"
df47[order(df47$S_socialBees_Comfrey),] 

#Solitary bees (without honeybees) / no !
#df48 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "Solitary" & taxon != "Apis mellifera" & phytometer_plant == "Comfrey" & nectar_robber != 1, FUN=function(x)length(unique(x)))
#colnames(df48)[2] <- "S_solitaryBees_Comfrey"
#df48[order(df48$S_solitaryBees_Comfrey),] 

#Other Aculeata / no !
#df49 <- aggregate(taxon ~ Id, data=df2, subset = pollinator_group == "Wasps" & phytometer_plant == "Comfrey" & nectar_robber != 1, FUN=function(x)length(unique(x)))
#colnames(df49)[2] <- "S_otherAculeata_Comfrey"
#df49[order(df49$S_otherAculeata_Comfrey),]

#Hoverflies / no !
#df50 <- aggregate(taxon ~ Id, data=df2, subset = family == "Syrphidae" & phytometer_plant == "Comfrey" & nectar_robber != 1, FUN=function(x)length(unique(x)))
#colnames(df50)[2] <- "S_Syrphidae_Comfrey"
#df50[order(df50$S_Syrphidae_Comfrey),]

#Beetles / no !
#df51 <- aggregate(taxon ~ Id, data=df2, subset = order == "Coleoptera" & phytometer_plant == "Comfrey" & nectar_robber != 1, FUN=function(x)length(unique(x)))
#colnames(df51)[2] <- "S_Coleoptera_Comfrey"
#df51[order(df51$S_Coleoptera_Comfrey),]

#All pollinators
df52 <- aggregate(taxon ~ Id, data=df2, subset = phytometer_plant == "Comfrey" & nectar_robber != 1, FUN=function(x)length(unique(x)))
colnames(df52)[2] <- "S_allPollinators_Comfrey"
df52[order(df52$S_allPollinators_Comfrey),]


########################################################################################################

### Combine the datasets
dfx <- Reduce(function(x, y) merge(x, y, all=TRUE, by="Id", sort=TRUE), list(df3,df4,df6,df7,df8,df9,df10,df11,df12,df13,df14,df16,df17,df18, df19,df20,df21,df22,df23,df24,df26,df27,df28,df29,df30,df31,df32,df33,df34,df36,df37,df38,df42,df43,df44,df46,df47,df52))
head(dfx) #With "all=TRUE" we also have the zeros, which written as NA's but should be transformed in 0's

#Define a function and apply it to the dataframe
NA_to_0 <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}

dfx <- as.data.frame(apply(dfx, MARGIN=2, FUN=NA_to_0))
head(dfx)

write.table(dfx, file = "cleaned_data/pollinator_SR_aggregated_garden_phytometer.txt", 
            append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

