#####################################################################
#####################################################################
#####################################################################
###
### 3. Calculate pollinator abundances per garden phytometer
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

traits <-read.csv("raw_data/RejiChacko_etal_2025_EnviDat/06_trait_data/individual_traits.csv", header = TRUE, sep = ",")
checklist <- read.csv("raw_data/RejiChacko_etal_2025_EnviDat/04_taxonomic_data/taxa_checklist.csv")

traits <- merge(traits, checklist, by = "taxon")
rm(checklist)
names(traits)

####Prepare the dataset: where do we have NA's?
summary(traits)
table(traits$family) # exclude empties
table(traits$genus) #exlucde empties, 3 more than famili = undetermiend syrphids
traits <- traits[traits$family!="",]
traits <- traits[traits$genus!="",]
summary(traits)

###############################
#### All plants per garden ####
###############################

#How many taxa do we have?
unique(traits$taxon)#163!
sort(table(traits$taxon), decreasing = TRUE) #List of observations per species.  

#All bees (Anthophila)
bees <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & nectar_robber != 1, sum)
colnames(bees)[2] <- "A_allBees"
bees[order(bees$A_allBees),] 

#All bees (Anthophila) without honeybees
nohoneybees <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & taxon != "Apis mellifera"& nectar_robber != 1, sum)
colnames(nohoneybees)[2] <- "A_allBees_noApis"
nohoneybees[order(nohoneybees$A_allBees_noApis),] 

#Honeybees
honeybees <- aggregate(adult ~ Id, data=traits, subset = taxon == "Apis mellifera", sum)
colnames(honeybees)[2] <- "A_Apis"
honeybees[order(honeybees$A_Apis),] 

#Bumblebees
bumble <- aggregate(adult ~ Id, data=traits, subset = genus == "Bombus" & nectar_robber != 1, sum)
colnames(bumble)[2] <- "A_Bombus"
bumble[order(bumble$A_Bombus),] 

#Social bees (without honeybees)
wildsocial <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "social_Bees" & taxon != "Apis mellifera"& nectar_robber != 1, sum)
colnames(wildsocial)[2] <- "A_socialBees"
wildsocial[order(wildsocial$A_socialBees),] 

#Solitary bees (without honeybees)
wildsolitary <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & pollinator_group_sociality == "solitary_Bees" & taxon != "Apis mellifera"& nectar_robber != 1, sum)
colnames(wildsolitary)[2] <- "A_solitaryBees"
wildsolitary[order(wildsolitary$A_solitaryBees),] 

#Other Aculeata
other_acu <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Wasps", sum)
colnames(other_acu)[2] <- "A_otherAculeata"
other_acu[order(other_acu$A_otherAculeata),]

#Hoverflies
hoverflies <- aggregate(adult ~ Id, data=traits, subset = family == "Syrphidae"& nectar_robber != 1, sum)
colnames(hoverflies)[2] <- "A_Syrphidae"
hoverflies[order(hoverflies$A_Syrphidae),]

#Beetles
beetles <- aggregate(adult ~ Id, data=traits, subset = order == "Coleoptera", sum)
colnames(beetles)[2] <- "A_Coleoptera"
beetles[order(beetles$A_Coleoptera),]

#All pollinators
all <- aggregate(adult ~ Id, data=traits, subset = nectar_robber != 1, sum)
colnames(all)[2] <- "A_allPollinators"
all[order(all$A_allPollinators),]

###################################
#### Per phytometer per garden ####
###################################

####################
#### Carrot ########
####################

#All bees (Anthophila)
df13 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & phytometer_plant == "Carrot", sum)
colnames(df13)[2] <- "A_allBees_Carrot"
df13[order(df13$A_allBees_Carrot),] 

#All bees (Anthophila) without honeybees
df14 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & taxon != "Apis mellifera" & phytometer_plant == "Carrot", sum)
colnames(df14)[2] <- "A_allBees_noApis_Carrot"
df14[order(df14$A_allBees_noApis_Carrot),] 

#Honeybees
df15 <- aggregate(adult ~ Id, data=traits, subset = taxon == "Apis mellifera" & phytometer_plant == "Carrot", sum)
colnames(df15)[2] <- "A_Apis_Carrot"
df15[order(df15$A_Apis_Carrot),] 

#Bumblebees
df16 <- aggregate(adult ~ Id, data=traits, subset = genus == "Bombus" & phytometer_plant == "Carrot", sum)
colnames(df16)[2] <- "A_Bombus_Carrot"
df16[order(df16$A_Bombus_Carrot),] 

#Social bees (without honeybees)
df17 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" &pollinator_group_sociality == "social_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Carrot", sum)
colnames(df17)[2] <- "A_socialBees_Carrot"
df17[order(df17$A_socialBees_Carrot),] 

#Solitary bees (without honeybees)
df18 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" &pollinator_group_sociality == "solitary_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Carrot", sum)
colnames(df18)[2] <- "A_solitaryBees_Carrot"
df18[order(df18$A_solitaryBees_Carrot),] 

#Other Aculeata
df19 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Wasps" & phytometer_plant == "Carrot", sum)
colnames(df19)[2] <- "A_otherAculeata_Carrot"
df19[order(df19$A_otherAculeata_Carrot),]

#Hoverflies
df20 <- aggregate(adult ~ Id, data=traits, subset = family == "Syrphidae" & phytometer_plant == "Carrot", sum)
colnames(df20)[2] <- "A_Syrphidae_Carrot"
df20[order(df20$A_Syrphidae_Carrot),]

#Beetles
df21 <- aggregate(adult ~ Id, data=traits, subset = order == "Coleoptera" & phytometer_plant == "Carrot", sum)
colnames(df21)[2] <- "A_Coleoptera_Carrot"
df21[order(df21$A_Coleoptera_Carrot),]

#All pollinators
df22 <- aggregate(adult ~ Id, data=traits, subset = phytometer_plant == "Carrot", sum)
colnames(df22)[2] <- "A_allPollinators_Carrot"
df22[order(df22$A_allPollinators_Carrot),]

######################
###### Radish ########
######################

#All bees (Anthophila)
df23 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & phytometer_plant == "Radish", sum)
colnames(df23)[2] <- "A_allBees_Radish"
df23[order(df23$A_allBees_Radish),] 

#All bees (Anthophila) without honeybees
df24 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & taxon != "Apis mellifera" & phytometer_plant == "Radish", sum)
colnames(df24)[2] <- "A_allBees_noApis_Radish"
df24[order(df24$A_allBees_noApis_Radish),] 

#Honeybees
df25 <- aggregate(adult ~ Id, data=traits, subset = taxon == "Apis mellifera" & phytometer_plant == "Radish", sum)
colnames(df25)[2] <- "A_Apis_Radish"
df25[order(df25$A_Apis_Radish),] 

#Bumblebees
df26 <- aggregate(adult ~ Id, data=traits, subset = genus == "Bombus" & phytometer_plant == "Radish", sum)
colnames(df26)[2] <- "A_Bombus_Radish"
df26[order(df26$A_Bombus_Radish),] 

#Social bees (without honeybees)
df27 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" &pollinator_group_sociality == "social_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Radish", sum)
colnames(df27)[2] <- "A_socialBees_Radish"
df27[order(df27$A_socialBees_Radish),] 

#Solitary bees (without honeybees)
df28 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" &pollinator_group_sociality == "solitary_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Radish", sum)
colnames(df28)[2] <- "A_solitaryBees_Radish"
df28[order(df28$A_solitaryBees_Radish),] 

#Other Aculeata
df29 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Wasps" & phytometer_plant == "Radish", sum)
colnames(df29)[2] <- "A_otherAculeata_Radish"
df29[order(df29$A_otherAculeata_Radish),]

#Hoverflies
df30 <- aggregate(adult ~ Id, data=traits, subset = family == "Syrphidae" & phytometer_plant == "Radish", sum)
colnames(df30)[2] <- "A_Syrphidae_Radish"
df30[order(df30$A_Syrphidae_Radish),]

#Beetles
df31 <- aggregate(adult ~ Id, data=traits, subset = order == "Coleoptera" & phytometer_plant == "Radish", sum)
colnames(df31)[2] <- "A_Coleoptera_Radish"
df31[order(df31$A_Coleoptera_Radish),]

#All pollinators
df32 <- aggregate(adult ~ Id, data=traits, subset = phytometer_plant == "Radish", sum)
colnames(df32)[2] <- "A_allPollinators_Radish"
df32[order(df32$A_allPollinators_Radish),]

########################
###### Sainfoin ########
########################

#All bees (Anthophila)
df33 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & phytometer_plant == "Sainfoin" & nectar_robber != 1, sum)
colnames(df33)[2] <- "A_allBees_Sainfoin"
df33[order(df33$A_allBees_Sainfoin),] 

#All bees (Anthophila) without honeybees
df34 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & taxon != "Apis mellifera" & phytometer_plant == "Sainfoin"& nectar_robber != 1, sum)
colnames(df34)[2] <- "A_allBees_noApis_Sainfoin"
df34[order(df34$A_allBees_noApis_Sainfoin),] 

#Honeybees
df35 <- aggregate(adult ~ Id, data=traits, subset = taxon == "Apis mellifera" & phytometer_plant == "Sainfoin"& nectar_robber != 1, sum)
colnames(df35)[2] <- "A_Apis_Sainfoin"
df35[order(df35$A_Apis_Sainfoin),] 

#Bumblebees
df36 <- aggregate(adult ~ Id, data=traits, subset = genus == "Bombus" & phytometer_plant == "Sainfoin"& nectar_robber != 1, sum)
colnames(df36)[2] <- "A_Bombus_Sainfoin"
df36[order(df36$A_Bombus_Sainfoin),] 

#Social bees (without honeybees)
df37 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" &pollinator_group_sociality == "social_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Sainfoin"& nectar_robber != 1, sum)
colnames(df37)[2] <- "A_socialBees_Sainfoin"
df37[order(df37$A_socialBees_Sainfoin),] 

#Solitary bees (without honeybees)
df38 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" &pollinator_group_sociality == "solitary_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Sainfoin"& nectar_robber != 1, sum)
colnames(df38)[2] <- "A_solitaryBees_Sainfoin"
df38[order(df38$A_solitaryBees_Sainfoin),] 

#Other Aculeata / no !
#df39 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "other_Aculeata" & phytometer_plant == "Sainfoin"& nectar_robber != 1, sum)
#colnames(df39)[2] <- "A_otherAculeata_Sainfoin"
#df39[order(df39$A_otherAculeata_Sainfoin),]

#Hoverflies
#df40 <- aggregate(adult ~ Id, data=traits, subset = family == "Syrphidae" & phytometer_plant == "Sainfoin"& nectar_robber != 1, sum)
#colnames(df40)[2] <- "A_Syrphidae_Sainfoin"
#df40[order(df40$A_Syrphidae_Sainfoin),]

#Beetles / no !
#df41 <- aggregate(adult ~ Id, data=traits, subset = order == "Coleoptera" & phytometer_plant == "Sainfoin"& nectar_robber != 1, sum)
#colnames(df41)[2] <- "A_Coleoptera_Sainfoin"
#df41[order(df41$A_Coleoptera_Sainfoin),]

#All pollinators
df42 <- aggregate(adult ~ Id, data=traits, subset = phytometer_plant == "Sainfoin"& nectar_robber != 1, sum)
colnames(df42)[2] <- "A_allPollinators_Sainfoin"
df42[order(df42$A_allPollinators_Sainfoin),]

########################
###### Comfrey ########
########################

names(traits)

#All bees (Anthophila)
df43 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & phytometer_plant == "Comfrey" & nectar_robber != 1, sum)
colnames(df43)[2] <- "A_allBees_Comfrey"
df43[order(df43$A_allBees_Comfrey),] 

#All bees (Anthophila) without honeybees
df44 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" & taxon != "Apis mellifera" & phytometer_plant == "Comfrey" & nectar_robber != 1, sum)
colnames(df44)[2] <- "A_allBees_noApis_Comfrey"
df44[order(df44$A_allBees_noApis_Comfrey),] 

#Honeybees / no !
#df45 <- aggregate(adult ~ Id, data=traits, subset = taxon == "Apis mellifera" & phytometer_plant == "Comfrey" & nectar_robber != 1, sum)
#colnames(df45)[2] <- "A_Apis_Comfrey"
#df45[order(df45$A_Apis_Comfrey),] 

#Bumblebees
df46 <- aggregate(adult ~ Id, data=traits, subset = genus == "Bombus" & phytometer_plant == "Comfrey" & nectar_robber != 1, sum)
colnames(df46)[2] <- "A_Bombus_Comfrey"
df46[order(df46$A_Bombus_Comfrey),] 

#Social bees (without honeybees)
df47 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" &pollinator_group_sociality == "social_Bees" & taxon != "Apis mellifera" & phytometer_plant == "Comfrey" & nectar_robber != 1, sum)
colnames(df47)[2] <- "A_socialBees_Comfrey"
df47[order(df47$A_socialBees_Comfrey),] 

#Solitary bees (without honeybees) / no !
#df48 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "Anthophila" &pollinator_group_sociality == "Solitary" & taxon != "Apis mellifera" & phytometer_plant == "Comfrey" & nectar_robber != 1, sum)
#colnames(df48)[2] <- "A_solitaryBees_Comfrey"
#df48[order(df48$A_solitaryBees_Comfrey),] 

#Other Aculeata / no !
#df49 <- aggregate(adult ~ Id, data=traits, subset = pollinator_group == "other_Aculeata" & phytometer_plant == "Comfrey" & nectar_robber != 1, sum)
#colnames(df49)[2] <- "A_otherAculeata_Comfrey"
#df49[order(df49$A_otherAculeata_Comfrey),]

#Hoverflies / no !
#df50 <- aggregate(adult ~ Id, data=traits, subset = family == "Syrphidae" & phytometer_plant == "Comfrey" & nectar_robber != 1, sum)
#colnames(df50)[2] <- "A_Syrphidae_Comfrey"
#df50[order(df50$A_Syrphidae_Comfrey),]

#Beetles / no !
#df51 <- aggregate(adult ~ Id, data=traits, subset = order == "Coleoptera" & phytometer_plant == "Comfrey" & nectar_robber != 1, sum)
#colnames(df51)[2] <- "A_Coleoptera_Comfrey"
#df51[order(df51$A_Coleoptera_Comfrey),]

#All pollinators
df52 <- aggregate(adult ~ Id, data=traits, subset = phytometer_plant == "Comfrey" & nectar_robber != 1, sum)
colnames(df52)[2] <- "A_allPollinators_Comfrey"
df52[order(df52$A_allPollinators_Comfrey),]


########################################################################################################

### Combine the datasets
dfx <- Reduce(function(x, y) merge(x, y, all=TRUE, by="Id", sort=TRUE), list(bees,nohoneybees,honeybees,bumble,wildsocial,wildsolitary,other_acu,hoverflies,beetles,all,df13,df14,df15,df16,df17,df18, df19,df20,df21,df22,df23,df24,df25,df26,df27,df28,df29,df30,df31,df32,df33,df34,df35,df36,df37,df38,df42,df43,df44,df46,df47,df52))
head(dfx) #With "all=TRUE" we also have the zeros, which written as NA's but should be transformed in 0's

#Define a function and apply it to the dataframe
NA_to_0 <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}

dfx <- as.data.frame(apply(dfx, MARGIN=2, FUN=NA_to_0))
head(dfx)

write.table(dfx, file = "cleaned_data/pollinator_abundance_aggregated_garden_phytometer.txt", 
            append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


##################################################################################################################################
##################################################################################################################################

