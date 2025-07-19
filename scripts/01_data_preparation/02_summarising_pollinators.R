#####################################################################
#####################################################################
#####################################################################
###
### 2. Summary stats for pollinators
### Code by: David Frey and Merin Reji Chacko
### contact: merin.rejichacko@gmail.com
### Last edited: 18/07/2025
### 
#####################################################################
#####################################################################
#####################################################################
###
### Prep

rm(list=ls())

#Load the pollinator data
trait <-read.csv("raw_data/RejiChacko_etal_2025_EnviDat/06_trait_data/individual_traits.csv", header = TRUE, sep = ",")
head(trait)
dim(trait)
names(trait)

#Exclude lost samples and flies other than Hoverflies ("not a hoverfly", N = 25) and get rid of empty factor levels
trait <- droplevels(subset(trait, comment != "not a hoverfly" & comment != "sample lost")) 
dim(trait)

head(trait)
names(trait)

# we need pollinator_group data 

checklist <- read.csv("raw_data/RejiChacko_etal_2025_EnviDat/04_taxonomic_data/taxa_checklist.csv")
table(trait$taxon %in% checklist$taxon)
trait <- merge(trait, checklist, by = "taxon")
rm(checklist)

#Make summary statistics 
library(dplyr)
#Check: https://www3.nd.edu/~steve/computing_with_data/24_dplyr/dplyr.html

#Nectar thieves
thieves <- trait  %>% filter(nectar_robber==1) %>% group_by(pollinator_group, phytometer_plant) %>% 
  select(pollinator_group:phytometer_plant, adult, taxon) %>% 
  summarise(abundance = sum(adult), S = n_distinct(taxon)) %>% 
  as.data.frame() %>% setNames(c("pollinator_group","phytometer_plant","abundance","species_richness"))
thieves

sum(thieves$abundance)
sum(thieves$abundance)/dim(trait)[1]*100 #%nectar thieves of entire dataset

sum(thieves[which(thieves$phytometer_plant=="Comfrey"),]$abundance)
sum(thieves[which(thieves$phytometer_plant=="Sainfoin"),]$abundance)

#All bees (Anthophila)
allbees <- trait %>% filter(nectar_robber!=1 & pollinator_group=="Anthophila") %>% 
  group_by(pollinator_group, phytometer_plant)  %>% 
  select(pollinator_group:phytometer_plant, adult, taxon) %>% 
  summarise(abundance = sum(adult), S = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","phytometer_plant","abundance","species_richness"))

allbees

sum(allbees$abundance)/dim(trait[which(trait$nectar_robber!=1),])[1]*100 #%bees

#All bees (Anthophila) no Honeybees
nohoneybees <- trait %>%  
  filter(nectar_robber!=1 & pollinator_group=="Anthophila" & taxon!="Apis mellifera") %>% 
  group_by(pollinator_group, phytometer_plant)  %>% 
  select(pollinator_group:phytometer_plant, adult, taxon) %>% 
  summarise(abundance = sum(adult), S = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","phytometer_plant","abundance","species_richness"))

nohoneybees$pollinator_group <- rep("Anthophila_noApis",4)
nohoneybees

#Honeybees
honeybees <- trait %>%  
  filter(nectar_robber!=1 & taxon=="Apis mellifera") %>% 
  group_by(taxon, phytometer_plant)  %>% 
  select(taxon:phytometer_plant, adult, taxon) %>%
  summarise(abundance = sum(adult), S = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","phytometer_plant","abundance","species_richness"))

honeybees

sum(honeybees$abundance)/dim(trait[which(trait$nectar_robber!=1),])[1]*100 #%bees

#Bumblebees
bumblebees <- trait %>%  
  filter(nectar_robber!=1 & genus=="Bombus") %>% 
  group_by(genus, phytometer_plant)  %>% 
  select(genus:phytometer_plant, adult, taxon) %>% 
  summarise(abundance = sum(adult), S = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","phytometer_plant","abundance","species_richness"))

bumblebees

#pollinator_group_sociality
social <- trait %>% 
  filter(nectar_robber!=1 & pollinator_group=="Anthophila" & taxon!="Apis mellifera") %>% 
  group_by(pollinator_group_sociality, phytometer_plant) %>% 
  select(pollinator_group_sociality:phytometer_plant, adult, taxon) %>% 
  summarise(abundance = sum(adult), S = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","phytometer_plant","abundance","species_richness"))

social

# note that "no" means polymorphic, which we will drop in some analyses

# wasps
other_aculeata <- trait %>% 
  filter(nectar_robber!=1 & pollinator_group=="Wasps") %>% 
  group_by(pollinator_group, phytometer_plant) %>% 
  select(pollinator_group:phytometer_plant, adult, taxon) %>% 
  summarise(abundance = sum(adult), S = n_distinct(taxon)) %>% 
  as.data.frame()%>% 
  setNames(c("pollinator_group","phytometer_plant","abundance","species_richness"))

other_aculeata

#Syrphids
hoverflies <- trait %>% 
  filter(nectar_robber!=1 & family=="Syrphidae") %>% 
  group_by(family, phytometer_plant) %>% 
  select(family:phytometer_plant, adult, taxon) %>% 
  summarise(abundance = sum(adult), S = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","phytometer_plant","abundance","species_richness"))

hoverflies

sum(hoverflies$abundance)/dim(trait[which(trait$nectar_robber!=1),])[1]*100 #%Syrphids

#Beetles
beetles <- trait %>% 
  filter(nectar_robber!=1 & order=="Coleoptera") %>% 
  group_by(order, phytometer_plant) %>% 
  select(order:phytometer_plant, adult, taxon) %>% 
  summarise(abundance = sum(adult), S = n_distinct(taxon)) %>% 
  as.data.frame()%>% 
  setNames(c("pollinator_group","phytometer_plant","abundance","species_richness"))

beetles

#All
all <- trait %>% 
  filter(nectar_robber!=1) %>% 
  group_by(phytometer_plant) %>% 
  select(phytometer_plant, adult, taxon) %>% 
  summarise(abundance = sum(adult), S = n_distinct(taxon)) %>% 
  as.data.frame()%>% 
  setNames(c("phytometer_plant","abundance","species_richness"))

all
all <- mutate(all, pollinator_group=rep("All_pollinators",4))
all
all <- all[,c(4,1,2,3)]
all

#Combine them:
combined <- rbind(allbees,nohoneybees,honeybees,bumblebees,social,other_aculeata,hoverflies,beetles,all)
combined
combined <- subset(combined, pollinator_group != "no") # these are polymorphic
combined

str(combined)

#Split them into phytometer species
dat_carrot <- combined[which(combined$phytometer_plant=="Carrot"),]
dat_carrot <- dat_carrot[,c(-2)]
dat_carrot

dat_radish <- combined[which(combined$phytometer_plant=="Radish"),]
dat_radish <- dat_radish[,c(-2)]
dat_radish

dat_sainfoin<- combined[which(combined$phytometer_plant=="Sainfoin"),]
dat_sainfoin <- dat_sainfoin[,c(-2)]
dat_sainfoin

dat_comfrey<- combined[which(combined$phytometer_plant=="Comfrey"),]
dat_comfrey <- dat_comfrey[,c(-2)]
dat_comfrey

#Merge them in a transposed table
dat_combined <- Reduce(function(x, y) merge(x, y, all=TRUE, by="pollinator_group"), list(dat_carrot,dat_radish,dat_sainfoin,dat_comfrey))
# THE ERROR IS NORMAL! WE JUST NEED TO RENAME EVERYTHING 

str(dat_combined)
dat_combined

colnames(dat_combined)[2]<-"A_carrot"
colnames(dat_combined)[3]<-"S_carrot"
colnames(dat_combined)[4]<-"A_radish"
colnames(dat_combined)[5]<-"S_radish"
colnames(dat_combined)[6]<-"A_sainfoin"
colnames(dat_combined)[7]<-"S_sainfoin"
colnames(dat_combined)[8]<-"A_comfrey"
colnames(dat_combined)[9]<-"S_comfrey"

dat_combined

#Replace NAs with 0
library(tidyr)
dat_combined <- dat_combined %>% replace_na(list(A_sainfoin = 0, S_sainfoin = 0, A_comfrey = 0, S_comfrey = 0))

dat_combined <- mutate(dat_combined, A_all = A_carrot+A_radish+A_sainfoin+A_comfrey)
dat_combined

#Overall species richness
#All bees (Anthophila)
allbees <- trait %>% 
  filter(nectar_robber!=1 & pollinator_group=="Anthophila") %>% 
  group_by(pollinator_group)  %>% 
  select(pollinator_group, taxon) %>% 
  summarise(S_all = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","S_all"))

honeybees <- trait %>% 
  filter(nectar_robber!=1 & taxon=="Apis mellifera") %>% 
  group_by(taxon)  %>% 
  select(taxon, taxon) %>% 
  summarise(S_all = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","S_all"))

nohoneybees <- trait %>% 
  filter(nectar_robber!=1 & pollinator_group=="Anthophila" & taxon!="Apis mellifera") %>% 
  group_by(pollinator_group)  %>% 
  select(pollinator_group, taxon) %>% 
  summarise(S_all = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","S_all"))

bumblebees <- trait %>% 
  filter(nectar_robber!=1 & genus=="Bombus")  %>% 
  group_by(genus)  %>% 
  select(genus, taxon) %>% 
  summarise(S_all = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","S_all"))

social_wild <- trait %>% 
  filter(nectar_robber!=1 & pollinator_group=="Anthophila" & taxon!="Apis mellifera"& pollinator_group_sociality == "social_Bees") %>% 
  group_by(pollinator_group_sociality) %>% 
  select(pollinator_group_sociality, taxon) %>% 
  summarise(S_all = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","S_all"))

solitary_wild <- trait %>% 
  filter(nectar_robber!=1 & pollinator_group=="Anthophila" & taxon!="Apis mellifera"& pollinator_group_sociality == "solitary_Bees") %>% 
  group_by(pollinator_group_sociality) %>% 
  select(pollinator_group_sociality, taxon) %>% 
  summarise(S_all = n_distinct(taxon)) %>% 
  as.data.frame() %>% 
  setNames(c("pollinator_group","S_all"))

other_aculeata <- trait %>% 
  filter(nectar_robber!=1 & pollinator_group=="Wasps") %>% 
  group_by(pollinator_group) %>% 
  select(pollinator_group, taxon) %>% 
  summarise(S_all = n_distinct(taxon)) %>% 
  as.data.frame()%>% 
  setNames(c("pollinator_group","S_all"))

hoverflies <- trait %>% 
  filter(nectar_robber!=1 & family=="Syrphidae") %>% 
  group_by(family) %>% 
  select(family,taxon) %>% 
  summarise(S_all = n_distinct(taxon)) %>% 
  as.data.frame()%>% 
  setNames(c("pollinator_group","S_all"))

beetles <- trait %>% 
  filter(nectar_robber!=1 & order=="Coleoptera") %>% 
  group_by(order) %>% 
  select(order,taxon) %>% 
  summarise(S_all = n_distinct(taxon)) %>% 
  as.data.frame()%>% 
  setNames(c("pollinator_group","S_all"))

#All together
together <- trait %>% 
  filter(nectar_robber!=1) %>% 
  select(taxon) %>% 
  summarise(S_all = n_distinct(taxon)) %>% 
  as.data.frame() %>% setNames(c("S_all"))

allbees
honeybees
nohoneybees
nohoneybees$pollinator_group<-c("Anthophila_noApis")
nohoneybees
bumblebees
social_wild
solitary_wild
other_aculeata
hoverflies
beetles

together
together$pollinator_group <- c("All_pollinators")
together <- together[,c(2,1)]
together

#Bind them together
full <- rbind(allbees, honeybees, nohoneybees, bumblebees, social_wild, solitary_wild, other_aculeata, hoverflies, beetles, together)
full

#Merge abundance with species richness data of all phytometer plants
dat_combined <- Reduce(function(x, y) merge(x, y, all=TRUE, by="pollinator_group"), list( full,dat_combined))
dat_combined

#Reorder rows & columns
dat_combined <- dat_combined[c(2,3,4,5,8,9,7,10,6,1),c(1,3:11,2)]
dat_combined

write.table(dat_combined, file = "cleaned_data/pollinators_sum_stat.txt", append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

###########################################################################################################################################################
###########################################################################################################################################################


