#####################################################################
#####################################################################
#####################################################################
###
### 1.Summarise fieldwork data
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

raw <-readxl::read_xlsx("raw_data/RejiChacko_etal_2025_EnviDat/05_field_data/raw_sampling_data.xlsx")

head(raw)
unique(raw$capture_period) #22 sampling dates
str(raw)

###########################################################################
### 1.1 Calculating the number of flowers surveyed per date and garden  ###
###########################################################################

#Extract the number of flowers for each observation round (by phytometer species + garden + date). 
#Unique because there is only one unique value (one estimate of N flowers) for each sampling day, but the number is written for every time window. 

flower_agg <- aggregate(cbind(n_flowers_carrot, n_flowers_radish, n_flowers_sainfoin, n_flowers_comfrey) ~ Id + capture_period, data=raw, FUN=function(x) unique(x))

dim(flower_agg) 
flower_agg[order(flower_agg$Id),]
head(flower_agg)

#Check briefly the distribution
par(mfrow=c(2,2))
hist(flower_agg$n_flowers_radish)
hist(flower_agg$n_flowers_sainfoin)
hist(flower_agg$n_flowers_carrot)
hist(flower_agg$n_flowers_comfrey)

boxplot(flower_agg$n_flowers_radish)
boxplot(flower_agg$n_flowers_sainfoin)
boxplot(flower_agg$n_flowers_carrot)
boxplot(flower_agg$n_flowers_comfrey)#There are gardens with very rich flower sets in all species. 
par(mfrow=c(1,1))

#######################################################################
### 1.2 Calculating the number fieldwork hours per date and garden  ###
#######################################################################

#Sampling effort in minutes fieldwork for each plant species separately (take into account the number of operators: for example, 2 persons observing for two h = 2 h of sampling effort) 
#Condition: phytometer plant must be in flower (>0 flowers), otherwise no observation and thus sampling effort = 0

effort_per_date <- aggregate(total_sampling_effort_min ~ Id + capture_period , data=raw, FUN=function(x) sum(x))
effort_per_date[order(effort_per_date$Id),]
effort_per_date

#Join with garden and capture_period
effort_flowers <- merge(flower_agg, effort_per_date, by=c("Id","capture_period"))

dim(effort_flowers)
effort_flowers[order(effort_flowers$total_sampling_effort_min),]
head(effort_flowers)

###############################################################################################
### 1.3 Calculating the number of fieldwork hours per date and garden and phytometer plant  ###
###############################################################################################


#Sum sampling effort of each plant where N flowers >0

carrot <- aggregate(total_sampling_effort_min ~ Id + capture_period, data=effort_flowers, subset = n_flowers_carrot>0, FUN=function(x) sum(x))
raphanus <- aggregate(total_sampling_effort_min ~ Id + capture_period, data=effort_flowers, subset = n_flowers_radish>0, FUN=function(x) sum(x))
onobrychis <- aggregate(total_sampling_effort_min ~ Id + capture_period, data=effort_flowers, subset = n_flowers_sainfoin>0, FUN=function(x) sum(x))
symphytum <- aggregate(total_sampling_effort_min ~ Id + capture_period, data=effort_flowers, subset = n_flowers_comfrey>0, FUN=function(x) sum(x))

dim(carrot)
dim(raphanus)
dim(onobrychis)
dim(symphytum)

colnames(carrot)[3]<-"sampling_effort_min_carrot"
colnames(raphanus)[3]<-"sampling_effort_min_radish"
colnames(onobrychis)[3]<-"sampling_effort_min_sainfoin"
colnames(symphytum)[3]<-"sampling_effort_min_comfrey"

####################################################
## Combine all the (sampling-date-wise) data sets ##
####################################################

effort_all <- Reduce(function(x, y) merge(x, y, all=TRUE, by=c("Id","capture_period")), list(effort_flowers, carrot, raphanus, onobrychis, symphytum))

str(effort_all) 
summary(effort_all)
dim(effort_all) #91 observation rounds!(non-empty ones)

#Replace NAs with 0
library(tidyr)
effort_all <- effort_all %>% replace_na(list(sampling_effort_min_carrot = 0, sampling_effort_min_radish = 0, sampling_effort_min_sainfoin = 0, sampling_effort_min_comfrey = 0))

summary(effort_all)
str(effort_all)

#Write 
write.table(effort_all, file = "cleaned_data/sampling_effort_gardens_plants_aggregated_per_date.txt", append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

##############
#GARDEN SCALE#
##############

########################################################################
### 1.4 Calculating the total number of flowers surveyed per garden  ###
########################################################################

#The total number of flowers surveyed per garden and phytometer plant 
#We use flower_agg, where the number of flower is given at its highest resolution (date-garden)

effort_all <- aggregate(cbind(n_flowers_radish, n_flowers_sainfoin, n_flowers_carrot, n_flowers_comfrey) ~ Id, data=flower_agg, FUN=function(x) sum(x))
effort_all

#######################################################################
### 1.5 Calculating the total number of fieldwork hours per garden  ###
#######################################################################

#Total sampling effort in min fieldwork (taking into account all persons)
#We use the rawdata of raw since the resolution is time-window-wise:

effort_garden <- aggregate(total_sampling_effort_min ~ Id, data=raw, FUN=function(x) sum(x))
head(effort_garden)

#Plot the data
hist(effort_garden$total_sampling_effort_min)
plot(effort_garden$total_sampling_effort_min~as.factor(effort_garden$Id)) #Garden 39 should be omitted

#The total min fieldworks (with and without garden 39)
sum(effort_garden$total_sampling_effort_min)/60 #1243.833
sum(effort_garden[which(effort_garden$Id != 39),]$total_sampling_effort_min)/60 #1235.333

##########################################################################################################
### 1.6 Calculating the total number of fieldwork hours per garden and phytometer plant (if in flower) ###
##########################################################################################################

#We use effort_flowers where the number of flowers + sampling hours is aggregated already at the date-garden scale:
carrot <- aggregate(total_sampling_effort_min ~ Id, data=effort_flowers, subset = n_flowers_carrot>0, FUN=function(x) sum(x))
raphanus <- aggregate(total_sampling_effort_min ~ Id, data=effort_flowers, subset = n_flowers_radish>0, FUN=function(x) sum(x))
onobrychis <- aggregate(total_sampling_effort_min ~ Id, data=effort_flowers, subset = n_flowers_sainfoin>0, FUN=function(x) sum(x))
symphytum <- aggregate(total_sampling_effort_min ~ Id, data=effort_flowers, subset = n_flowers_comfrey>0, FUN=function(x) sum(x))

dim(carrot)
dim(raphanus)
dim(onobrychis)
dim(symphytum)

colnames(carrot)[2]<-"sampling_effort_min_carrot"
colnames(raphanus)[2]<-"sampling_effort_min_radish"
colnames(onobrychis)[2]<-"sampling_effort_min_sainfoin"
colnames(symphytum)[2]<-"sampling_effort_min_comfrey"


#############################################
## Combine all the (garden-wise) data sets ##
#############################################

all <- Reduce(function(x, y) merge(x, y, all=TRUE, by="Id", sort=TRUE), list(effort_all, effort_garden, carrot, raphanus, onobrychis, symphytum))
all 

write.table(all, file = "cleaned_data/sampling_effort_gardens_plants_aggregated.txt", append = FALSE, quote = TRUE, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
