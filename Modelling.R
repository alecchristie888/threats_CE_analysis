#loading necessary libraries
library(sf)
library(ggplot2)
library(dplyr)
library(tmap)
#install.packages("~/maptools_1.1-8/maptools", repos = NULL, type="source")
library(maptools)
library(ggrepel)
library(readxl)
library(glm2)
library(emmeans)
library(car)
library(DHARMa)
library(glmmTMB)
library(MASS)
library(MuMIn)
library(arm)

#set working directory
#setwd("")

#reading in files from working directory
studies_data <- read.csv("Data/CE_data_complete.csv")
non_english <- read_excel("Data/CE_non_english_data_clean.xlsx")
threats_data <- st_read("Data/Grid_mammal_amph_threat_predictions.shp")
EPI <- read_excel('Data/EPI.xlsx')
GDP <- read_excel('Data/GDP.xlsx')
lit_rates <- read_excel('Data/LiteracyRates.xlsx')
NBI <- read_excel('Data/NBI.xlsx')
WGI <- read_excel('Data/WGI.xlsx')

#cleaning up
studies_data <- studies_data[!is.na(studies_data$long),]
studies_data <- filter(studies_data, long >= -180 & long <= 180)
studies_data <- filter(studies_data, lat >= -90 & lat <= 90)
studies_data <- filter(studies_data, !(lat == 0 & long == 0))

#converting studies data to a shapefile
studies_data <- st_as_sf(studies_data, coords = c("long", "lat"), crs = 4326)

#setting studies data to same crs as threats data 
studies_data <- st_transform(studies_data, crs = st_crs(threats_data))


#################################### AMPHIBIANS ########################
#################################### AMPHIBIANS ########################


#### AMBHIBIANS AND AGRICULTURE - A2 ####

## getting threats specific to amphibians and agriculture 
A2_threats <- threats_data[!is.na(threats_data$A2),]

# applying quality filter to threats data 
A2_threats <- A2_threats %>% filter(ASR > 10)

# calculating risk
A2_threats$A2_risk <- A2_threats$A2 * A2_threats$ASR

## getting studies specific to amphibians and agriculture 
A_studies <- studies_data[studies_data$taxa == "Amphibians",]
A2_studies <- A_studies[A_studies$threattype == "Agriculture & aquaculture",]

## similarly for non-english studies 
A2_non_english <- filter(non_english, taxa == "Amphibians" & threattype == "Agriculture & aquaculture")


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
A2_studies <- subset(A2_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
A2_studies <- st_set_geometry(A2_studies, NULL)
A2_studies <- A2_studies[!is.na(A2_studies$country),]
A2_studies <- distinct(A2_studies)
A2_counts <- count(A2_studies, country)
#A2_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
A2_non_english_counts <- count(A2_non_english, country)

#adding the two to get a complete count of studies per country 
A2_counts <- merge(A2_counts, A2_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
A2_counts[is.na(A2_counts)] <- 0
A2_counts$n <- A2_counts$n.x + A2_counts$n.y
A2_counts <- subset(A2_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
A2_centroids <- st_centroid(A2_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean A2 threat value for each country
A2_centroids_sp <- as(A2_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
A2_country <- over(world_sp, A2_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## creating a complete dataframe (this contains all countries for which there is threat data)
A2_risk_country <- A2_country$A2_risk
A2_threat_country <- A2_country$A2
A2_world_threats <- mutate(world, threats = A2_threat_country, risk = A2_risk_country)
A2_world_threats <- A2_world_threats[order(A2_world_threats$NAME),]#ordering alphabetically according to country name

A2_counts <- rename(A2_counts, NAME = country)
A2_world_studies <- merge(world, A2_counts, by = "NAME", all = TRUE)

A2_world_complete <- mutate(A2_world_threats, nstudies = A2_world_studies$n)
A2_world_complete$nstudies[is.na(A2_world_complete$nstudies)] <- 0
A2_world_complete <- A2_world_complete[!is.na(A2_world_complete$threats),]
A2_world_complete <- rename(A2_world_complete, Country = NAME)
A2_world_complete <- subset(A2_world_complete, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
A2_world_complete <- st_set_geometry(A2_world_complete, NULL)

#adding socioeconomic variables
A2_LM_data <- merge(A2_world_complete, EPI, by = "Country", all.x = TRUE)
A2_LM_data <- merge(A2_LM_data, GDP, by = "Country", all.x = TRUE)
A2_LM_data <- merge(A2_LM_data, NBI, by = "Country", all.x = TRUE)
A2_LM_data <- merge(A2_LM_data, WGI, by = "Country", all.x = TRUE)
A2_LM_data <- merge(A2_LM_data, lit_rates, by = "Country", all.x = TRUE)
A2_LM_data <- subset(A2_LM_data, select = -c(dataYear, pop2022))
A2_LM_data$WGI_GE_Estimate <- as.numeric(A2_LM_data$WGI_GE_Estimate)
A2_LM_data$WGI_GE_Rank <- as.numeric(A2_LM_data$WGI_GE_Rank)
A2_LM_data$literacy_rate <- as.numeric(A2_LM_data$literacy_rate)

#### AMPHIBIANS AND POLLUTION - A9 ####

## getting threats specific to amphibians and pollution
A9_threats <- threats_data[!is.na(threats_data$A9),]

# applying quality filter to threats data 
A9_threats <- A9_threats %>% filter(ASR > 10)

# calculating risk
A9_threats$A9_risk <- A9_threats$A9 * A9_threats$ASR

## getting studies specific to amphibians and pollution
A_studies <- studies_data[studies_data$taxa == "Amphibians",]
A9_studies <- A_studies[A_studies$threattype == "Pollution",]

## similarly for non english studies
A9_non_english <- filter(non_english, taxa == "Amphibians" & threattype == "Pollution")


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
A9_studies <- subset(A9_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
A9_studies <- st_set_geometry(A9_studies, NULL)
A9_studies <- A9_studies[!is.na(A9_studies$country),]
A9_studies <- distinct(A9_studies)
A9_counts <- count(A9_studies, country)
#A9_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
A9_non_english_counts <- count(A9_non_english, country)

#adding the two to get a complete count of studies per country 
A9_counts <- merge(A9_counts, A9_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
A9_counts[is.na(A9_counts)] <- 0
A9_counts$n <- A9_counts$n.x + A9_counts$n.y
A9_counts <- subset(A9_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
A9_centroids <- st_centroid(A9_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean A9 threat value for each country
A9_centroids_sp <- as(A9_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
A9_country <- over(world_sp, A9_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## creating a complete dataframe (this contains all countries for which there is threat data)
A9_threat_country <- A9_country$A9
A9_risk_country <- A9_country$A9_risk
A9_world_threats <- mutate(world, threats = A9_threat_country, risk = A9_risk_country)
A9_world_threats <- A9_world_threats[order(A9_world_threats$NAME),]#ordering alphabetically according to country name

A9_counts <- rename(A9_counts, NAME = country)
A9_world_studies <- merge(world, A9_counts, by = "NAME", all = TRUE)

A9_world_complete <- mutate(A9_world_threats, nstudies = A9_world_studies$n)
A9_world_complete$nstudies[is.na(A9_world_complete$nstudies)] <- 0
A9_world_complete <- A9_world_complete[!is.na(A9_world_complete$threats),]
A9_world_complete <- rename(A9_world_complete, Country = NAME)
A9_world_complete <- subset(A9_world_complete, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
A9_world_complete <- st_set_geometry(A9_world_complete, NULL)

#adding socioeconomic variables
A9_LM_data <- merge(A9_world_complete, EPI, by = "Country", all.x = TRUE)
A9_LM_data <- merge(A9_LM_data, GDP, by = "Country", all.x = TRUE)
A9_LM_data <- merge(A9_LM_data, NBI, by = "Country", all.x = TRUE)
A9_LM_data <- merge(A9_LM_data, WGI, by = "Country", all.x = TRUE)
A9_LM_data <- merge(A9_LM_data, lit_rates, by = "Country", all.x = TRUE)
A9_LM_data <- subset(A9_LM_data, select = -c(dataYear, pop2022))
A9_LM_data$WGI_GE_Estimate <- as.numeric(A9_LM_data$WGI_GE_Estimate)
A9_LM_data$WGI_GE_Rank <- as.numeric(A9_LM_data$WGI_GE_Rank)
A9_LM_data$literacy_rate <- as.numeric(A9_LM_data$literacy_rate)

##################################################################################
################################### amphibian modelling ##########################
##################################################################################

A9_LM_data$Pollution <- A9_LM_data$risk
A2_LM_data$Agriculture <- A2_LM_data$risk

A_LM_data <- cbind(A9_LM_data,Agriculture=A2_LM_data$Agriculture)
A_LM_data <- na.omit(A_LM_data)


### scale variables
scale_0.5SD <- function(x){
  (x - mean(x)) / (2*sd(x))
}

A_LM_data$Agriculturescaled <- scale_0.5SD(A_LM_data$Agriculture)
A_LM_data$Pollutionscaled <- scale_0.5SD(A_LM_data$Pollution)
A_LM_data$riskscaled <- scale_0.5SD(A_LM_data$risk)
A_LM_data$EPI_newscaled <- scale_0.5SD(A_LM_data$EPI_new)
A_LM_data$GDP_2020scaled <- scale_0.5SD(A_LM_data$GDP_2020)
A_LM_data$NBIscaled <- scale_0.5SD(A_LM_data$NBI)
A_LM_data$WGI_GE_Rankscaled <- scale_0.5SD(A_LM_data$WGI_GE_Rank)
A_LM_data$literacy_ratescaled <- scale_0.5SD(A_LM_data$literacy_rate)

### poisson - results show tax model passes DHARMa tests - so can use poisson distribution, zero inflation not an issue.
glm_A <- glm(nstudies ~ Agriculturescaled + Pollutionscaled + EPI_newscaled + GDP_2020scaled + NBIscaled + WGI_GE_Rankscaled + literacy_ratescaled, family = "poisson", data = A_LM_data)
summary(glm_A)

head(A_LM_data)

testDispersion(glm_A)
simulationOutput <- simulateResiduals(fittedModel = glm_A, plot = F)
residuals(simulationOutput)
plot(simulationOutput)

###### create global model and dredge to find best model
options(na.action = "na.fail") #Must run this code once to use dredge
all_models_a <- dredge(glm_A)
#write.csv(data.frame(all_models_a),"all_models_a.csv")

all_models_a
nrow(all_models_a)

#several models with deltaAICc <2
m.top.models.2aic <- get.models(all_models_a, subset = delta <2)
length(m.top.models.2aic)
m.top.models.4aic <- get.models(all_models_a, subset = delta <4)
length(m.top.models.4aic)
m.top.models.95ci <- get.models(all_models_a, cumsum(weight) <= 0.95)
length(m.top.models.95ci)

all_models_a_avg <- summary(model.avg(m.top.models.2aic))

all_models_a_avg
sw(all_models_a_avg)
confint(all_models_a_avg,full=TRUE)

# write.csv(summary(all_models_a_avg)[9],"glm_Anb_summaryupdate.csv")
# write.csv(sw(all_models_a_avg),"glm_Anb_relimp.csv")
# write.csv(confint(all_models_a_avg,full=TRUE),"glm_Anb_confint.csv")



########################################### BIRDS ###########################################
########################################### BIRDS ###########################################

#### BIRDS AND AGRICULTURE - B2 ####

## getting threats specific to birds and agriculture
B2_threats <- threats_data[!is.na(threats_data$B2),]

# applying quality filter to threats data 
B2_threats <- B2_threats %>% filter(BSR > 10)

# calculating risk 
B2_threats$B2_risk <- B2_threats$B2 * B2_threats$BSR

## getting studies specific to birds and agriculture 
B_studies <- studies_data[studies_data$taxa == "Birds",]
B2_studies <- B_studies[B_studies$threattype == "Agriculture & aquaculture",]
B2_studies <- B2_studies[B2_studies$threatname == "Annual/perennial non-timber crops" | B2_studies$threatname == "Livestock farm/ranch" | 
                           B2_studies$threatname == "Wood/pulp plantations" | is.na(B2_studies$threatname),]

## similarly for non-english studies 
B2_non_english <- filter(non_english, taxa == "Birds" & threattype == "Agriculture & aquaculture")


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
B2_studies <- subset(B2_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
B2_studies <- st_set_geometry(B2_studies, NULL)
B2_studies <- B2_studies[!is.na(B2_studies$country),]
B2_studies <- distinct(B2_studies)
B2_counts <- count(B2_studies, country)
#B2_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
B2_non_english_counts <- count(B2_non_english, country)

#adding the two to get a complete count of studies per country 
B2_counts <- merge(B2_counts, B2_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
B2_counts[is.na(B2_counts)] <- 0
B2_counts$n <- B2_counts$n.x + B2_counts$n.y
B2_counts <- subset(B2_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
B2_centroids <- st_centroid(B2_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean B2 threat value for each country
B2_centroids_sp <- as(B2_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
B2_country <- over(world_sp, B2_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## creating a complete dataframe (this contains all countries for which there is threat data)
B2_threat_country <- B2_country$B2
B2_risk_country <- B2_country$B2_risk
B2_world_threats <- mutate(world, threats = B2_threat_country, risk = B2_risk_country)
B2_world_threats <- B2_world_threats[order(B2_world_threats$NAME),]#ordering alphabetically according to country name

B2_counts <- rename(B2_counts, NAME = country)
B2_world_studies <- merge(world, B2_counts, by = "NAME", all = TRUE)

setdiff(B2_world_studies$NAME,B2_world_threats$NAME)

B2_world_complete <- mutate(B2_world_threats, nstudies = B2_world_studies$n)
B2_world_complete$nstudies[is.na(B2_world_complete$nstudies)] <- 0
B2_world_complete <- B2_world_complete[!is.na(B2_world_complete$threats),]
B2_world_complete <- rename(B2_world_complete, Country = NAME)
B2_world_complete <- subset(B2_world_complete, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
B2_world_complete <- st_set_geometry(B2_world_complete, NULL)

#adding socioeconomic variables
B2_LM_data <- merge(B2_world_complete, EPI, by = "Country", all.x = TRUE)
B2_LM_data <- merge(B2_LM_data, GDP, by = "Country", all.x = TRUE)
B2_LM_data <- merge(B2_LM_data, NBI, by = "Country", all.x = TRUE)
B2_LM_data <- merge(B2_LM_data, WGI, by = "Country", all.x = TRUE)
B2_LM_data <- merge(B2_LM_data, lit_rates, by = "Country", all.x = TRUE)
B2_LM_data <- subset(B2_LM_data, select = -c(dataYear, pop2022))
B2_LM_data$WGI_GE_Estimate <- as.numeric(B2_LM_data$WGI_GE_Estimate)
B2_LM_data$WGI_GE_Rank <- as.numeric(B2_LM_data$WGI_GE_Rank)
B2_LM_data$literacy_rate <- as.numeric(B2_LM_data$literacy_rate)

#### BIRDS AND POLLUTION - B9 ####

## getting threats specific to birds and pollution
B9_threats <- threats_data[!is.na(threats_data$B9),]

# applying quality filter to threats data 
B9_threats <- B9_threats %>% filter(BSR > 10)

# calculating risk
B9_threats$B9_risk <- B9_threats$B9 * B9_threats$BSR

## getting studies specific to birds and pollution 
B_studies <- studies_data[studies_data$taxa == "Birds",]
B9_studies <- B_studies[B_studies$threattype == "Pollution",]

## similarly for non-english studies
B9_non_english <- filter(non_english, taxa == "Birds" & threattype == "Pollution")


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
B9_studies <- subset(B9_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
B9_studies <- st_set_geometry(B9_studies, NULL)
B9_studies <- B9_studies[!is.na(B9_studies$country),]
B9_studies <- distinct(B9_studies)
B9_counts <- count(B9_studies, country)
#B9_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
B9_non_english_counts <- count(B9_non_english, country)

#adding the two to get a complete count of studies per country 
B9_counts <- merge(B9_counts, B9_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
B9_counts[is.na(B9_counts)] <- 0
B9_counts$n <- B9_counts$n.x + B9_counts$n.y
B9_counts <- subset(B9_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
B9_centroids <- st_centroid(B9_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean B9 threat value for each country
B9_centroids_sp <- as(B9_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
B9_country <- over(world_sp, B9_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## creating a complete dataframe (this contains all countries for which there is threat data)
B9_threat_country <- B9_country$B9
B9_risk_country <- B9_country$B9_risk
B9_world_threats <- mutate(world, threats = B9_threat_country, risk = B9_risk_country)
B9_world_threats <- B9_world_threats[order(B9_world_threats$NAME),]#ordering alphabetically according to country name

B9_counts <- rename(B9_counts, NAME = country)
B9_world_studies <- merge(world, B9_counts, by = "NAME", all = TRUE)

B9_world_complete <- mutate(B9_world_threats, nstudies = B9_world_studies$n)
B9_world_complete$nstudies[is.na(B9_world_complete$nstudies)] <- 0
B9_world_complete <- B9_world_complete[!is.na(B9_world_complete$threats),]
B9_world_complete <- rename(B9_world_complete, Country = NAME)
B9_world_complete <- subset(B9_world_complete, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
B9_world_complete <- st_set_geometry(B9_world_complete, NULL)

#adding socioeconomic variables
B9_LM_data <- merge(B9_world_complete, EPI, by = "Country", all.x = TRUE)
B9_LM_data <- merge(B9_LM_data, GDP, by = "Country", all.x = TRUE)
B9_LM_data <- merge(B9_LM_data, NBI, by = "Country", all.x = TRUE)
B9_LM_data <- merge(B9_LM_data, WGI, by = "Country", all.x = TRUE)
B9_LM_data <- merge(B9_LM_data, lit_rates, by = "Country", all.x = TRUE)
B9_LM_data <- subset(B9_LM_data, select = -c(dataYear, pop2022))
B9_LM_data$WGI_GE_Estimate <- as.numeric(B9_LM_data$WGI_GE_Estimate)
B9_LM_data$WGI_GE_Rank <- as.numeric(B9_LM_data$WGI_GE_Rank)
B9_LM_data$literacy_rate <- as.numeric(B9_LM_data$literacy_rate)

####  BIRDS AND INVASIVES - B8_1 ####

## getting threats specific to birds and invasives
B8_1_threats <- threats_data[!is.na(threats_data$B8_1),]

# applying quality filter to threats data 
B8_1_threats <- B8_1_threats %>% filter(BSR > 10)

# calculating risk
B8_1_threats$B8_1_risk <- B8_1_threats$B8_1 * B8_1_threats$BSR

## getting studies specific to birds and invasives
B_studies <- studies_data[studies_data$taxa == "Birds",]
B8_1_studies <- B_studies[!is.na(B_studies$threatname),]
B8_1_studies <- B8_1_studies[B8_1_studies$threatname == "Invasive non-native species",]

## similarly for non-english studies
B8_1_non_english <- filter(non_english, taxa == "Birds" & threatname == "Invasive non-native species")


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
B8_1_studies <- subset(B8_1_studies, select = -c(region, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
B8_1_studies <- st_set_geometry(B8_1_studies, NULL)
B8_1_studies <- B8_1_studies[!is.na(B8_1_studies$country),]
B8_1_studies <- distinct(B8_1_studies)
B8_1_counts <- count(B8_1_studies, country)
#B8_1_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
B8_1_non_english_counts <- count(B8_1_non_english, country)

#adding the two to get a complete count of studies per country 
B8_1_counts <- merge(B8_1_counts, B8_1_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
B8_1_counts[is.na(B8_1_counts)] <- 0
B8_1_counts$n <- B8_1_counts$n.x + B8_1_counts$n.y
B8_1_counts <- subset(B8_1_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
B8_1_centroids <- st_centroid(B8_1_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean B8_1 threat value for each country
B8_1_centroids_sp <- as(B8_1_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
B8_1_country <- over(world_sp, B8_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## creating a complete dataframe (this contains all countries for which there is threat data)
B8_1_threat_country <- B8_1_country$B8_1
B8_1_risk_country <- B8_1_country$B8_1_risk
B8_1_world_threats <- mutate(world, threats = B8_1_threat_country, risk = B8_1_risk_country)
B8_1_world_threats <- B8_1_world_threats[order(B8_1_world_threats$NAME),]#ordering alphabetically according to country name

B8_1_counts <- rename(B8_1_counts, NAME = country)
B8_1_world_studies <- merge(world, B8_1_counts, by = "NAME", all = TRUE)

B8_1_world_complete <- mutate(B8_1_world_threats, nstudies = B8_1_world_studies$n)
B8_1_world_complete$nstudies[is.na(B8_1_world_complete$nstudies)] <- 0
B8_1_world_complete <- B8_1_world_complete[!is.na(B8_1_world_complete$threats),]
B8_1_world_complete <- rename(B8_1_world_complete, Country = NAME)
B8_1_world_complete <- subset(B8_1_world_complete, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
B8_1_world_complete <- st_set_geometry(B8_1_world_complete, NULL)

#adding socioeconomic variables
B8_1_LM_data <- merge(B8_1_world_complete, EPI, by = "Country", all.x = TRUE)
B8_1_LM_data <- merge(B8_1_LM_data, GDP, by = "Country", all.x = TRUE)
B8_1_LM_data <- merge(B8_1_LM_data, NBI, by = "Country", all.x = TRUE)
B8_1_LM_data <- merge(B8_1_LM_data, WGI, by = "Country", all.x = TRUE)
B8_1_LM_data <- merge(B8_1_LM_data, lit_rates, by = "Country", all.x = TRUE)
B8_1_LM_data <- subset(B8_1_LM_data, select = -c(dataYear, pop2022))
B8_1_LM_data$WGI_GE_Estimate <- as.numeric(B8_1_LM_data$WGI_GE_Estimate)
B8_1_LM_data$WGI_GE_Rank <- as.numeric(B8_1_LM_data$WGI_GE_Rank)
B8_1_LM_data$literacy_rate <- as.numeric(B8_1_LM_data$literacy_rate)

#### BIRDS AND LOGGING - B5_3 ####

## getting threats specific to birds and logging
B5_3_threats <- threats_data[!is.na(threats_data$B5_3),]

# applying quality filter to threats data 
B5_3_threats <- B5_3_threats %>% filter(BSR > 10)

# calculating risk
B5_3_threats$B5_3_risk <- B5_3_threats$B5_3 * B5_3_threats$BSR

## getting studies specific to birds and logging
B_studies <- studies_data[studies_data$taxa == "Birds",]
B5_3_studies <- B_studies[!is.na(B_studies$threatname),]
B5_3_studies <- B5_3_studies[B5_3_studies$threatname == "Logging/wood harvesting",]

## similarly for non-english studies 
B5_3_non_english <- filter(non_english, taxa == "Birds" & threatname == "Logging/wood harvesting")


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
B5_3_studies <- subset(B5_3_studies, select = -c(region, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
B5_3_studies <- st_set_geometry(B5_3_studies, NULL)
B5_3_studies <- B5_3_studies[!is.na(B5_3_studies$country),]
B5_3_studies <- distinct(B5_3_studies)
B5_3_counts <- count(B5_3_studies, country)
#B5_3_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
B5_3_non_english_counts <- count(B5_3_non_english, country)

#adding the two to get a complete count of studies per country 
B5_3_counts <- merge(B5_3_counts, B5_3_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
B5_3_counts[is.na(B5_3_counts)] <- 0
B5_3_counts$n <- B5_3_counts$n.x + B5_3_counts$n.y
B5_3_counts <- subset(B5_3_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
B5_3_centroids <- st_centroid(B5_3_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean B5_3 threat value for each country
B5_3_centroids_sp <- as(B5_3_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
B5_3_country <- over(world_sp, B5_3_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## creating a complete dataframe (this contains all countries for which there is threat data)
B5_3_threat_country <- B5_3_country$B5_3
B5_3_risk_country <- B5_3_country$B5_3_risk
B5_3_world_threats <- mutate(world, threats = B5_3_threat_country, risk = B5_3_risk_country)
B5_3_world_threats <- B5_3_world_threats[order(B5_3_world_threats$NAME),]#ordering alphabetically according to country name

B5_3_counts <- rename(B5_3_counts, NAME = country)
B5_3_world_studies <- merge(world, B5_3_counts, by = "NAME", all = TRUE)

B5_3_world_complete <- mutate(B5_3_world_threats, nstudies = B5_3_world_studies$n)
B5_3_world_complete$nstudies[is.na(B5_3_world_complete$nstudies)] <- 0
B5_3_world_complete <- B5_3_world_complete[!is.na(B5_3_world_complete$threats),]
B5_3_world_complete <- rename(B5_3_world_complete, Country = NAME)
B5_3_world_complete <- subset(B5_3_world_complete, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
B5_3_world_complete <- st_set_geometry(B5_3_world_complete, NULL)

#adding socioeconomic variables
B5_3_LM_data <- merge(B5_3_world_complete, EPI, by = "Country", all.x = TRUE)
B5_3_LM_data <- merge(B5_3_LM_data, GDP, by = "Country", all.x = TRUE)
B5_3_LM_data <- merge(B5_3_LM_data, NBI, by = "Country", all.x = TRUE)
B5_3_LM_data <- merge(B5_3_LM_data, WGI, by = "Country", all.x = TRUE)
B5_3_LM_data <- merge(B5_3_LM_data, lit_rates, by = "Country", all.x = TRUE)
B5_3_LM_data <- subset(B5_3_LM_data, select = -c(dataYear, pop2022))
B5_3_LM_data$WGI_GE_Estimate <- as.numeric(B5_3_LM_data$WGI_GE_Estimate)
B5_3_LM_data$WGI_GE_Rank <- as.numeric(B5_3_LM_data$WGI_GE_Rank)
B5_3_LM_data$literacy_rate <- as.numeric(B5_3_LM_data$literacy_rate)


###################################################################################
#################### birds modelling ##############################################
###################################################################################

B9_LM_data$Pollution <- B9_LM_data$risk
B2_LM_data$Agriculture <- B2_LM_data$risk
B8_1_LM_data$Invasives <- B8_1_LM_data$risk
B5_3_LM_data$Logging <- B5_3_LM_data$risk

B_LM_data <- cbind(B5_3_LM_data,Invasives=B8_1_LM_data$Invasives,Pollution=B9_LM_data$Pollution,Agriculture=B2_LM_data$Agriculture)
B_LM_data <- na.omit(B_LM_data)

# scaling variables
scale_0.5SD <- function(x){
  (x - mean(x)) / (2*sd(x))
}

B_LM_data$Invasivesscaled <- scale_0.5SD(B_LM_data$Invasives)
B_LM_data$Pollutionscaled <- scale_0.5SD(B_LM_data$Pollution)
B_LM_data$Agriculturescaled <- scale_0.5SD(B_LM_data$Agriculture)
B_LM_data$Loggingscaled <- scale_0.5SD(B_LM_data$Logging)
B_LM_data$riskscaled <- scale_0.5SD(B_LM_data$risk)
B_LM_data$EPI_newscaled <- scale_0.5SD(B_LM_data$EPI_new)
B_LM_data$GDP_2020scaled <- scale_0.5SD(B_LM_data$GDP_2020)
B_LM_data$NBIscaled <- scale_0.5SD(B_LM_data$NBI)
B_LM_data$WGI_GE_Rankscaled <- scale_0.5SD(B_LM_data$WGI_GE_Rank)
B_LM_data$literacy_ratescaled <- scale_0.5SD(B_LM_data$literacy_rate)

### poisson - results show tax model can use poisson distribution, passes DHARMa tests, zero inflation not an issue
glm_B <- glm(nstudies ~ Pollutionscaled + Agriculturescaled + Invasivesscaled + Loggingscaled + EPI_newscaled + GDP_2020scaled + NBIscaled + WGI_GE_Rankscaled + literacy_ratescaled, family = "poisson", data = B_LM_data)
summary(glm_B)

head(B_LM_data)

testDispersion(glm_B)
simulationOutput <- simulateResiduals(fittedModel = glm_B, plot = F)
residuals(simulationOutput)
plot(simulationOutput)

###### create global model and dredge to find best model
options(na.action = "na.fail") #Must run this code once to use dredge
all_models_b <- dredge(glm_B)
#write.csv(data.frame(all_models_b),"all_models_b.csv")

all_models_b
nrow(all_models_b)

#several models with deltaAICc <2
top.models.2aic <- get.models(all_models_b, subset = delta <2)
length(top.models.2aic)
top.models.4aic <- get.models(all_models_b, subset = delta <4)
length(top.models.4aic)
top.models.95ci <- get.models(all_models_b, cumsum(weight) <= 0.95)
length(top.models.95ci)

all_models_b_avg <- summary(model.avg(top.models.2aic))

all_models_b_avg
sw(all_models_b_avg)
confint(all_models_b_avg, full = T)

#write.csv(summary(all_models_b_avg)[9],"glm_Bnb_summaryupdate.csv")
#write.csv(sw(all_models_b_avg),"glm_Bnb_relimp.csv")
#write.csv(confint(all_models_b_avg,full=TRUE),"glm_Bnb_confint.csv")


########################################### MAMMALS ###########################################
########################################### MAMMALS ###########################################



#### MAMMALS AND AGRICULTURE - M2 ####

## getting threats specific to mammals and agriculture
M2_threats <- threats_data[!is.na(threats_data$M2),]

# applying quality filter to threats data 
M2_threats <- M2_threats %>% filter(MSR > 10)

# calculating risk
M2_threats$M2_risk <- M2_threats$M2 * M2_threats$MSR

## getting studies specific to mammals and agriculture 
M_studies <- studies_data[studies_data$taxa == "Mammals",]
M2_studies <- M_studies[M_studies$threattype == "Agriculture & aquaculture",]
M2_studies <- M2_studies[M2_studies$threatname == "Annual/perennial non-timber crops" | M2_studies$threatname == "Livestock farm/ranch" | 
                           M2_studies$threatname == "Wood/pulp plantations" | is.na(M2_studies$threatname),]

## similarly for non-english studies
M2_non_english <- filter(non_english, taxa == "Mammals" & threattype == "Agriculture & aquaculture")


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
M2_studies <- subset(M2_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
M2_studies <- st_set_geometry(M2_studies, NULL)
M2_studies <- M2_studies[!is.na(M2_studies$country),]
M2_studies <- distinct(M2_studies)
M2_counts <- count(M2_studies, country)
#M2_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
M2_non_english_counts <- count(M2_non_english, country)

#adding the two to get a complete count of studies per country 
M2_counts <- merge(M2_counts, M2_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
M2_counts[is.na(M2_counts)] <- 0
M2_counts$n <- M2_counts$n.x + M2_counts$n.y
M2_counts <- subset(M2_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
M2_centroids <- st_centroid(M2_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean M2 threat value for each country
M2_centroids_sp <- as(M2_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
M2_country <- over(world_sp, M2_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## creating a complete dataframe (this contains all countries for which there is threat data)
M2_threat_country <- M2_country$M2
M2_risk_country <- M2_country$M2_risk
M2_world_threats <- mutate(world, threats = M2_threat_country, risk = M2_risk_country)
M2_world_threats <- M2_world_threats[order(M2_world_threats$NAME),]#ordering alphabetically according to country name

M2_counts <- rename(M2_counts, NAME = country)
M2_world_studies <- merge(world, M2_counts, by = "NAME", all = TRUE)

M2_world_complete <- mutate(M2_world_threats, nstudies = M2_world_studies$n)
M2_world_complete$nstudies[is.na(M2_world_complete$nstudies)] <- 0
M2_world_complete <- M2_world_complete[!is.na(M2_world_complete$threats),]
M2_world_complete <- rename(M2_world_complete, Country = NAME)
M2_world_complete <- subset(M2_world_complete, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
M2_world_complete <- st_set_geometry(M2_world_complete, NULL)

#adding socioeconomic variables
M2_LM_data <- merge(M2_world_complete, EPI, by = "Country", all.x = TRUE)
M2_LM_data <- merge(M2_LM_data, GDP, by = "Country", all.x = TRUE)
M2_LM_data <- merge(M2_LM_data, NBI, by = "Country", all.x = TRUE)
M2_LM_data <- merge(M2_LM_data, WGI, by = "Country", all.x = TRUE)
M2_LM_data <- merge(M2_LM_data, lit_rates, by = "Country", all.x = TRUE)
M2_LM_data <- subset(M2_LM_data, select = -c(dataYear, pop2022))
M2_LM_data$WGI_GE_Estimate <- as.numeric(M2_LM_data$WGI_GE_Estimate)
M2_LM_data$WGI_GE_Rank <- as.numeric(M2_LM_data$WGI_GE_Rank)
M2_LM_data$literacy_rate <- as.numeric(M2_LM_data$literacy_rate)

#### MAMMALS AND POLLUTION - M9 ####

## getting threats specific to mammals and pollution
M9_threats <- threats_data[!is.na(threats_data$M9),]

# applying quality filter to threats data 
M9_threats <- M9_threats %>% filter(MSR > 10)

# calculating risk
M9_threats$M9_risk <- M9_threats$M9 * M9_threats$MSR

## getting studies specific to mammals and pollution
M_studies <- studies_data[studies_data$taxa == "Mammals",]
M9_studies <- M_studies[M_studies$threattype == "Pollution",]


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
M9_studies <- subset(M9_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
M9_studies <- st_set_geometry(M9_studies, NULL)
M9_studies <- M9_studies[!is.na(M9_studies$country),]
M9_studies <- distinct(M9_studies)
M9_counts <- count(M9_studies, country)
#M9_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
M9_centroids <- st_centroid(M9_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean M9 threat value for each country
M9_centroids_sp <- as(M9_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
M9_country <- over(world_sp, M9_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## creating a complete dataframe (this contains all countries for which there is threat data)
M9_threat_country <- M9_country$M9
M9_risk_country <- M9_country$M9_risk
M9_world_threats <- mutate(world, threats = M9_threat_country, risk = M9_risk_country)
M9_world_threats <- M9_world_threats[order(M9_world_threats$NAME),]#ordering alphabetically according to country name

M9_counts <- rename(M9_counts, NAME = country)
M9_world_studies <- merge(world, M9_counts, by = "NAME", all = TRUE)

M9_world_complete <- mutate(M9_world_threats, nstudies = M9_world_studies$n)
M9_world_complete$nstudies[is.na(M9_world_complete$nstudies)] <- 0
M9_world_complete <- M9_world_complete[!is.na(M9_world_complete$threats),]
M9_world_complete <- rename(M9_world_complete, Country = NAME)
M9_world_complete <- subset(M9_world_complete, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
M9_world_complete <- st_set_geometry(M9_world_complete, NULL)

#adding socioeconomic variables
M9_LM_data <- merge(M9_world_complete, EPI, by = "Country", all.x = TRUE)
M9_LM_data <- merge(M9_LM_data, GDP, by = "Country", all.x = TRUE)
M9_LM_data <- merge(M9_LM_data, NBI, by = "Country", all.x = TRUE)
M9_LM_data <- merge(M9_LM_data, WGI, by = "Country", all.x = TRUE)
M9_LM_data <- merge(M9_LM_data, lit_rates, by = "Country", all.x = TRUE)
M9_LM_data <- subset(M9_LM_data, select = -c(dataYear, pop2022))
M9_LM_data$WGI_GE_Estimate <- as.numeric(M9_LM_data$WGI_GE_Estimate)
M9_LM_data$WGI_GE_Rank <- as.numeric(M9_LM_data$WGI_GE_Rank)
M9_LM_data$literacy_rate <- as.numeric(M9_LM_data$literacy_rate)


#### MAMMALS AND HUNTING - M5_1 ####

## getting threats specific to mammals and hunting
M5_1_threats <- threats_data[!is.na(threats_data$M5_1),]

# applying quality filter to threats data 
M5_1_threats <- M5_1_threats %>% filter(MSR > 10)

# calculating risk
M5_1_threats$M5_1_risk <- M5_1_threats$M5_1 * M5_1_threats$MSR

## getting studies specific to mammals and hunting
M_studies <- studies_data[studies_data$taxa == "Mammals",]
M5_1_studies <- M_studies[!is.na(M_studies$threatname),]
M5_1_studies <- M5_1_studies[M5_1_studies$threatname == "Hunting/trapping terrestrial animals",]

## similarly for non-english studies
M5_1_non_english <- filter(non_english, taxa == "Mammals" & threatname == "Hunting/trapping terrestrial animals")


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
M5_1_studies <- subset(M5_1_studies, select = -c(region, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
M5_1_studies <- st_set_geometry(M5_1_studies, NULL)
M5_1_studies <- M5_1_studies[!is.na(M5_1_studies$country),]
M5_1_studies <- distinct(M5_1_studies)
M5_1_counts <- count(M5_1_studies, country)
#M5_1_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
M5_1_non_english_counts <- count(M5_1_non_english, country)

#adding the two to get a complete count of studies per country 
M5_1_counts <- merge(M5_1_counts, M5_1_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
M5_1_counts[is.na(M5_1_counts)] <- 0
M5_1_counts$n <- M5_1_counts$n.x + M5_1_counts$n.y
M5_1_counts <- subset(M5_1_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
M5_1_centroids <- st_centroid(M5_1_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean M5_1 threat value for each country
M5_1_centroids_sp <- as(M5_1_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
M5_1_country <- over(world_sp, M5_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## creating a complete dataframe (this contains all countries for which there is threat data)
M5_1_threat_country <- M5_1_country$M5_1
M5_1_risk_country <- M5_1_country$M5_1_risk
M5_1_world_threats <- mutate(world, threats = M5_1_threat_country, risk = M5_1_risk_country)
M5_1_world_threats <- M5_1_world_threats[order(M5_1_world_threats$NAME),]#ordering alphabetically according to country name

M5_1_counts <- rename(M5_1_counts, NAME = country)
M5_1_world_studies <- merge(world, M5_1_counts, by = "NAME", all = TRUE)

M5_1_world_complete <- mutate(M5_1_world_threats, nstudies = M5_1_world_studies$n)
M5_1_world_complete$nstudies[is.na(M5_1_world_complete$nstudies)] <- 0
M5_1_world_complete <- M5_1_world_complete[!is.na(M5_1_world_complete$threats),]
M5_1_world_complete <- rename(M5_1_world_complete, Country = NAME)
M5_1_world_complete <- subset(M5_1_world_complete, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
M5_1_world_complete <- st_set_geometry(M5_1_world_complete, NULL)

#adding socioeconomic variables
M5_1_LM_data <- merge(M5_1_world_complete, EPI, by = "Country", all.x = TRUE)
M5_1_LM_data <- merge(M5_1_LM_data, GDP, by = "Country", all.x = TRUE)
M5_1_LM_data <- merge(M5_1_LM_data, NBI, by = "Country", all.x = TRUE)
M5_1_LM_data <- merge(M5_1_LM_data, WGI, by = "Country", all.x = TRUE)
M5_1_LM_data <- merge(M5_1_LM_data, lit_rates, by = "Country", all.x = TRUE)
M5_1_LM_data <- subset(M5_1_LM_data, select = -c(dataYear, pop2022))
M5_1_LM_data$WGI_GE_Estimate <- as.numeric(M5_1_LM_data$WGI_GE_Estimate)
M5_1_LM_data$WGI_GE_Rank <- as.numeric(M5_1_LM_data$WGI_GE_Rank)
M5_1_LM_data$literacy_rate <- as.numeric(M5_1_LM_data$literacy_rate)




#### MAMMALS AND LOGGING - M5_3 ####

## getting threats specific to mammals and logging
M5_3_threats <- threats_data[!is.na(threats_data$M5_3),]

# applying quality filter to threats data 
M5_3_threats <- M5_3_threats %>% filter(MSR > 10)

# calculating risk
M5_3_threats$M5_3_risk <- M5_3_threats$M5_3 * M5_3_threats$MSR

## getting studies specific to mammals and logging
M_studies <- studies_data[studies_data$taxa == "Mammals",]
M5_3_studies <- M_studies[!is.na(M_studies$threatname),]
M5_3_studies <- M5_3_studies[M5_3_studies$threatname == "Logging/wood harvesting",]

## similarly for non-english studies 
M5_3_non_english <- filter(non_english, taxa == "Mammals" & threatname == "Logging/wood harvesting")


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
M5_3_studies <- subset(M5_3_studies, select = -c(region, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
M5_3_studies <- st_set_geometry(M5_3_studies, NULL)
M5_3_studies <- M5_3_studies[!is.na(M5_3_studies$country),]
M5_3_studies <- distinct(M5_3_studies)
M5_3_counts <- count(M5_3_studies, country)
#M5_3_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
M5_3_non_english_counts <- count(M5_3_non_english, country)

#adding the two to get a complete count of studies per country 
M5_3_counts <- merge(M5_3_counts, M5_3_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
M5_3_counts[is.na(M5_3_counts)] <- 0
M5_3_counts$n <- M5_3_counts$n.x + M5_3_counts$n.y
M5_3_counts <- subset(M5_3_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
M5_3_centroids <- st_centroid(M5_3_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean M5_3 threat value for each country
M5_3_centroids_sp <- as(M5_3_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
M5_3_country <- over(world_sp, M5_3_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## creating a complete dataframe (this contains all countries for which there is threat data)
M5_3_threat_country <- M5_3_country$M5_3
M5_3_risk_country <- M5_3_country$M5_3_risk
M5_3_world_threats <- mutate(world, threats = M5_3_threat_country, risk = M5_3_risk_country)
M5_3_world_threats <- M5_3_world_threats[order(M5_3_world_threats$NAME),]#ordering alphabetically according to country name

M5_3_counts <- rename(M5_3_counts, NAME = country)
M5_3_world_studies <- merge(world, M5_3_counts, by = "NAME", all = TRUE)

M5_3_world_complete <- mutate(M5_3_world_threats, nstudies = M5_3_world_studies$n)
M5_3_world_complete$nstudies[is.na(M5_3_world_complete$nstudies)] <- 0
M5_3_world_complete <- M5_3_world_complete[!is.na(M5_3_world_complete$threats),]
M5_3_world_complete <- rename(M5_3_world_complete, Country = NAME)
M5_3_world_complete <- subset(M5_3_world_complete, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
M5_3_world_complete <- st_set_geometry(M5_3_world_complete, NULL)

#adding socioeconomic variables
M5_3_LM_data <- merge(M5_3_world_complete, EPI, by = "Country", all.x = TRUE)
M5_3_LM_data <- merge(M5_3_LM_data, GDP, by = "Country", all.x = TRUE)
M5_3_LM_data <- merge(M5_3_LM_data, NBI, by = "Country", all.x = TRUE)
M5_3_LM_data <- merge(M5_3_LM_data, WGI, by = "Country", all.x = TRUE)
M5_3_LM_data <- merge(M5_3_LM_data, lit_rates, by = "Country", all.x = TRUE)
M5_3_LM_data <- subset(M5_3_LM_data, select = -c(dataYear, pop2022))
M5_3_LM_data$WGI_GE_Estimate <- as.numeric(M5_3_LM_data$WGI_GE_Estimate)
M5_3_LM_data$WGI_GE_Rank <- as.numeric(M5_3_LM_data$WGI_GE_Rank)
M5_3_LM_data$literacy_rate <- as.numeric(M5_3_LM_data$literacy_rate)

#######################################################################
############################# mammals modelling #######################
#######################################################################

M9_LM_data$Pollution <- M9_LM_data$risk
M2_LM_data$Agriculture <- M2_LM_data$risk
M5_1_LM_data$Hunting <- M5_1_LM_data$risk
M5_3_LM_data$Logging <- M5_3_LM_data$risk

M_LM_data <- cbind(M5_3_LM_data,Hunting=M5_1_LM_data$Hunting,Pollution=M9_LM_data$Pollution,Agriculture=M2_LM_data$Agriculture)
M_LM_data <- na.omit(M_LM_data)

#scale variables
scale_0.5SD <- function(x){
  (x - mean(x)) / (2*sd(x))
}

M_LM_data$Pollutionscaled <- scale_0.5SD(M_LM_data$Pollution)
M_LM_data$Agriculturescaled <- scale_0.5SD(M_LM_data$Agriculture)
M_LM_data$Loggingscaled <- scale_0.5SD(M_LM_data$Logging)
M_LM_data$Huntingscaled <- scale_0.5SD(M_LM_data$Hunting)
M_LM_data$riskscaled <- scale_0.5SD(M_LM_data$risk)
M_LM_data$EPI_newscaled <- scale_0.5SD(M_LM_data$EPI_new)
M_LM_data$GDP_2020scaled <- scale_0.5SD(M_LM_data$GDP_2020)
M_LM_data$NBIscaled <- scale_0.5SD(M_LM_data$NBI)
M_LM_data$WGI_GE_Rankscaled <- scale_0.5SD(M_LM_data$WGI_GE_Rank)
M_LM_data$literacy_ratescaled <- scale_0.5SD(M_LM_data$literacy_rate)

### poisson - results show tax model cannot use poisson distribution, does not pass DHARMa tests, overdispersion/zero inflation an issue
glm_M <- glm(nstudies ~ Huntingscaled + Pollutionscaled + Agriculturescaled + Loggingscaled + EPI_newscaled + GDP_2020scaled + NBIscaled + WGI_GE_Rankscaled + literacy_ratescaled, family = "poisson", data = M_LM_data)
summary(glm_M)

head(M_LM_data)
testDispersion(glm_M)
simulationOutput <- simulateResiduals(fittedModel = glm_M, plot = F)
residuals(simulationOutput)
plot(simulationOutput)

### quasi-POISSON REGRESSION
### Tests for the quasipoisson distribution for Tax model - not symmetrical
glm_Mqp <- glm(nstudies ~ Huntingscaled + Pollutionscaled + Agriculturescaled + Loggingscaled + EPI_newscaled + GDP_2020scaled + NBIscaled + WGI_GE_Rankscaled + literacy_ratescaled, family = "quasipoisson", data = M_LM_data)
summary(glm_Mqp)
dev_residuals <- residuals(glm_Mqp, type = "deviance")
# Plot deviance residuals against predicted values
plot(fitted(glm_Mqp), dev_residuals, main = "Deviance Residuals vs. Fitted Values")

# use negative binomial that accounts for zeron inflation instead as poisson and quasi-poisson not suitable
glm_Mnb <- glmmTMB(nstudies ~ Huntingscaled + Pollutionscaled + Agriculturescaled + Loggingscaled + EPI_newscaled + GDP_2020scaled + NBIscaled + WGI_GE_Rankscaled + literacy_ratescaled, family=nbinom2, ziformula = ~0, data = M_LM_data)
summary(glm_Mnb)

# create global model and dredge to find best model
options(na.action = "na.fail") #Must run this code once to use dredge
all_models_m <- dredge(glm_Mnb)
#write.csv(data.frame(all_models_m),"all_models_m.csv")

all_models_m
nrow(all_models_m)

#several models with deltaAICc <2
m.top.models.2aic <- get.models(all_models_m, subset = delta <2)
length(m.top.models.2aic)
m.top.models.4aic <- get.models(all_models_m, subset = delta <4)
length(m.top.models.4aic)
m.top.models.95ci <- get.models(all_models_m, cumsum(weight) <= 0.95)
length(m.top.models.95ci)

all_models_m_avg <- summary(model.avg(m.top.models.2aic))

all_models_m_avg
sw(all_models_m_avg)
confint(all_models_m_avg,full=TRUE)


# write.csv(summary(all_models_m_avg)[9],"glm_Mnb_summaryupdate.csv")
# write.csv(sw(all_models_m_avg),"glm_Mnb_relimp.csv")
# write.csv(confint(all_models_m_avg,full=TRUE),"glm_Mnb_confint.csv")
