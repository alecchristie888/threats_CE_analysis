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
library(performance)
library(biscale)
library(cowplot)

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
A2_country_ncells <- over(world_sp, A2_centroids_sp, fn = length)
A2_country$ncells <- A2_country_ncells$A2_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
A2_risk_country <- A2_country$A2_risk
A2_threat_country <- A2_country$A2
A2_threat_ncells <- A2_country$ncells
A2_world_threats <- mutate(world, threats = A2_threat_country, risk = A2_risk_country, ncells=A2_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
A2_world_threats$cellstosize <- as.numeric((A2_world_threats$ncells*2500)/(A2_world_threats$AREA*10))
A2_world_threats[order(A2_world_threats$cellstosize),c("NAME","cellstosize")]
A2_world_threats <- A2_world_threats[!is.na(A2_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
A2_world_threats <- A2_world_threats[A2_world_threats$cellstosize>=0.10,]

A2_world_complete <- merge(A2_world_threats, A2_counts[,c("country","n")], by.x="NAME",by.y="country",all=TRUE)
A2_world_complete$n[is.na(A2_world_complete$n)] <- 0
A2_world_complete <- A2_world_complete[!is.na(A2_world_complete$threats),]
unique(A2_world_complete$NAME)

A2_world_complete$nstudies <- A2_world_complete$n
#bivariate maps
A2_world_complete$threats_bin <- cut(A2_world_complete$risk, breaks = c(0, 5, 10, 15, 38), include.lowest = TRUE)
A2_world_complete$nstudies_bin <- cut(A2_world_complete$n, breaks = c(0, 0.5, 6, 15, 30), include.lowest = TRUE)
classes <- bi_class(A2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(A2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A2 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("A2_risk.svg", height = 10, width = 12, bg = "white",device="svg")


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
A9_country_ncells <- over(world_sp, A9_centroids_sp, fn = length)
A9_country$ncells <- A9_country_ncells$A9_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
A9_risk_country <- A9_country$A9_risk
A9_threat_country <- A9_country$A9
A9_threat_ncells <- A9_country$ncells
A9_world_threats <- mutate(world, threats = A9_threat_country, risk = A9_risk_country, ncells=A9_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
A9_world_threats$cellstosize <- as.numeric((A9_world_threats$ncells*2500)/(A9_world_threats$AREA*10))
as.data.frame(A9_world_threats[order(A9_world_threats$cellstosize),c("NAME","cellstosize")])
A9_world_threats <- A9_world_threats[!is.na(A9_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
A9_world_threats <- A9_world_threats[A9_world_threats$cellstosize>=0.10,]

A9_world_complete <- merge(A9_world_threats, A9_counts[,c("country","n")], by.x="NAME", by.y="country",all=TRUE)
A9_world_complete$n[is.na(A9_world_complete$n)] <- 0
A9_world_complete <- A9_world_complete[!is.na(A9_world_complete$threats),]
unique(A9_world_complete$NAME)

A9_world_complete$nstudies <- A9_world_complete$n
#bivariate maps
A9_world_complete$threats_bin <- cut(A9_world_complete$risk, breaks = c(0, 0.50, 1, 3, 7), include.lowest = TRUE)
A9_world_complete$nstudies_bin <- cut(A9_world_complete$nstudies, breaks = c(0, 0.5, 1, 5, 9), include.lowest = TRUE)
classes <- bi_class(A9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(A9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A9 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("A9_risk.svg", height = 10, width = 12, bg = "white",device="svg")


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

unique(A_LM_data$Country)
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
write.csv(data.frame(all_models_a),"all_models_a.csv")

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

#calculate the VIF for each predictor variable in the final models
vifmodel <- glm(nstudies ~ Agriculturescaled + Pollutionscaled + EPI_newscaled + GDP_2020scaled + NBIscaled + literacy_ratescaled, family = "poisson", data = A_LM_data)

check_collinearity(
     vifmodel,
     component = c('all') # 'all' shows both conditional and zi components
   ) #all below 5, low-moderate multicollinearity, not severe.

 write.csv(summary(all_models_a_avg)[9],"glm_Anb_summaryupdate.csv")
 write.csv(sw(all_models_a_avg),"glm_Anb_relimp.csv")
 write.csv(confint(all_models_a_avg,full=TRUE),"glm_Anb_confint.csv")

#robustness check
backward_model <- step(glm_A, direction = "backward")

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
B2_country_ncells <- over(world_sp, B2_centroids_sp, fn = length)
B2_country$ncells <- B2_country_ncells$B2_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
B2_risk_country <- B2_country$B2_risk
B2_threat_country <- B2_country$B2
B2_threat_ncells <- B2_country$ncells
B2_world_threats <- mutate(world, threats = B2_threat_country, risk = B2_risk_country, ncells=B2_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
B2_world_threats$cellstosize <- as.numeric((B2_world_threats$ncells*2500)/(B2_world_threats$AREA*10))
as.data.frame(B2_world_threats[order(B2_world_threats$cellstosize),c("NAME","cellstosize")])
B2_world_threats <- B2_world_threats[!is.na(B2_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
B2_world_threats <- B2_world_threats[B2_world_threats$cellstosize>=0.10,]

B2_world_complete <- merge(B2_world_threats, B2_counts[,c("country","n")], by.x="NAME", by.y="country",all=TRUE)
B2_world_complete$n[is.na(B2_world_complete$n)] <- 0
B2_world_complete <- B2_world_complete[!is.na(B2_world_complete$threats),]
unique(B2_world_complete$NAME)

B2_world_complete$nstudies <- B2_world_complete$n
#bivariate maps
B2_world_complete$threats_bin <- cut(B2_world_complete$risk, breaks = c(0, 5, 20, 60, 115), include.lowest = TRUE)
B2_world_complete$nstudies_bin <- cut(B2_world_complete$nstudies, breaks = c(0, 0.5, 20, 75, 139), include.lowest = TRUE)
classes <- bi_class(B2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(B2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_B2 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("B2_risk.svg", height = 10, width = 12, bg = "white",device="svg")


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
B9_country_ncells <- over(world_sp, B9_centroids_sp, fn = length)
B9_country$ncells <- B9_country_ncells$B9_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
B9_risk_country <- B9_country$B9_risk
B9_threat_country <- B9_country$B9
B9_threat_ncells <- B9_country$ncells
B9_world_threats <- mutate(world, threats = B9_threat_country, risk = B9_risk_country, ncells=B9_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
B9_world_threats$cellstosize <- as.numeric((B9_world_threats$ncells*2500)/(B9_world_threats$AREA*10))
as.data.frame(B9_world_threats[order(B9_world_threats$cellstosize),c("NAME","cellstosize")])
B9_world_threats <- B9_world_threats[!is.na(B9_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
B9_world_threats <- B9_world_threats[B9_world_threats$cellstosize>=0.10,]

B9_world_complete <- merge(B9_world_threats, B9_counts[,c("country","n")], by.x="NAME", by.y="country",all=TRUE)
B9_world_complete$n[is.na(B9_world_complete$n)] <- 0
B9_world_complete <- B9_world_complete[!is.na(B9_world_complete$threats),]
unique(B9_world_complete$NAME)

B9_world_complete$nstudies <- B9_world_complete$n
#bivariate maps
B9_world_complete$threats_bin <- cut(B9_world_complete$risk, breaks = c(0, 0.50, 2, 5, 8), include.lowest = TRUE)
B9_world_complete$nstudies_bin <- cut(B9_world_complete$nstudies, breaks = c(0, 0.5, 5, 12, 24), include.lowest = TRUE)
classes <- bi_class(B9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(B9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_B9 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("B9_risk.svg", height = 10, width = 12, bg = "white",device="svg")


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
B8_1_country_ncells <- over(world_sp, B8_1_centroids_sp, fn = length)
B8_1_country$ncells <- B8_1_country_ncells$B8_1_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
B8_1_risk_country <- B8_1_country$B8_1_risk
B8_1_threat_country <- B8_1_country$B8_1
B8_1_threat_ncells <- B8_1_country$ncells
B8_1_world_threats <- mutate(world, threats = B8_1_threat_country, risk = B8_1_risk_country, ncells=B8_1_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
B8_1_world_threats$cellstosize <- as.numeric((B8_1_world_threats$ncells*2500)/(B8_1_world_threats$AREA*10))
as.data.frame(B8_1_world_threats[order(B8_1_world_threats$cellstosize),c("NAME","cellstosize")])
B8_1_world_threats <- B8_1_world_threats[!is.na(B8_1_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
B8_1_world_threats <- B8_1_world_threats[B8_1_world_threats$cellstosize>=0.10,]

B8_1_world_complete <- merge(B8_1_world_threats, B8_1_counts[,c("country","n")], by.x="NAME", by.y="country",all=TRUE)
B8_1_world_complete$n[is.na(B8_1_world_complete$n)] <- 0
B8_1_world_complete <- B8_1_world_complete[!is.na(B8_1_world_complete$threats),]
unique(B8_1_world_complete$NAME)

B8_1_world_complete$nstudies <- B8_1_world_complete$n
#bivariate maps
B8_1_world_complete$threats_bin <- cut(B8_1_world_complete$risk, breaks = c(0, 5, 10, 20, 51), include.lowest = TRUE)
B8_1_world_complete$nstudies_bin <- cut(B8_1_world_complete$nstudies, breaks = c(0, 0.5, 1, 3, 5), include.lowest = TRUE)
classes <- bi_class(B8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(B8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_B8_1 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("B8_1_risk.svg", height = 10, width = 12, bg = "white",device="svg")


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
B5_3_country_ncells <- over(world_sp, B5_3_centroids_sp, fn = length)
B5_3_country$ncells <- B5_3_country_ncells$B5_3_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
B5_3_risk_country <- B5_3_country$B5_3_risk
B5_3_threat_country <- B5_3_country$B5_3
B5_3_threat_ncells <- B5_3_country$ncells
B5_3_world_threats <- mutate(world, threats = B5_3_threat_country, risk = B5_3_risk_country, ncells=B5_3_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
B5_3_world_threats$cellstosize <- as.numeric((B5_3_world_threats$ncells*2500)/(B5_3_world_threats$AREA*10))
as.data.frame(B5_3_world_threats[order(B5_3_world_threats$cellstosize),c("NAME","cellstosize")])
B5_3_world_threats <- B5_3_world_threats[!is.na(B5_3_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
B5_3_world_threats <- B5_3_world_threats[B5_3_world_threats$cellstosize>=0.10,]

B5_3_world_complete <- merge(B5_3_world_threats, B5_3_counts[,c("country","n")], by.x="NAME", by.y="country",all=TRUE)
B5_3_world_complete$n[is.na(B5_3_world_complete$n)] <- 0
B5_3_world_complete <- B5_3_world_complete[!is.na(B5_3_world_complete$threats),]
unique(B5_3_world_complete$NAME)

B5_3_world_complete$nstudies <- B5_3_world_complete$n
#bivariate maps
B5_3_world_complete$threats_bin <- cut(B5_3_world_complete$risk, breaks = c(0, 5, 15, 50, 112), include.lowest = TRUE)
B5_3_world_complete$nstudies_bin <- cut(B5_3_world_complete$nstudies, breaks = c(0, 0.5, 2, 6, 15), include.lowest = TRUE)
classes <- bi_class(B5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(B5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_B5_3 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("B5_3_risk.svg", height = 10, width = 12, bg = "white",device="svg")


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

#calculate the VIF for each predictor variable in the final models
vifmodel <- glm(nstudies ~ Pollutionscaled + Invasivesscaled + EPI_newscaled + GDP_2020scaled + NBIscaled + WGI_GE_Rankscaled + literacy_ratescaled, family = "poisson", data = B_LM_data)
check_collinearity(
  vifmodel,
  component = c('all') # 'all' shows both conditional and zi components
) #all below 5, low-moderate multicollinearity, not severe.

#write.csv(summary(all_models_b_avg)[9],"glm_Bnb_summaryupdate.csv")
#write.csv(sw(all_models_b_avg),"glm_Bnb_relimp.csv")
#write.csv(confint(all_models_b_avg,full=TRUE),"glm_Bnb_confint.csv")

#robustness check
backward_model <- step(glm_B, direction = "backward")

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
M2_country_ncells <- over(world_sp, M2_centroids_sp, fn = length)
M2_country$ncells <- M2_country_ncells$M2_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
M2_risk_country <- M2_country$M2_risk
M2_threat_country <- M2_country$M2
M2_threat_ncells <- M2_country$ncells
M2_world_threats <- mutate(world, threats = M2_threat_country, risk = M2_risk_country, ncells=M2_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
M2_world_threats$cellstosize <- as.numeric((M2_world_threats$ncells*2500)/(M2_world_threats$AREA*10))
as.data.frame(M2_world_threats[order(M2_world_threats$cellstosize),c("NAME","cellstosize")])
M2_world_threats <- M2_world_threats[!is.na(M2_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
M2_world_threats <- M2_world_threats[M2_world_threats$cellstosize>=0.10,]

M2_world_complete <- merge(M2_world_threats, M2_counts[,c("country","n")], by.x="NAME", by.y="country",all=TRUE)
M2_world_complete$n[is.na(M2_world_complete$n)] <- 0
M2_world_complete <- M2_world_complete[!is.na(M2_world_complete$threats),]
unique(M2_world_complete$NAME)

M2_world_complete$nstudies <- M2_world_complete$n

summary(M2_world_complete$risk)
summary(M2_world_complete$nstudies)
#bivariate maps
M2_world_complete$threats_bin <- cut(M2_world_complete$risk, breaks = c(0, 5, 20, 50, 84), include.lowest = TRUE)
M2_world_complete$nstudies_bin <- cut(M2_world_complete$nstudies, breaks = c(0, 0.5,5, 50, 113), include.lowest = TRUE)
classes <- bi_class(M2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(M2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_M2 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("M2_risk.svg", height = 10, width = 12, bg = "white",device="svg")


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
M9_country_ncells <- over(world_sp, M9_centroids_sp, fn = length)
M9_country$ncells <- M9_country_ncells$M9_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
M9_risk_country <- M9_country$M9_risk
M9_threat_country <- M9_country$M9
M9_threat_ncells <- M9_country$ncells
M9_world_threats <- mutate(world, threats = M9_threat_country, risk = M9_risk_country, ncells=M9_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
M9_world_threats$cellstosize <- as.numeric((M9_world_threats$ncells*2500)/(M9_world_threats$AREA*10))
as.data.frame(M9_world_threats[order(M9_world_threats$cellstosize),c("NAME","cellstosize")])
M9_world_threats <- M9_world_threats[!is.na(M9_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
M9_world_threats <- M9_world_threats[M9_world_threats$cellstosize>=0.10,]

M9_world_complete <- merge(M9_world_threats, M9_counts[,c("country","n")], by.x="NAME", by.y="country",all=TRUE)
M9_world_complete$n[is.na(M9_world_complete$n)] <- 0
M9_world_complete <- M9_world_complete[!is.na(M9_world_complete$threats),]
unique(M9_world_complete$NAME)

M9_world_complete$nstudies <- M9_world_complete$n
#bivariate maps
M9_world_complete$threats_bin <- cut(M9_world_complete$risk, breaks = c(0, 0.5, 1, 2, 5), include.lowest = TRUE)
M9_world_complete$nstudies_bin <- cut(M9_world_complete$nstudies, breaks = c(0, 0.5, 3, 10, 16), include.lowest = TRUE)
classes <- bi_class(M9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(M9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_M9 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("M9_risk.svg", height = 10, width = 12, bg = "white",device="svg")


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
M5_1_country_ncells <- over(world_sp, M5_1_centroids_sp, fn = length)
M5_1_country$ncells <- M5_1_country_ncells$M5_1_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
M5_1_risk_country <- M5_1_country$M5_1_risk
M5_1_threat_country <- M5_1_country$M5_1
M5_1_threat_ncells <- M5_1_country$ncells
M5_1_world_threats <- mutate(world, threats = M5_1_threat_country, risk = M5_1_risk_country, ncells=M5_1_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
M5_1_world_threats$cellstosize <- as.numeric((M5_1_world_threats$ncells*2500)/(M5_1_world_threats$AREA*10))
as.data.frame(M5_1_world_threats[order(M5_1_world_threats$cellstosize),c("NAME","cellstosize")])
M5_1_world_threats <- M5_1_world_threats[!is.na(M5_1_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
M5_1_world_threats <- M5_1_world_threats[M5_1_world_threats$cellstosize>=0.10,]

M5_1_world_complete <- merge(M5_1_world_threats, M5_1_counts[,c("country","n")], by.x="NAME", by.y="country",all=TRUE)
M5_1_world_complete$n[is.na(M5_1_world_complete$n)] <- 0
M5_1_world_complete <- M5_1_world_complete[!is.na(M5_1_world_complete$threats),]
unique(M5_1_world_complete$NAME)

M5_1_world_complete$nstudies <- M5_1_world_complete$n
#bivariate maps
M5_1_world_complete$threats_bin <- cut(M5_1_world_complete$risk, breaks = c(0, 5, 12, 25, 46), include.lowest = TRUE)
M5_1_world_complete$nstudies_bin <- cut(M5_1_world_complete$nstudies, breaks = c(0, 0.5, 5, 20, 41), include.lowest = TRUE)
classes <- bi_class(M5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(M5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_M5_1 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("M5_1_risk.svg", height = 10, width = 12, bg = "white",device="svg")


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
M5_3_country_ncells <- over(world_sp, M5_3_centroids_sp, fn = length)
M5_3_country$ncells <- M5_3_country_ncells$M5_3_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
M5_3_risk_country <- M5_3_country$M5_3_risk
M5_3_threat_country <- M5_3_country$M5_3
M5_3_threat_ncells <- M5_3_country$ncells
M5_3_world_threats <- mutate(world, threats = M5_3_threat_country, risk = M5_3_risk_country, ncells=M5_3_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
M5_3_world_threats$cellstosize <- as.numeric((M5_3_world_threats$ncells*2500)/(M5_3_world_threats$AREA*10))
as.data.frame(M5_3_world_threats[order(M5_3_world_threats$cellstosize),c("NAME","cellstosize")])
M5_3_world_threats <- M5_3_world_threats[!is.na(M5_3_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
M5_3_world_threats <- M5_3_world_threats[M5_3_world_threats$cellstosize>=0.10,]

M5_3_world_complete <- merge(M5_3_world_threats, M5_3_counts[,c("country","n")], by.x="NAME", by.y="country",all=TRUE)
M5_3_world_complete$n[is.na(M5_3_world_complete$n)] <- 0
M5_3_world_complete <- M5_3_world_complete[!is.na(M5_3_world_complete$threats),]
unique(M5_3_world_complete$NAME)

M5_3_world_complete$nstudies <- M5_3_world_complete$n
#bivariate maps
M5_3_world_complete$threats_bin <- cut(M5_3_world_complete$risk, breaks = c(0, 3, 12, 35, 88), include.lowest = TRUE)
M5_3_world_complete$nstudies_bin <- cut(M5_3_world_complete$nstudies, breaks = c(0, 0.5, 5, 15, 30), include.lowest = TRUE)
classes <- bi_class(M5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(M5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_M5_3 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
ggsave("M5_3_risk.svg", height = 10, width = 12, bg = "white",device="svg")


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

# use negative binomial that accounts for zero inflation instead as poisson and quasi-poisson not suitable
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


#calculate the VIF for each predictor variable in the final models
vifmodel <- glmmTMB(nstudies ~ Huntingscaled + Pollutionscaled + Agriculturescaled + Loggingscaled + EPI_newscaled + GDP_2020scaled + NBIscaled + WGI_GE_Rankscaled + literacy_ratescaled, family=nbinom2, ziformula = ~0, data = M_LM_data)
check_collinearity(
  vifmodel,
  component = c('all') # 'all' shows both conditional and zi components
) #all non conservation risk below 5, low-moderate multicollinearity, not severe, except for agriculture and logging.
#agriculture is of low relative importance so remove.

vifmodel <- glmmTMB(nstudies ~ Huntingscaled + Pollutionscaled + Loggingscaled + EPI_newscaled + GDP_2020scaled + NBIscaled + WGI_GE_Rankscaled + literacy_ratescaled, family=nbinom2, ziformula = ~0, data = M_LM_data)
check_collinearity(
  vifmodel,
  component = c('all') # 'all' shows both conditional and zi components
) #all non conservation risk below 5, low-moderate multicollinearity, not severe.


# use negative binomial that accounts for zeron inflation instead as poisson and quasi-poisson not suitable
glm_Mnb <- glmmTMB(nstudies ~ Huntingscaled + Pollutionscaled + Loggingscaled + EPI_newscaled + GDP_2020scaled + NBIscaled + WGI_GE_Rankscaled + literacy_ratescaled, family=nbinom2, ziformula = ~0, data = M_LM_data)
summary(glm_Mnb)

# create global model and dredge to find best model
options(na.action = "na.fail") #Must run this code once to use dredge
all_models_m <- dredge(glm_Mnb)
#write.csv(data.frame(all_models_m),"all_models_m_withoutagri.csv")

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


# write.csv(summary(all_models_m_avg)[9],"glm_Mnb_summaryupdate_withoutagri.csv")
# write.csv(sw(all_models_m_avg),"glm_Mnb_relimp_withoutagri.csv")
# write.csv(confint(all_models_m_avg,full=TRUE),"glm_Mnb_confint_withoutagri.csv")

#robustness check
backward_model <- step(glm_Mnb, direction = "backward")














########## other threat maps ###########################################################



#### AMBHIBIANS AND CLIMATE CHANGE - A11 ####

## getting threats specific to amphibians and climate change 
A11_threats <- threats_data[!is.na(threats_data$A11),]

# applying quality filter to threats data 
A11_threats <- A11_threats %>% filter(ASR > 10)

# calculating risk
A11_threats$A11_risk <- A11_threats$A11 * A11_threats$ASR

## getting studies specific to amphibians and climate change
A_studies <- studies_data[studies_data$taxa == "Amphibians",]
A11_studies <- A_studies[A_studies$threattype == "Climate change & severe weather",]

## similarly for non-english studies - NO STUDIES 
#A11_non_english <- filter(non_english, taxa == "Amphibians" & threattype == "Climate change & severe weather")

## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
A11_studies <- subset(A11_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
A11_studies <- st_set_geometry(A11_studies, NULL)
A11_studies <- A11_studies[!is.na(A11_studies$country),]
A11_studies <- distinct(A11_studies)
A11_counts <- count(A11_studies, country)
#A11_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies - NO STUDIES
#A11_non_english_counts <- count(A11_non_english, country)

## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
A11_centroids <- st_centroid(A11_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean A11 threat value for each country
A11_centroids_sp <- as(A11_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
A11_country <- over(world_sp, A11_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!
A11_country_ncells <- over(world_sp, A11_centroids_sp, fn = length)
A11_country$ncells <- A11_country_ncells$A11_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
A11_risk_country <- A11_country$A11_risk
A11_threat_country <- A11_country$A11
A11_threat_ncells <- A11_country$ncells
A11_world_threats <- mutate(world, threats = A11_threat_country, risk = A11_risk_country, ncells=A11_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
A11_world_threats$cellstosize <- as.numeric((A11_world_threats$ncells*2500)/(A11_world_threats$AREA*10))
A11_world_threats[order(A11_world_threats$cellstosize),c("NAME","cellstosize")]
A11_world_threats <- A11_world_threats[!is.na(A11_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
A11_world_threats <- A11_world_threats[A11_world_threats$cellstosize>=0.10,]

A11_world_complete <- merge(A11_world_threats, A11_counts[,c("country","n")], by.x="NAME",by.y="country",all=TRUE)
A11_world_complete$n[is.na(A11_world_complete$n)] <- 0
A11_world_complete <- A11_world_complete[!is.na(A11_world_complete$threats),]
unique(A11_world_complete$NAME)

A11_world_complete$nstudies <- A11_world_complete$n
#bivariate maps
A11_world_complete$threats_bin <- cut(A11_world_complete$risk, breaks = c(0, 0.2, 0.5, 1, 3), include.lowest = TRUE)
A11_world_complete$nstudies_bin <- cut(A11_world_complete$n, breaks = c(0, 0.5, 1, 3, 6), include.lowest = TRUE)
classes <- bi_class(A11_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(A11_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A11 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("A11_risk.svg", height = 10, width = 12, bg = "white",device="svg")




#### AMBHIBIANS AND INVASIVES - A8_1 ####

## getting threats specific to amphibians and invasives
A8_1_threats <- threats_data[!is.na(threats_data$A8_1),]

# applying quality filter to threats data 
A8_1_threats <- A8_1_threats %>% filter(ASR > 10)

# calculating risk
A8_1_threats$A8_1_risk <- A8_1_threats$A8_1 * A8_1_threats$ASR

## getting studies specific to amphibians and invasives 
A_studies <- studies_data[studies_data$taxa == "Amphibians",]
A8_1_studies <- A_studies[!is.na(A_studies$threatname),]
A8_1_studies <- A8_1_studies[A8_1_studies$threatname == "Invasive non-native species",]

## similarly for non-english studies
A8_1_non_english <- filter(non_english, taxa == "Amphibians" & threatname == "Invasive non-native species")


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
A8_1_studies <- subset(A8_1_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
A8_1_studies <- st_set_geometry(A8_1_studies, NULL)
A8_1_studies <- A8_1_studies[!is.na(A8_1_studies$country),]
A8_1_studies <- distinct(A8_1_studies)
A8_1_counts <- count(A8_1_studies, country)
#A8_1_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
A8_1_non_english_counts <- count(A8_1_non_english, country)

#adding the two to get a complete count of studies per country 
A8_1_counts <- merge(A8_1_counts, A8_1_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
A8_1_counts[is.na(A8_1_counts)] <- 0
A8_1_counts$n <- A8_1_counts$n.x + A8_1_counts$n.y
A8_1_counts <- subset(A8_1_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
A8_1_centroids <- st_centroid(A8_1_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean A8_1 threat value for each country
A8_1_centroids_sp <- as(A8_1_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
A8_1_country <- over(world_sp, A8_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!
A8_1_country_ncells <- over(world_sp, A8_1_centroids_sp, fn = length)
A8_1_country$ncells <- A8_1_country_ncells$A8_1_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
A8_1_risk_country <- A8_1_country$A8_1_risk
A8_1_threat_country <- A8_1_country$A8_1
A8_1_threat_ncells <- A8_1_country$ncells
A8_1_world_threats <- mutate(world, threats = A8_1_threat_country, risk = A8_1_risk_country, ncells=A8_1_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
A8_1_world_threats$cellstosize <- as.numeric((A8_1_world_threats$ncells*2500)/(A8_1_world_threats$AREA*10))
A8_1_world_threats[order(A8_1_world_threats$cellstosize),c("NAME","cellstosize")]
A8_1_world_threats <- A8_1_world_threats[!is.na(A8_1_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
A8_1_world_threats <- A8_1_world_threats[A8_1_world_threats$cellstosize>=0.10,]

A8_1_world_complete <- merge(A8_1_world_threats, A8_1_counts[,c("country","n")], by.x="NAME",by.y="country",all=TRUE)
A8_1_world_complete$n[is.na(A8_1_world_complete$n)] <- 0
A8_1_world_complete <- A8_1_world_complete[!is.na(A8_1_world_complete$threats),]
unique(A8_1_world_complete$NAME)

A8_1_world_complete$nstudies <- A8_1_world_complete$n
#bivariate maps
A8_1_world_complete$threats_bin <- cut(A8_1_world_complete$risk, breaks = c(0, 1, 5, 15, 31), include.lowest = TRUE)
A8_1_world_complete$nstudies_bin <- cut(A8_1_world_complete$n, breaks = c(0, 0.5, 1, 3, 6), include.lowest = TRUE)
classes <- bi_class(A8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(A8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A8_1 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("A8_1_risk.svg", height = 10, width = 12, bg = "white",device="svg")




#### AMPHIBIANS AND HUNTING - A5_1 ####

## getting threats specific to amphibians and hunting
A5_1_threats <- threats_data[!is.na(threats_data$A5_1),]

# applying quality filter to threats data 
A5_1_threats <- A5_1_threats %>% filter(ASR > 10)

# calculating risk
A5_1_threats$A5_1_risk <- A5_1_threats$A5_1 * A5_1_threats$ASR

## getting studies specific to amphibians and hunting
A_studies <- studies_data[studies_data$taxa == "Amphibians",]
A5_1_studies <- A_studies[!is.na(A_studies$threatname),]
A5_1_studies <- A5_1_studies[A5_1_studies$threatname == "Hunting/trapping terrestrial animals",]

## similarly for non-english studies
A5_1_non_english <- filter(non_english, taxa == "Amphibians" & threatname == "Hunting/trapping terrestrial animals")


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
A5_1_studies <- subset(A5_1_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
A5_1_studies <- st_set_geometry(A5_1_studies, NULL)
A5_1_studies <- A5_1_studies[!is.na(A5_1_studies$country),]
A5_1_studies <- distinct(A5_1_studies)
A5_1_counts <- count(A5_1_studies, country)
#A5_1_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
A5_1_non_english_counts <- count(A5_1_non_english, country)

#adding the two to get a complete count of studies per country 
A5_1_counts <- merge(A5_1_counts, A5_1_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
A5_1_counts[is.na(A5_1_counts)] <- 0
A5_1_counts$n <- A5_1_counts$n.x + A5_1_counts$n.y
A5_1_counts <- subset(A5_1_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
A5_1_centroids <- st_centroid(A5_1_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean A5_1 threat value for each country
A5_1_centroids_sp <- as(A5_1_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
A5_1_country <- over(world_sp, A5_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!
A5_1_country_ncells <- over(world_sp, A5_1_centroids_sp, fn = length)
A5_1_country$ncells <- A5_1_country_ncells$A5_1_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
A5_1_risk_country <- A5_1_country$A5_1_risk
A5_1_threat_country <- A5_1_country$A5_1
A5_1_threat_ncells <- A5_1_country$ncells
A5_1_world_threats <- mutate(world, threats = A5_1_threat_country, risk = A5_1_risk_country, ncells=A5_1_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
A5_1_world_threats$cellstosize <- as.numeric((A5_1_world_threats$ncells*2500)/(A5_1_world_threats$AREA*10))
A5_1_world_threats[order(A5_1_world_threats$cellstosize),c("NAME","cellstosize")]
A5_1_world_threats <- A5_1_world_threats[!is.na(A5_1_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
A5_1_world_threats <- A5_1_world_threats[A5_1_world_threats$cellstosize>=0.10,]

A5_1_world_complete <- merge(A5_1_world_threats, A5_1_counts[,c("country","n")], by.x="NAME",by.y="country",all=TRUE)
A5_1_world_complete$n[is.na(A5_1_world_complete$n)] <- 0
A5_1_world_complete <- A5_1_world_complete[!is.na(A5_1_world_complete$threats),]
unique(A5_1_world_complete$NAME)

A5_1_world_complete$nstudies <- A5_1_world_complete$n
#bivariate maps
A5_1_world_complete$threats_bin <- cut(A5_1_world_complete$risk, breaks = c(0, 1, 3, 9), include.lowest = TRUE)
A5_1_world_complete$nstudies_bin <- cut(A5_1_world_complete$n, breaks = c(0, 0.5, 1, 2), include.lowest = TRUE)
classes <- bi_class(A5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 3)
breaks1 <- bi_class_breaks(A5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 3, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 3) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 3, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A5_1 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("A5_1_risk.svg", height = 10, width = 12, bg = "white",device="svg")


#### AMPHIBIANS AND LOGGING - A5_3 ####

## getting threats specific to amphibians and logging
A5_3_threats <- threats_data[!is.na(threats_data$A5_3),]

# applying quality filter to threats data 
A5_3_threats <- A5_3_threats %>% filter(ASR > 10)

# calculating risk
A5_3_threats$A5_3_risk <- A5_3_threats$A5_3 * A5_3_threats$ASR

## getting studies specific to amphibians and logging
A_studies <- studies_data[studies_data$taxa == "Amphibians",]
A5_3_studies <- A_studies[!is.na(A_studies$threatname),]
A5_3_studies <- A5_3_studies[A5_3_studies$threatname == "Logging/wood harvesting",]

## similarly for non-english studies - NO STUDIES
#A5_3_non_english <- filter(non_english, taxa == "Amphibians" & threatname == "Logging/wood harvesting")

## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
A5_3_studies <- subset(A5_3_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
A5_3_studies <- st_set_geometry(A5_3_studies, NULL)
A5_3_studies <- A5_3_studies[!is.na(A5_3_studies$country),]
A5_3_studies <- distinct(A5_3_studies)
A5_3_counts <- count(A5_3_studies, country)
#A5_3_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies - NO STUDIES
#A5_3_non_english_counts <- count(A5_3_non_english, country)


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
A5_3_centroids <- st_centroid(A5_3_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean A5_3 threat value for each country
A5_3_centroids_sp <- as(A5_3_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
A5_3_country <- over(world_sp, A5_3_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!
A5_3_country_ncells <- over(world_sp, A5_3_centroids_sp, fn = length)
A5_3_country$ncells <- A5_3_country_ncells$A5_3_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
A5_3_risk_country <- A5_3_country$A5_3_risk
A5_3_threat_country <- A5_3_country$A5_3
A5_3_threat_ncells <- A5_3_country$ncells
A5_3_world_threats <- mutate(world, threats = A5_3_threat_country, risk = A5_3_risk_country, ncells=A5_3_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
A5_3_world_threats$cellstosize <- as.numeric((A5_3_world_threats$ncells*2500)/(A5_3_world_threats$AREA*10))
A5_3_world_threats[order(A5_3_world_threats$cellstosize),c("NAME","cellstosize")]
A5_3_world_threats <- A5_3_world_threats[!is.na(A5_3_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
A5_3_world_threats <- A5_3_world_threats[A5_3_world_threats$cellstosize>=0.10,]

A5_3_world_complete <- merge(A5_3_world_threats, A5_3_counts[,c("country","n")], by.x="NAME",by.y="country",all=TRUE)
A5_3_world_complete$n[is.na(A5_3_world_complete$n)] <- 0
A5_3_world_complete <- A5_3_world_complete[!is.na(A5_3_world_complete$threats),]
unique(A5_3_world_complete$NAME)

A5_3_world_complete$nstudies <- A5_3_world_complete$n
#bivariate maps
A5_3_world_complete$threats_bin <- cut(A5_3_world_complete$risk, breaks = c(0, 2, 6, 15, 39), include.lowest = TRUE)
A5_3_world_complete$nstudies_bin <- cut(A5_3_world_complete$n, breaks = c(0, 0.5, 5, 15, 26), include.lowest = TRUE)
classes <- bi_class(A5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(A5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A5_3 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("A5_3_risk.svg", height = 10, width = 12, bg = "white",device="svg")




#### BIRDS AND CLIMATE CHANGE - B11 ####

## getting threats specific to birds and climate change 
B11_threats <- threats_data[!is.na(threats_data$B11),]

# applying quality filter to threats data 
B11_threats <- B11_threats %>% filter(BSR > 10)

# calculating risk
B11_threats$B11_risk <- B11_threats$B11 * B11_threats$BSR

## getting studies specific to birds and climate change 
B_studies <- studies_data[studies_data$taxa == "Birds",]
B11_studies <- B_studies[B_studies$threattype == "Climate change & severe weather",]

## similarly for non-english studies - NO STUDIES
#B11_non_english <- filter(non_english, taxa == "Birds" & threattype == "Climate change & severe weather")

## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
B11_studies <- subset(B11_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
B11_studies <- st_set_geometry(B11_studies, NULL)
B11_studies <- B11_studies[!is.na(B11_studies$country),]
B11_studies <- distinct(B11_studies)
B11_counts <- count(B11_studies, country)
#B11_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies - NO STUDIES
#B11_non_english_counts <- count(B11_non_english, country)


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
B11_centroids <- st_centroid(B11_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean B11 threat value for each country
B11_centroids_sp <- as(B11_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
B11_country <- over(world_sp, B11_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!
B11_country_ncells <- over(world_sp, B11_centroids_sp, fn = length)
B11_country$ncells <- B11_country_ncells$B11_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
B11_risk_country <- B11_country$B11_risk
B11_threat_country <- B11_country$B11
B11_threat_ncells <- B11_country$ncells
B11_world_threats <- mutate(world, threats = B11_threat_country, risk = B11_risk_country, ncells=B11_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
B11_world_threats$cellstosize <- as.numeric((B11_world_threats$ncells*2500)/(B11_world_threats$AREA*10))
B11_world_threats[order(B11_world_threats$cellstosize),c("NAME","cellstosize")]
B11_world_threats <- B11_world_threats[!is.na(B11_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
B11_world_threats <- B11_world_threats[B11_world_threats$cellstosize>=0.10,]

B11_world_complete <- merge(B11_world_threats, B11_counts[,c("country","n")], by.x="NAME",by.y="country",all=TRUE)
B11_world_complete$n[is.na(B11_world_complete$n)] <- 0
B11_world_complete <- B11_world_complete[!is.na(B11_world_complete$threats),]
unique(B11_world_complete$NAME)

B11_world_complete$nstudies <- B11_world_complete$n
#bivariate maps
B11_world_complete$threats_bin <- cut(B11_world_complete$risk, breaks = c(0, 5, 20, 52), include.lowest = TRUE)
B11_world_complete$nstudies_bin <- cut(B11_world_complete$n, breaks = c(0, 0.5, 1, 2), include.lowest = TRUE)
classes <- bi_class(B11_world_complete, x = threats_bin, y = nstudies_bin, dim = 3)
breaks1 <- bi_class_breaks(B11_world_complete, x = threats_bin, y = nstudies_bin, dim = 3, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 3) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 3, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_B11 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("B11_risk.svg", height = 10, width = 12, bg = "white",device="svg")



#### BIRDS AND HUNTING - B5_1 ####

## getting threats specific to birds and hunting
B5_1_threats <- threats_data[!is.na(threats_data$B5_1),]

# applying quality filter to threats data 
B5_1_threats <- B5_1_threats %>% filter(BSR > 10)

# calculating risk
B5_1_threats$B5_1_risk <- B5_1_threats$B5_1 * B5_1_threats$BSR

## getting studies specific to birds and hunting
B_studies <- studies_data[studies_data$taxa == "Birds",]
B5_1_studies <- B_studies[!is.na(B_studies$threatname),]
B5_1_studies <- B5_1_studies[B5_1_studies$threatname == "Hunting/trapping terrestrial animals",]

## similarly for non-english studies
B5_1_non_english <- filter(non_english, taxa == "Birds" & threatname == "Hunting/trapping terrestrial animals")

## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
B5_1_studies <- subset(B5_1_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
B5_1_studies <- st_set_geometry(B5_1_studies, NULL)
B5_1_studies <- B5_1_studies[!is.na(B5_1_studies$country),]
B5_1_studies <- distinct(B5_1_studies)
B5_1_counts <- count(B5_1_studies, country)
#B5_1_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
B5_1_non_english_counts <- count(B5_1_non_english, country)

#adding the two to get a complete count of studies per country 
B5_1_counts <- merge(B5_1_counts, B5_1_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
B5_1_counts[is.na(B5_1_counts)] <- 0
B5_1_counts$n <- B5_1_counts$n.x + B5_1_counts$n.y
B5_1_counts <- subset(B5_1_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
B5_1_centroids <- st_centroid(B5_1_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean B5_1 threat value for each country
B5_1_centroids_sp <- as(B5_1_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
B5_1_country <- over(world_sp, B5_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!
B5_1_country_ncells <- over(world_sp, B5_1_centroids_sp, fn = length)
B5_1_country$ncells <- B5_1_country_ncells$B5_1_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
B5_1_risk_country <- B5_1_country$B5_1_risk
B5_1_threat_country <- B5_1_country$B5_1
B5_1_threat_ncells <- B5_1_country$ncells
B5_1_world_threats <- mutate(world, threats = B5_1_threat_country, risk = B5_1_risk_country, ncells=B5_1_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
B5_1_world_threats$cellstosize <- as.numeric((B5_1_world_threats$ncells*2500)/(B5_1_world_threats$AREA*10))
B5_1_world_threats[order(B5_1_world_threats$cellstosize),c("NAME","cellstosize")]
B5_1_world_threats <- B5_1_world_threats[!is.na(B5_1_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
B5_1_world_threats <- B5_1_world_threats[B5_1_world_threats$cellstosize>=0.10,]

B5_1_world_complete <- merge(B5_1_world_threats, B5_1_counts[,c("country","n")], by.x="NAME",by.y="country",all=TRUE)
B5_1_world_complete$n[is.na(B5_1_world_complete$n)] <- 0
B5_1_world_complete <- B5_1_world_complete[!is.na(B5_1_world_complete$threats),]
unique(B5_1_world_complete$NAME)

B5_1_world_complete$nstudies <- B5_1_world_complete$n
#bivariate maps
B5_1_world_complete$threats_bin <- cut(B5_1_world_complete$risk, breaks = c(0, 5, 15, 30, 46), include.lowest = TRUE)
B5_1_world_complete$nstudies_bin <- cut(B5_1_world_complete$n, breaks = c(0, 0.5, 1, 2, 3), include.lowest = TRUE)
classes <- bi_class(B5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(B5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_B5_1 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
ggsave("B5_1_risk.svg", height = 10, width = 12, bg = "white",device="svg")




#### MAMMALS AND CLIMATE CHANGE - M11 ####

## getting threats specific to mammals and climate change 
M11_threats <- threats_data[!is.na(threats_data$M11),]

# applying quality filter to threats data 
M11_threats <- M11_threats %>% filter(MSR > 10)

# calculating risk 
M11_threats$M11_risk <- M11_threats$M11 * M11_threats$MSR

## getting studies specific to mammals and climate change 
M_studies <- studies_data[studies_data$taxa == "Mammals",]
M11_studies <- M_studies[M_studies$threattype == "Climate change & severe weather",]


## similarly for non-english studies - NO STUDIES
#M11_non_english <- filter(non_english, taxa == "Mammals" & threattype == "Climate change & severe weather")

## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
M11_studies <- subset(M11_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
M11_studies <- st_set_geometry(M11_studies, NULL)
M11_studies <- M11_studies[!is.na(M11_studies$country),]
M11_studies <- distinct(M11_studies)
M11_counts <- count(M11_studies, country)
#M11_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies - NO STUDIES
#M11_non_english_counts <- count(M11_non_english, country)


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
M11_centroids <- st_centroid(M11_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean M11 threat value for each country
M11_centroids_sp <- as(M11_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
M11_country <- over(world_sp, M11_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!
M11_country_ncells <- over(world_sp, M11_centroids_sp, fn = length)
M11_country$ncells <- M11_country_ncells$M11_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
M11_risk_country <- M11_country$M11_risk
M11_threat_country <- M11_country$M11
M11_threat_ncells <- M11_country$ncells
M11_world_threats <- mutate(world, threats = M11_threat_country, risk = M11_risk_country, ncells=M11_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
M11_world_threats$cellstosize <- as.numeric((M11_world_threats$ncells*2500)/(M11_world_threats$AREA*10))
M11_world_threats[order(M11_world_threats$cellstosize),c("NAME","cellstosize")]
M11_world_threats <- M11_world_threats[!is.na(M11_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
M11_world_threats <- M11_world_threats[M11_world_threats$cellstosize>=0.10,]

M11_world_complete <- merge(M11_world_threats, M11_counts[,c("country","n")], by.x="NAME",by.y="country",all=TRUE)
M11_world_complete$n[is.na(M11_world_complete$n)] <- 0
M11_world_complete <- M11_world_complete[!is.na(M11_world_complete$threats),]
unique(M11_world_complete$NAME)

M11_world_complete$nstudies <- M11_world_complete$n
#bivariate maps
M11_world_complete$threats_bin <- cut(M11_world_complete$risk, breaks = c(0, 0.2, 1, 3, 6), include.lowest = TRUE)
M11_world_complete$nstudies_bin <- cut(M11_world_complete$n, breaks = c(0, 0.5, 1, 3, 5), include.lowest = TRUE)
classes <- bi_class(M11_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(M11_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_M11 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
#ggsave("M11_risk.svg", height = 10, width = 12, bg = "white",device="svg")



#### MAMMALS AND INVASIVES - M8_1 ####

## getting threats specific to mammals and invasives
M8_1_threats <- threats_data[!is.na(threats_data$M8_1),]

# applying quality filter to threats data 
M8_1_threats <- M8_1_threats %>% filter(MSR > 10)

# calculating risk
M8_1_threats$M8_1_risk <- M8_1_threats$M8_1 * M8_1_threats$MSR

## getting studies specific to mammals and invasives
M_studies <- studies_data[studies_data$taxa == "Mammals",]
M8_1_studies <- M_studies[!is.na(M_studies$threatname),]
M8_1_studies <- M8_1_studies[M8_1_studies$threatname == "Invasive non-native species",]

## similarly for non-english studies
M8_1_non_english <- filter(non_english, taxa == "Mammals" & threatname == "Invasive non-native species")

## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
M8_1_studies <- subset(M8_1_studies, select = -c(region, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
M8_1_studies <- st_set_geometry(M8_1_studies, NULL)
M8_1_studies <- M8_1_studies[!is.na(M8_1_studies$country),]
M8_1_studies <- distinct(M8_1_studies)
M8_1_counts <- count(M8_1_studies, country)
#M8_1_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 

#for non-english studies
M8_1_non_english_counts <- count(M8_1_non_english, country)

#adding the two to get a complete count of studies per country 
M8_1_counts <- merge(M8_1_counts, M8_1_non_english_counts, by = "country", all.x = TRUE, all.y = TRUE)
M8_1_counts[is.na(M8_1_counts)] <- 0
M8_1_counts$n <- M8_1_counts$n.x + M8_1_counts$n.y
M8_1_counts <- subset(M8_1_counts, select = -c(n.x, n.y))


## now let's look at threat intensity per country

#converting our polygon threats dataframe to a spatial points dataframe 
M8_1_centroids <- st_centroid(M8_1_threats)

#getting a sf object of world countries
data("wrld_simpl")
world <- st_as_sf(wrld_simpl)
world <- st_transform(world, crs = st_crs(threats_data))

#now using the world dataframe and the threat points dataframe to generate a mean M8_1 threat value for each country
M8_1_centroids_sp <- as(M8_1_centroids, Class = "Spatial")
world_sp <- as(world, Class = "Spatial")
M8_1_country <- over(world_sp, M8_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!
M8_1_country_ncells <- over(world_sp, M8_1_centroids_sp, fn = length)
M8_1_country$ncells <- M8_1_country_ncells$M8_1_risk

## creating a complete dataframe (this contains all countries for which there is threat data)
M8_1_risk_country <- M8_1_country$M8_1_risk
M8_1_threat_country <- M8_1_country$M8_1
M8_1_threat_ncells <- M8_1_country$ncells
M8_1_world_threats <- mutate(world, threats = M8_1_threat_country, risk = M8_1_risk_country, ncells=M8_1_threat_ncells)

#need to multiply area by 10 to get km2, cells are 50x50km so 2500km2
M8_1_world_threats$cellstosize <- as.numeric((M8_1_world_threats$ncells*2500)/(M8_1_world_threats$AREA*10))
M8_1_world_threats[order(M8_1_world_threats$cellstosize),c("NAME","cellstosize")]
M8_1_world_threats <- M8_1_world_threats[!is.na(M8_1_world_threats$threats),]
#remove countries with less than 10% of their area covered by cells
#remove, e.g., Russia, Iran, Sweden etc.
M8_1_world_threats <- M8_1_world_threats[M8_1_world_threats$cellstosize>=0.10,]

M8_1_world_complete <- merge(M8_1_world_threats, M8_1_counts[,c("country","n")], by.x="NAME",by.y="country",all=TRUE)
M8_1_world_complete$n[is.na(M8_1_world_complete$n)] <- 0
M8_1_world_complete <- M8_1_world_complete[!is.na(M8_1_world_complete$threats),]
unique(M8_1_world_complete$NAME)

M8_1_world_complete$nstudies <- M8_1_world_complete$n
#bivariate maps
M8_1_world_complete$threats_bin <- cut(M8_1_world_complete$risk, breaks = c(0, 2, 5, 13), include.lowest = TRUE)
M8_1_world_complete$nstudies_bin <- cut(M8_1_world_complete$n, breaks = c(0, 0.5, 1, 2), include.lowest = TRUE)
classes <- bi_class(M8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 3)
breaks1 <- bi_class_breaks(M8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 3, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + 
  geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + 
  bi_scale_fill(pal = "GrPink2", dim = 3) + 
  bi_theme()

legend <- bi_legend(pal = "GrPink2", dim = 3, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_M8_1 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
ggsave("M8_1_risk.svg", height = 10, width = 12, bg = "white",device="svg")


