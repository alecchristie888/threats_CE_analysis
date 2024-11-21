#loading necessary libraries
library(sf)
library(ggplot2)
library(dplyr)
library(maptools)
library(biscale)
library(cowplot)
library(ggpattern)
library(ggrepel)
library(stats)
library(readxl)


#set working directory
#setwd("")


#reading in files
studies_data <- read.csv("Data/CE_data_complete.csv")
non_english <- read_excel("Data/CE_non_english_data_clean.xlsx")
threats_data <- st_read("Data/Grid_mammal_amph_threat_predictions.shp")

#cleaning up
studies_data <- studies_data[!is.na(studies_data$long),]
studies_data <- filter(studies_data, long >= -180 & long <= 180)
studies_data <- filter(studies_data, lat >= -90 & lat <= 90)
studies_data <- filter(studies_data, !(lat == 0 & long == 0))

#converting studies data to a shapefile
studies_data <- st_as_sf(studies_data, coords = c("long", "lat"), crs = 4326)

#setting studies data to same crs as threats data 
studies_data <- st_transform(studies_data, crs = st_crs(threats_data))


#### AMBHIBIANS AND AGRICULTURE - A2 ####

## getting threats specific to amphibians and agriculture 
A2_threats <- threats_data[!is.na(threats_data$A2),]

# applying quality filter to threats data 
A2_threats <- A2_threats %>% filter(ASR > 10)

# calculating risk 
A2_threats$A2 <- A2_threats$A2 * A2_threats$ASR

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
A2_threat_country <- over(world_sp, A2_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
A2_threat_country <- A2_threat_country$A2
A2_world_threats <- mutate(world, threats = A2_threat_country)
A2_world_threats <- A2_world_threats[order(A2_world_threats$NAME),]#ordering alphabetically according to country name
A2_world_threats_only <- A2_world_threats[!is.na(A2_world_threats$threats),]
summary(A2_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = A2_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                                                       legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
A2_threats_rank <- subset(A2_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(A2_threats_rank) <- NULL
A2_threats_rank <- A2_threats_rank[order(-A2_threats_rank$threats),]

##creating a complete dataframe 
A2_counts <- rename(A2_counts, NAME = country)
A2_world_studies <- merge(world, A2_counts, by = "NAME", all = TRUE)

A2_world_complete <- mutate(A2_world_threats, nstudies = A2_world_studies$n)
A2_world_complete <- A2_world_complete[!is.na(A2_world_complete$threats),]
A2_world_complete$nstudies[is.na(A2_world_complete$nstudies)] <- 0

## correlation between threats and studies 
cor.test(A2_world_complete$threats, A2_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
A2_world_complete$threats_bin <- cut(A2_world_complete$threats, breaks = c(0, 5, 10, 15, 38), include.lowest = TRUE)
A2_world_complete$nstudies_bin <- cut(A2_world_complete$nstudies, breaks = c(0, 0.5, 6, 15, 30), include.lowest = TRUE)
classes <- bi_class(A2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(A2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A2 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot_A2



#### AMBHIBIANS AND POLLUTION - A9 ####

## getting threats specific to amphibians and pollution 
A9_threats <- threats_data[!is.na(threats_data$A9),]

# applying quality filter to threats data 
A9_threats <- A9_threats %>% filter(ASR > 10)

# calculating risk
A9_threats$A9 <- A9_threats$A9 * A9_threats$ASR

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
A9_threat_country <- over(world_sp, A9_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
A9_threat_country <- A9_threat_country$A9
A9_world_threats <- mutate(world, threats = A9_threat_country)
A9_world_threats <- A9_world_threats[order(A9_world_threats$NAME),]#ordering alphabetically according to country name
A9_world_threats_only <- A9_world_threats[!is.na(A9_world_threats$threats),]
summary(A9_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = A9_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
A9_threats_rank <- subset(A9_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(A9_threats_rank) <- NULL
A9_threats_rank <- A9_threats_rank[order(-A9_threats_rank$threats),]


## creating a complete dataframe (this contains all countries for which there is threat data)
A9_counts <- rename(A9_counts, NAME = country)
A9_world_studies <- merge(world, A9_counts, by = "NAME", all = TRUE)

A9_world_complete <- mutate(A9_world_threats, nstudies = A9_world_studies$n)
A9_world_complete <- A9_world_complete[!is.na(A9_world_complete$threats),]
A9_world_complete$nstudies[is.na(A9_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(A9_world_complete$threats, A9_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
A9_world_complete$threats_bin <- cut(A9_world_complete$threats, breaks = c(0, 0.50, 1, 3, 7), include.lowest = TRUE)
A9_world_complete$nstudies_bin <- cut(A9_world_complete$nstudies, breaks = c(0, 0.5, 1, 5, 9), include.lowest = TRUE)
classes <- bi_class(A9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(A9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A9 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot_A9


#### AMBHIBIANS AND CLIMATE CHANGE - A11 ####

## getting threats specific to amphibians and climate change 
A11_threats <- threats_data[!is.na(threats_data$A11),]

# applying quality filter to threats data 
A11_threats <- A11_threats %>% filter(ASR > 10)

# calculating risk
A11_threats$A11 <- A11_threats$A11 * A11_threats$ASR

## getting studies specific to amphibians and climate change
A_studies <- studies_data[studies_data$taxa == "Amphibians",]
A11_studies <- A_studies[A_studies$threattype == "Climate change & severe weather",]


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
A11_threat_country <- over(world_sp, A11_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
A11_threat_country <- A11_threat_country$A11
A11_world_threats <- mutate(world, threats = A11_threat_country)
A11_world_threats <- A11_world_threats[order(A11_world_threats$NAME),]#ordering alphabetically according to country name
A11_world_threats_only <- A11_world_threats[!is.na(A11_world_threats$threats),]
summary(A11_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = A11_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
A11_threats_rank <- subset(A11_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(A11_threats_rank) <- NULL
A11_threats_rank <- A11_threats_rank[order(-A11_threats_rank$threats),]


## creating a complete dataframe (this contains all countries for which there is threat data)
A11_counts <- rename(A11_counts, NAME = country)
A11_world_studies <- merge(world, A11_counts, by = "NAME", all = TRUE)

A11_world_complete <- mutate(A11_world_threats, nstudies = A11_world_studies$n)
A11_world_complete <- A11_world_complete[!is.na(A11_world_complete$threats),]
A11_world_complete$nstudies[is.na(A11_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(A11_world_complete$threats, A11_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
A11_world_complete$threats_bin <- cut(A11_world_complete$threats, breaks = c(0, 0.2, 0.5, 1, 3), include.lowest = TRUE)
A11_world_complete$nstudies_bin <- cut(A11_world_complete$nstudies, breaks = c(0, 0.5, 1, 3, 6), include.lowest = TRUE)
classes <- bi_class(A11_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(A11_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A11 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot_A9


#### AMBHIBIANS AND INVASIVES - A8_1 ####

## getting threats specific to amphibians and invasives
A8_1_threats <- threats_data[!is.na(threats_data$A8_1),]

# applying quality filter to threats data 
A8_1_threats <- A8_1_threats %>% filter(ASR > 10)

# calculating risk
A8_1_threats$A8_1 <- A8_1_threats$A8_1 * A8_1_threats$ASR

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
A8_1_studies <- subset(A8_1_studies, select = -c(region, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
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
A8_1_threat_country <- over(world_sp, A8_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
A8_1_threat_country <- A8_1_threat_country$A8_1
A8_1_world_threats <- mutate(world, threats = A8_1_threat_country)
A8_1_world_threats <- A8_1_world_threats[order(A8_1_world_threats$NAME),]#ordering alphabetically according to country name
A8_1_world_threats_only <- A8_1_world_threats[!is.na(A8_1_world_threats$threats),]
summary(A8_1_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = A8_1_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
A8_1_threats_rank <- subset(A8_1_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(A8_1_threats_rank) <- NULL
A8_1_threats_rank <- A8_1_threats_rank[order(-A8_1_threats_rank$threats),]

## creating a complete dataframe (this contains all countries for which there is threat data)
A8_1_counts <- rename(A8_1_counts, NAME = country)
A8_1_world_studies <- merge(world, A8_1_counts, by = "NAME", all = TRUE)

A8_1_world_complete <- mutate(A8_1_world_threats, nstudies = A8_1_world_studies$n)
A8_1_world_complete <- A8_1_world_complete[!is.na(A8_1_world_complete$threats),]
A8_1_world_complete$nstudies[is.na(A8_1_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(A8_1_world_complete$threats, A8_1_world_complete$nstudies, method = 'spearman')

## bivariate heatmap
A8_1_world_complete$threats_bin <- cut(A8_1_world_complete$threats, breaks = c(0, 1, 5, 15, 31), include.lowest = TRUE)
A8_1_world_complete$nstudies_bin <- cut(A8_1_world_complete$nstudies, breaks = c(0, 0.5, 1, 3, 6), include.lowest = TRUE)
classes <- bi_class(A8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(A8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A8_1 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot_A8_1



#### AMPHIBIANS AND HUNTING - A5_1 ####

## getting threats specific to amphibians and hunting
A5_1_threats <- threats_data[!is.na(threats_data$A5_1),]

# applying quality filter to threats data 
A5_1_threats <- A5_1_threats %>% filter(ASR > 10)

# calculating risk
A5_1_threats$A5_1 <- A5_1_threats$A5_1 * A5_1_threats$ASR

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
A5_1_studies <- subset(A5_1_studies, select = -c(region, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
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
A5_1_threat_country <- over(world_sp, A5_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
A5_1_threat_country <- A5_1_threat_country$A5_1
A5_1_world_threats <- mutate(world, threats = A5_1_threat_country)
A5_1_world_threats <- A5_1_world_threats[order(A5_1_world_threats$NAME),]#ordering alphabetically according to country name
A5_1_world_threats_only <- A5_1_world_threats[!is.na(A5_1_world_threats$threats),]
summary(A5_1_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = A5_1_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
A5_1_threats_rank <- subset(A5_1_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(A5_1_threats_rank) <- NULL
A5_1_threats_rank <- A5_1_threats_rank[order(-A5_1_threats_rank$threats),]


## creating a complete dataframe (this contains all countries for which there is threat data)
A5_1_counts <- rename(A5_1_counts, NAME = country)
A5_1_world_studies <- merge(world, A5_1_counts, by = "NAME", all = TRUE)

A5_1_world_complete <- mutate(A5_1_world_threats, nstudies = A5_1_world_studies$n)
A5_1_world_complete <- A5_1_world_complete[!is.na(A5_1_world_complete$threats),]
A5_1_world_complete$nstudies[is.na(A5_1_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(A5_1_world_complete$threats, A5_1_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
A5_1_world_complete$threats_bin <- cut(A5_1_world_complete$threats, breaks = c(0, 1, 3, 9), include.lowest = TRUE)
A5_1_world_complete$nstudies_bin <- cut(A5_1_world_complete$nstudies, breaks = c(0, 0.5, 1, 2), include.lowest = TRUE)
classes <- bi_class(A5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 3)
breaks1 <- bi_class_breaks(A5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 3, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink", dim = 3) + bi_theme()
legend <- bi_legend(pal = "GrPink", dim = 3, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A5_1 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot_A5_1



#### AMPHIBIANS AND LOGGING - A5_3 ####

## getting threats specific to amphibians and logging
A5_3_threats <- threats_data[!is.na(threats_data$A5_3),]

# applying quality filter to threats data 
A5_3_threats <- A5_3_threats %>% filter(ASR > 10)

# calculating risk
A5_3_threats$A5_3 <- A5_3_threats$A5_3 * A5_3_threats$ASR

## getting studies specific to amphibians and logging
A_studies <- studies_data[studies_data$taxa == "Amphibians",]
A5_3_studies <- A_studies[!is.na(A_studies$threatname),]
A5_3_studies <- A5_3_studies[A5_3_studies$threatname == "Logging/wood harvesting",]


## number of studies per country 

#first need to deal with multiple rows corresponding to the same study - 
#will drop the columns that vary even within studies so I can then get one row for each study, and then count
#according to country
A5_3_studies <- subset(A5_3_studies, select = -c(region, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
A5_3_studies <- st_set_geometry(A5_3_studies, NULL)
A5_3_studies <- A5_3_studies[!is.na(A5_3_studies$country),]
A5_3_studies <- distinct(A5_3_studies)
A5_3_counts <- count(A5_3_studies, country)
#A5_3_counts here functions as a good list of number of studies per country as some studies take place across multiple countries 


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
A5_3_threat_country <- over(world_sp, A5_3_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
A5_3_threat_country <- A5_3_threat_country$A5_3
A5_3_world_threats <- mutate(world, threats = A5_3_threat_country)
A5_3_world_threats <- A5_3_world_threats[order(A5_3_world_threats$NAME),]#ordering alphabetically according to country name
A5_3_world_threats_only <- A5_3_world_threats[!is.na(A5_3_world_threats$threats),]
summary(A5_3_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = A5_3_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
A5_3_threats_rank <- subset(A5_3_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(A5_3_threats_rank) <- NULL
A5_3_threats_rank <- A5_3_threats_rank[order(-A5_3_threats_rank$threats),]


## creating a complete dataframe (this contains all countries for which there is threat data)
A5_3_counts <- rename(A5_3_counts, NAME = country)
A5_3_world_studies <- merge(world, A5_3_counts, by = "NAME", all = TRUE)

A5_3_world_complete <- mutate(A5_3_world_threats, nstudies = A5_3_world_studies$n)
A5_3_world_complete <- A5_3_world_complete[!is.na(A5_3_world_complete$threats),]
A5_3_world_complete$nstudies[is.na(A5_3_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(A5_3_world_complete$threats, A5_3_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
A5_3_world_complete$threats_bin <- cut(A5_3_world_complete$threats, breaks = c(0, 2, 6, 15, 39), include.lowest = TRUE)
A5_3_world_complete$nstudies_bin <- cut(A5_3_world_complete$nstudies, breaks = c(0, 0.5, 5, 15, 26), include.lowest = TRUE)
classes <- bi_class(A5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(A5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot_A5_3 <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot_A5_3



#### BIRDS AND AGRICULTURE - B2 ####

## getting threats specific to birds and agriculture 
B2_threats <- threats_data[!is.na(threats_data$B2),]

# applying quality filter to threats data 
B2_threats <- B2_threats %>% filter(BSR > 10)

# calculating risk
B2_threats$B2 <- B2_threats$B2 * B2_threats$BSR

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
B2_threat_country <- over(world_sp, B2_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
B2_threat_country <- B2_threat_country$B2
B2_world_threats <- mutate(world, threats = B2_threat_country)
B2_world_threats <- B2_world_threats[order(B2_world_threats$NAME),]#ordering alphabetically according to country name
B2_world_threats_only <- B2_world_threats[!is.na(B2_world_threats$threats),]
summary(B2_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = B2_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')


## creating rank tables
B2_threats_rank <- subset(B2_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(B2_threats_rank) <- NULL
B2_threats_rank <- B2_threats_rank[order(-B2_threats_rank$threats),]


## creating a complete dataframe (this contains all countries for which there is threat data)
B2_counts <- rename(B2_counts, NAME = country)
B2_world_studies <- merge(world, B2_counts, by = "NAME", all = TRUE)

B2_world_studies <- B2_world_studies %>% filter(!NAME == "European Union")
B2_world_complete <- mutate(B2_world_threats, nstudies = B2_world_studies$n)
B2_world_complete <- B2_world_complete[!is.na(B2_world_complete$threats),]
B2_world_complete$nstudies[is.na(B2_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(B2_world_complete$threats, B2_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
B2_world_complete$threats_bin <- cut(B2_world_complete$threats, breaks = c(0, 5, 20, 60, 115), include.lowest = TRUE)
B2_world_complete$nstudies_bin <- cut(B2_world_complete$nstudies, breaks = c(0, 0.5, 20, 75, 139), include.lowest = TRUE)
classes <- bi_class(B2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(B2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot



#### BIRDS AND POLLUTION - B9 ####

## getting threats specific to birds and pollution 
B9_threats <- threats_data[!is.na(threats_data$B9),]

# applying quality filter to threats data 
B9_threats <- B9_threats %>% filter(BSR > 10)

# calculating risk
B9_threats$B9 <- B9_threats$B9 * B9_threats$BSR

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
B9_threat_country <- over(world_sp, B9_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
B9_threat_country <- B9_threat_country$B9
B9_world_threats <- mutate(world, threats = B9_threat_country)
B9_world_threats <- B9_world_threats[order(B9_world_threats$NAME),]#ordering alphabetically according to country name
B9_world_threats_only <- B9_world_threats[!is.na(B9_world_threats$threats),]
summary(B9_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = B9_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

##creating rank tables
B9_threats_rank <- subset(B9_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(B9_threats_rank) <- NULL
B9_threats_rank <- B9_threats_rank[order(-B9_threats_rank$threats),]



## creating a complete dataframe (this contains all countries for which there is threat data)
B9_counts <- rename(B9_counts, NAME = country)
B9_world_studies <- merge(world, B9_counts, by = "NAME", all = TRUE)

B9_world_complete <- mutate(B9_world_threats, nstudies = B9_world_studies$n)
B9_world_complete <- B9_world_complete[!is.na(B9_world_complete$threats),]
B9_world_complete$nstudies[is.na(B9_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(B9_world_complete$threats, B9_world_complete$nstudies, method = 'spearman')

## bivariate heatmap
B9_world_complete$threats_bin <- cut(B9_world_complete$threats, breaks = c(0, 0.50, 2, 5, 8), include.lowest = TRUE)
B9_world_complete$nstudies_bin <- cut(B9_world_complete$nstudies, breaks = c(0, 0.5, 5, 12, 24), include.lowest = TRUE)
classes <- bi_class(B9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(B9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot



#### BIRDS AND CLIMATE CHANGE - B11 ####

## getting threats specific to birds and climate change 
B11_threats <- threats_data[!is.na(threats_data$B11),]

# applying quality filter to threats data 
B11_threats <- B11_threats %>% filter(BSR > 10)

# calculating risk
B11_threats$B11 <- B11_threats$B11 * B11_threats$BSR

## getting studies specific to birds and climate change 
B_studies <- studies_data[studies_data$taxa == "Birds",]
B11_studies <- B_studies[B_studies$threattype == "Climate change & severe weather",]


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
B11_threat_country <- over(world_sp, B11_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
B11_threat_country <- B11_threat_country$B11
B11_world_threats <- mutate(world, threats = B11_threat_country)
B11_world_threats <- B11_world_threats[order(B11_world_threats$NAME),]#ordering alphabetically according to country name
B11_world_threats_only <- B11_world_threats[!is.na(B11_world_threats$threats),]
summary(B11_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = B11_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
B11_threats_rank <- subset(B11_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(B11_threats_rank) <- NULL
B11_threats_rank <- B11_threats_rank[order(-B11_threats_rank$threats),]



## creating a complete dataframe (this contains all countries for which there is threat data)
B11_counts <- rename(B11_counts, NAME = country)
B11_world_studies <- merge(world, B11_counts, by = "NAME", all = TRUE)

B11_world_complete <- mutate(B11_world_threats, nstudies = B11_world_studies$n)
B11_world_complete <- B11_world_complete[!is.na(B11_world_complete$threats),]
B11_world_complete$nstudies[is.na(B11_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(B11_world_complete$threats, B11_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
B11_world_complete$threats_bin <- cut(B11_world_complete$threats, breaks = c(0, 5, 20, 52), include.lowest = TRUE)
B11_world_complete$nstudies_bin <- cut(B11_world_complete$nstudies, breaks = c(0, 0.5, 1, 2), include.lowest = TRUE)
classes <- bi_class(B11_world_complete, x = threats_bin, y = nstudies_bin, dim = 3)
breaks1 <- bi_class_breaks(B11_world_complete, x = threats_bin, y = nstudies_bin, dim = 3, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink", dim = 3) + bi_theme()
legend <- bi_legend(pal = "GrPink", dim = 3, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot



#### BIRDS AND INVASIVES - B8_1 ####

## getting threats specific to birds and invasives
B8_1_threats <- threats_data[!is.na(threats_data$B8_1),]

# applying quality filter to threats data 
B8_1_threats <- B8_1_threats %>% filter(BSR > 10)

# calculating risk
B8_1_threats$B8_1 <- B8_1_threats$B8_1 * B8_1_threats$BSR

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
B8_1_threat_country <- over(world_sp, B8_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
B8_1_threat_country <- B8_1_threat_country$B8_1
B8_1_world_threats <- mutate(world, threats = B8_1_threat_country)
B8_1_world_threats <- B8_1_world_threats[order(B8_1_world_threats$NAME),]#ordering alphabetically according to country name
B8_1_world_threats_only <- B8_1_world_threats[!is.na(B8_1_world_threats$threats),]
summary(B8_1_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = B8_1_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
B8_1_threats_rank <- subset(B8_1_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(B8_1_threats_rank) <- NULL
B8_1_threats_rank <- B8_1_threats_rank[order(-B8_1_threats_rank$threats),]



## creating a complete dataframe (this contains all countries for which there is threat data)
B8_1_counts <- rename(B8_1_counts, NAME = country)
B8_1_world_studies <- merge(world, B8_1_counts, by = "NAME", all = TRUE)

B8_1_world_complete <- mutate(B8_1_world_threats, nstudies = B8_1_world_studies$n)
B8_1_world_complete <- B8_1_world_complete[!is.na(B8_1_world_complete$threats),]
B8_1_world_complete$nstudies[is.na(B8_1_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(B8_1_world_complete$threats, B8_1_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
B8_1_world_complete$threats_bin <- cut(B8_1_world_complete$threats, breaks = c(0, 5, 10, 20, 51), include.lowest = TRUE)
B8_1_world_complete$nstudies_bin <- cut(B8_1_world_complete$nstudies, breaks = c(0, 0.5, 1, 3, 5), include.lowest = TRUE)
classes <- bi_class(B8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(B8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot



#### BIRDS AND HUNTING - B5_1 ####

## getting threats specific to birds and hunting
B5_1_threats <- threats_data[!is.na(threats_data$B5_1),]

# applying quality filter to threats data 
B5_1_threats <- B5_1_threats %>% filter(BSR > 10)

# calculating risk
B5_1_threats$B5_1 <- B5_1_threats$B5_1 * B5_1_threats$BSR

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
B5_1_studies <- subset(B5_1_studies, select = -c(region, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
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
B5_1_threat_country <- over(world_sp, B5_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
B5_1_threat_country <- B5_1_threat_country$B5_1
B5_1_world_threats <- mutate(world, threats = B5_1_threat_country)
B5_1_world_threats <- B5_1_world_threats[order(B5_1_world_threats$NAME),]#ordering alphabetically according to country name
B5_1_world_threats_only <- B5_1_world_threats[!is.na(B5_1_world_threats$threats),]
summary(B5_1_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = B5_1_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')


## creating rank tables
B5_1_threats_rank <- subset(B5_1_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(B5_1_threats_rank) <- NULL
B5_1_threats_rank <- B5_1_threats_rank[order(-B5_1_threats_rank$threats),]


## creating a complete dataframe (this contains all countries for which there is threat data)
B5_1_counts <- rename(B5_1_counts, NAME = country)
B5_1_world_studies <- merge(world, B5_1_counts, by = "NAME", all = TRUE)

B5_1_world_complete <- mutate(B5_1_world_threats, nstudies = B5_1_world_studies$n)
B5_1_world_complete <- B5_1_world_complete[!is.na(B5_1_world_complete$threats),]
B5_1_world_complete$nstudies[is.na(B5_1_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(B5_1_world_complete$threats, B5_1_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
B5_1_world_complete$threats_bin <- cut(B5_1_world_complete$threats, breaks = c(0, 5, 15, 30, 46), include.lowest = TRUE)
B5_1_world_complete$nstudies_bin <- cut(B5_1_world_complete$nstudies, breaks = c(0, 0.5, 1, 2, 3), include.lowest = TRUE)
classes <- bi_class(B5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(B5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot



#### BIRDS AND LOGGING - B5_3 ####

## getting threats specific to birds and logging
B5_3_threats <- threats_data[!is.na(threats_data$B5_3),]

# applying quality filter to threats data 
B5_3_threats <- B5_3_threats %>% filter(BSR > 10)

# calculating risk
B5_3_threats$B5_3 <- B5_3_threats$B5_3 * B5_3_threats$BSR

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
B5_3_threat_country <- over(world_sp, B5_3_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
B5_3_threat_country <- B5_3_threat_country$B5_3
B5_3_world_threats <- mutate(world, threats = B5_3_threat_country)
B5_3_world_threats <- B5_3_world_threats[order(B5_3_world_threats$NAME),]#ordering alphabetically according to country name
B5_3_world_threats_only <- B5_3_world_threats[!is.na(B5_3_world_threats$threats),]
summary(B5_3_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = B5_3_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
B5_3_threats_rank <- subset(B5_3_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(B5_3_threats_rank) <- NULL
B5_3_threats_rank <- B5_3_threats_rank[order(-B5_3_threats_rank$threats),]



## creating a complete dataframe (this contains all countries for which there is threat data)
B5_3_counts <- rename(B5_3_counts, NAME = country)
B5_3_world_studies <- merge(world, B5_3_counts, by = "NAME", all = TRUE)

B5_3_world_complete <- mutate(B5_3_world_threats, nstudies = B5_3_world_studies$n)
B5_3_world_complete <- B5_3_world_complete[!is.na(B5_3_world_complete$threats),]
B5_3_world_complete$nstudies[is.na(B5_3_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(B5_3_world_complete$threats, B5_3_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
B5_3_world_complete$threats_bin <- cut(B5_3_world_complete$threats, breaks = c(0, 5, 15, 50, 112), include.lowest = TRUE)
B5_3_world_complete$nstudies_bin <- cut(B5_3_world_complete$nstudies, breaks = c(0, 0.5, 2, 6, 15), include.lowest = TRUE)
classes <- bi_class(B5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(B5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot



#### MAMMALS AND AGRICULTURE - M2 ####

## getting threats specific to mammals and agriculture 
M2_threats <- threats_data[!is.na(threats_data$M2),]

# applying quality filter to threats data 
M2_threats <- M2_threats %>% filter(MSR > 10)

# calculating risk
M2_threats$M2 <- M2_threats$M2 * M2_threats$MSR

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
M2_threat_country <- over(world_sp, M2_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
M2_threat_country <- M2_threat_country$M2
M2_world_threats <- mutate(world, threats = M2_threat_country)
M2_world_threats <- M2_world_threats[order(M2_world_threats$NAME),]#ordering alphabetically according to country name
M2_world_threats_only <- M2_world_threats[!is.na(M2_world_threats$threats),]
summary(M2_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = M2_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')


## creating rank tables
M2_threats_rank <- subset(M2_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(M2_threats_rank) <- NULL
M2_threats_rank <- M2_threats_rank[order(-M2_threats_rank$threats),]



## creating a complete dataframe (this contains all countries for which there is threat data)
M2_counts <- rename(M2_counts, NAME = country)
M2_world_studies <- merge(world, M2_counts, by = "NAME", all = TRUE)

M2_world_complete <- mutate(M2_world_threats, nstudies = M2_world_studies$n)
M2_world_complete <- M2_world_complete[!is.na(M2_world_complete$threats),]
M2_world_complete$nstudies[is.na(M2_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(M2_world_complete$threats, M2_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
M2_world_complete$threats_bin <- cut(M2_world_complete$threats, breaks = c(0, 5, 20, 50, 84), include.lowest = TRUE)
M2_world_complete$nstudies_bin <- cut(M2_world_complete$nstudies, breaks = c(0, 0.5,5, 50, 113), include.lowest = TRUE)
classes <- bi_class(M2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(M2_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot


#### MAMMALS AND POLLUTION - M9 ####

## getting threats specific to mammals and pollution 
M9_threats <- threats_data[!is.na(threats_data$M9),]

# applying quality filter to threats data 
M9_threats <- M9_threats %>% filter(MSR > 10)

# calculating risk
M9_threats$M9 <- M9_threats$M9 * M9_threats$MSR

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
M9_threat_country <- over(world_sp, M9_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
M9_threat_country <- M9_threat_country$M9
M9_world_threats <- mutate(world, threats = M9_threat_country)
M9_world_threats <- M9_world_threats[order(M9_world_threats$NAME),]#ordering alphabetically according to country name
M9_world_threats_only <- M9_world_threats[!is.na(M9_world_threats$threats),]
summary(M9_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = M9_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')


## creating rank tables
M9_threats_rank <- subset(M9_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(M9_threats_rank) <- NULL
M9_threats_rank <- M9_threats_rank[order(-M9_threats_rank$threats),]



## creating a complete dataframe (this contains all countries for which there is threat data)
M9_counts <- rename(M9_counts, NAME = country)
M9_world_studies <- merge(world, M9_counts, by = "NAME", all = TRUE)

M9_world_complete <- mutate(M9_world_threats, nstudies = M9_world_studies$n)
M9_world_complete <- M9_world_complete[!is.na(M9_world_complete$threats),]
M9_world_complete$nstudies[is.na(M9_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(M9_world_complete$threats, M9_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
M9_world_complete$threats_bin <- cut(M9_world_complete$threats, breaks = c(0, 0.5, 1, 2, 5), include.lowest = TRUE)
M9_world_complete$nstudies_bin <- cut(M9_world_complete$nstudies, breaks = c(0, 0.5, 3, 10, 16), include.lowest = TRUE)
classes <- bi_class(M9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(M9_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot


#### MAMMALS AND CLIMATE CHANGE - M11 ####

## getting threats specific to mammals and climate change 
M11_threats <- threats_data[!is.na(threats_data$M11),]

# applying quality filter to threats data 
M11_threats <- M11_threats %>% filter(MSR > 10)

# calculating risk 
M11_threats$M11 <- M11_threats$M11 * M11_threats$MSR

## getting studies specific to mammals and climate change 
M_studies <- studies_data[studies_data$taxa == "Mammals",]
M11_studies <- M_studies[M_studies$threattype == "Climate change & severe weather",]


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
M11_threat_country <- over(world_sp, M11_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
M11_threat_country <- M11_threat_country$M11
M11_world_threats <- mutate(world, threats = M11_threat_country)
M11_world_threats <- M11_world_threats[order(M11_world_threats$NAME),]#ordering alphabetically according to country name
M11_world_threats_only <- M11_world_threats[!is.na(M11_world_threats$threats),]
summary(M11_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = M11_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
M11_threats_rank <- subset(M11_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(M11_threats_rank) <- NULL
M11_threats_rank <- M11_threats_rank[order(-M11_threats_rank$threats),]



## creating a complete dataframe (this contains all countries for which there is threat data)
M11_counts <- rename(M11_counts, NAME = country)
M11_world_studies <- merge(world, M11_counts, by = "NAME", all = TRUE)

M11_world_complete <- mutate(M11_world_threats, nstudies = M11_world_studies$n)
M11_world_complete <- M11_world_complete[!is.na(M11_world_complete$threats),]
M11_world_complete$nstudies[is.na(M11_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(M11_world_complete$threats, M11_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
M11_world_complete$threats_bin <- cut(M11_world_complete$threats, breaks = c(0, 0.2, 1, 3, 6), include.lowest = TRUE)
M11_world_complete$nstudies_bin <- cut(M11_world_complete$nstudies, breaks = c(0, 0.5, 1, 3, 5), include.lowest = TRUE)
classes <- bi_class(M11_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(M11_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot



#### MAMMALS AND INVASIVES - M8_1 ####

## getting threats specific to mammals and invasives
M8_1_threats <- threats_data[!is.na(threats_data$M8_1),]

# applying quality filter to threats data 
M8_1_threats <- M8_1_threats %>% filter(MSR > 10)

# calculating risk
M8_1_threats$M8_1 <- M8_1_threats$M8_1 * M8_1_threats$MSR

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
M8_1_studies <- subset(M8_1_studies, select = -c(region, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
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
M8_1_threat_country <- over(world_sp, M8_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
M8_1_threat_country <- M8_1_threat_country$M8_1
M8_1_world_threats <- mutate(world, threats = M8_1_threat_country)
M8_1_world_threats <- M8_1_world_threats[order(M8_1_world_threats$NAME),]#ordering alphabetically according to country name
M8_1_world_threats_only <- M8_1_world_threats[!is.na(M8_1_world_threats$threats),]
summary(M8_1_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = M8_1_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
M8_1_threats_rank <- subset(M8_1_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(M8_1_threats_rank) <- NULL
M8_1_threats_rank <- M8_1_threats_rank[order(-M8_1_threats_rank$threats),]



## creating a complete dataframe (this contains all countries for which there is threat data)
M8_1_counts <- rename(M8_1_counts, NAME = country)
M8_1_world_studies <- merge(world, M8_1_counts, by = "NAME", all = TRUE)

M8_1_world_complete <- mutate(M8_1_world_threats, nstudies = M8_1_world_studies$n)
M8_1_world_complete <- M8_1_world_complete[!is.na(M8_1_world_complete$threats),]
M8_1_world_complete$nstudies[is.na(M8_1_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(M8_1_world_complete$threats, M8_1_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
M8_1_world_complete$threats_bin <- cut(M8_1_world_complete$threats, breaks = c(0, 2, 5, 13), include.lowest = TRUE)
M8_1_world_complete$nstudies_bin <- cut(M8_1_world_complete$nstudies, breaks = c(0, 0.5, 1, 2), include.lowest = TRUE)
classes <- bi_class(M8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 3)
breaks1 <- bi_class_breaks(M8_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 3, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink", dim = 3) + bi_theme()
legend <- bi_legend(pal = "GrPink", dim = 3, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot



#### MAMMALS AND HUNTING - M5_1 ####

## getting threats specific to mammals and hunting
M5_1_threats <- threats_data[!is.na(threats_data$M5_1),]

# applying quality filter to threats data 
M5_1_threats <- M5_1_threats %>% filter(MSR > 10)

# calculating risk
M5_1_threats$M5_1 <- M5_1_threats$M5_1 * M5_1_threats$MSR

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
M5_1_threat_country <- over(world_sp, M5_1_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
M5_1_threat_country <- M5_1_threat_country$M5_1
M5_1_world_threats <- mutate(world, threats = M5_1_threat_country)
M5_1_world_threats <- M5_1_world_threats[order(M5_1_world_threats$NAME),]#ordering alphabetically according to country name
M5_1_world_threats_only <- M5_1_world_threats[!is.na(M5_1_world_threats$threats),]
summary(M5_1_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = M5_1_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')

## creating rank tables
M5_1_threats_rank <- subset(M5_1_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(M5_1_threats_rank) <- NULL
M5_1_threats_rank <- M5_1_threats_rank[order(-M5_1_threats_rank$threats),]



## creating a complete dataframe (this contains all countries for which there is threat data)
M5_1_counts <- rename(M5_1_counts, NAME = country)
M5_1_world_studies <- merge(world, M5_1_counts, by = "NAME", all = TRUE)

M5_1_world_complete <- mutate(M5_1_world_threats, nstudies = M5_1_world_studies$n)
M5_1_world_complete <- M5_1_world_complete[!is.na(M5_1_world_complete$threats),]
M5_1_world_complete$nstudies[is.na(M5_1_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(M5_1_world_complete$threats, M5_1_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
M5_1_world_complete$threats_bin <- cut(M5_1_world_complete$threats, breaks = c(0, 5, 12, 25, 46), include.lowest = TRUE)
M5_1_world_complete$nstudies_bin <- cut(M5_1_world_complete$nstudies, breaks = c(0, 0.5, 5, 20, 41), include.lowest = TRUE)
classes <- bi_class(M5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(M5_1_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot



#### MAMMALS AND LOGGING - M5_3 ####

## getting threats specific to mammals and logging
M5_3_threats <- threats_data[!is.na(threats_data$M5_3),]

# applying quality filter to threats data 
M5_3_threats <- M5_3_threats %>% filter(MSR > 10)

# calculating risk
M5_3_threats$M5_3 <- M5_3_threats$M5_3 * M5_3_threats$MSR

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
M5_3_threat_country <- over(world_sp, M5_3_centroids_sp, fn = mean)
#this last command seems to return a dataframe with the 246 countries, and average values of each attribute from the threats file per country!


## plotting threats
M5_3_threat_country <- M5_3_threat_country$M5_3
M5_3_world_threats <- mutate(world, threats = M5_3_threat_country)
M5_3_world_threats <- M5_3_world_threats[order(M5_3_world_threats$NAME),]#ordering alphabetically according to country name
M5_3_world_threats_only <- M5_3_world_threats[!is.na(M5_3_world_threats$threats),]
summary(M5_3_world_threats_only$threats)

ggplot() + geom_sf(data = world, fill = 'dark grey', colour = 'black') + geom_sf(data = M5_3_world_threats_only, colour = 'black', aes(fill = threats)) + 
  scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red') + theme_void() + theme(legend.position = c(0.16, 0.3), 
                                                                                         legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=17)) + labs(fill = 'Conservation risk ')


## creating rank tables
M5_3_threats_rank <- subset(M5_3_world_threats_only, select = -c(FIPS, ISO2, ISO3, UN, AREA, POP2005, REGION, SUBREGION, LON, LAT))
st_geometry(M5_3_threats_rank) <- NULL
M5_3_threats_rank <- M5_3_threats_rank[order(-M5_3_threats_rank$threats),]



## creating a complete dataframe (this contains all countries for which there is threat data)
M5_3_counts <- rename(M5_3_counts, NAME = country)
M5_3_world_studies <- merge(world, M5_3_counts, by = "NAME", all = TRUE)

M5_3_world_complete <- mutate(M5_3_world_threats, nstudies = M5_3_world_studies$n)
M5_3_world_complete <- M5_3_world_complete[!is.na(M5_3_world_complete$threats),]
M5_3_world_complete$nstudies[is.na(M5_3_world_complete$nstudies)] <- 0


## correlation between threats and studies 
cor.test(M5_3_world_complete$threats, M5_3_world_complete$nstudies, method = 'spearman')


## bivariate heatmap
M5_3_world_complete$threats_bin <- cut(M5_3_world_complete$threats, breaks = c(0, 3, 12, 35, 88), include.lowest = TRUE)
M5_3_world_complete$nstudies_bin <- cut(M5_3_world_complete$nstudies, breaks = c(0, 0.5, 5, 15, 30), include.lowest = TRUE)
classes <- bi_class(M5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4)
breaks1 <- bi_class_breaks(M5_3_world_complete, x = threats_bin, y = nstudies_bin, dim = 4, dig_lab = c(x = 5, y = 2), split = TRUE)
map <- ggplot() + geom_sf(data = world, fill = 'white', colour = 'black') + geom_sf(data = classes, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) + bi_scale_fill(pal = "GrPink2", dim = 4) + bi_theme()
legend <- bi_legend(pal = "GrPink2", dim = 4, xlab = "Conservation risk", ylab = "No. of studies", size = 23, breaks = breaks1, arrows = FALSE)
finalPlot <- ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot

