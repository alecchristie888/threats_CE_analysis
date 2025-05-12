## loading necessary libraries
library(sf)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(readxl)
library(svglite)
library(ggpubr)

## reading in files
studies_data <- read.csv("Data/CE_data_complete.csv")
non_english <- read_excel("Data/CE_non_english_data_clean.xlsx")

length(unique(studies_data$pageid))
length(unique(non_english$PaperID))

## cleaning up data
data_clean <- studies_data[!is.na(studies_data$long),] #remove reviews that have no coordinates
data <- filter(data_clean, !(lat == 0 & long == 0)) #remove reviews that have no coordinates
#non_english[is.na(non_english$long),] #no reviews
#filter(non_english, (lat == 0 & long == 0)) #no reviews

## creating a dataframe of studies used
agri <- data[data$threattype == "Agriculture & aquaculture",]
agri <- agri[agri$threatname == "Annual/perennial non-timber crops" | agri$threatname == "Livestock farm/ranch" | 
                                       agri$threatname == "Wood/pulp plantations" | is.na(agri$threatname),]
agri <- cbind(agri,data.frame(threat="Agriculture"))

poll <- cbind(data[data$threattype == "Pollution",],data.frame(threat="Pollution"))

cc <- cbind(data[data$threattype == "Climate change & severe weather",],data.frame(threat="Climate Change"))
              
log <- cbind(data[data$threatname == "Logging/wood harvesting",],data.frame(threat="Logging"))

inv <- cbind(data[data$threatname == "Invasive non-native species",],data.frame(threat="Invasives"))

hunt <- cbind(data[data$threatname == "Hunting/trapping terrestrial animals",],data.frame(threat="Hunting"))

studies <- rbind(agri, poll, cc, log, inv, hunt)
    

## creating a dataframe of studies used for non-english studies
ne.agri <- non_english[non_english$threattype == "Agriculture & aquaculture",]
ne.agri <- ne.agri[ne.agri$threatname == "Annual/perennial non-timber crops" | ne.agri$threatname == "Livestock farm/ranch" | 
                  ne.agri$threatname == "Wood/pulp plantations" | is.na(ne.agri$threatname),]
ne.agri <- cbind(ne.agri,data.frame(threat="Agriculture"))

ne.poll <- cbind(non_english[non_english$threattype == "Pollution",],data.frame(threat="Pollution"))

#no studies
#ne.cc <- cbind(non_english[non_english$threattype == "Climate change & severe weather",],data.frame(threat="Climate Change"))

ne.log <- cbind(non_english[non_english$threatname == "Logging/wood harvesting",],data.frame(threat="Logging"))
ne.log <- ne.log[!is.na(ne.log$PaperID),]

ne.inv <- cbind(non_english[non_english$threatname == "Invasive non-native species",],data.frame(threat="Invasives"))
ne.inv <- ne.inv[!is.na(ne.inv$PaperID),]

ne.hunt <- cbind(non_english[non_english$threatname == "Hunting/trapping terrestrial animals",],data.frame(threat="Hunting"))
ne.hunt <- ne.hunt[!is.na(ne.hunt$PaperID),]

ne.studies.all <- rbind(ne.agri, ne.poll, ne.log, ne.inv, ne.hunt)
            
## finding total number of studies used
length(unique(studies$pageid)) #1025 studies
length(unique(ne.studies.all$PaperID)) #147 studies

unique(ne.studies.all$Language)

## finding rank tables of countries
# for english studies
studies_country <- subset(studies, select = -c(lat, long, syn, region, threattype, threatname, action_type, action_name, action_sub, eff_cat, eff_score, cert_score, harm_score, int, habitat, rowid))
studies_country <- studies_country[!is.na(studies_country$country),]
country_counts <- distinct(studies_country) %>% count(country,threat,taxa)
country_counts <- cbind(country_counts, data.frame(language="English"))


# for non-english studies 
non_english_country_counts <- subset(ne.studies.all, select = -c(lat, long, threattype, threatname)) %>% 
  distinct() %>%  
  count(country,threat,taxa)
non_english_country_counts <- cbind(non_english_country_counts, data.frame(language="Non-English"))
  
# merging the two
country_counts.all <- merge(country_counts, non_english_country_counts, by = c("country","threat","taxa","language"), all = TRUE)
country_counts.all$n.x[is.na(country_counts.all$n.x)] <- 0
country_counts.all$n.y[is.na(country_counts.all$n.y)] <- 0
country_counts.all$n <- country_counts.all$n.x + country_counts.all$n.y
country_counts.all.plot <- subset(country_counts.all, select = -c(n.x, n.y))
country_counts.all.plot.eng <- subset(country_counts.all.plot, language == "English")
country_counts.all.plot.non.eng <- subset(country_counts.all.plot, language == "Non-English")
country_counts.all.plot.order.eng <- country_counts.all.plot.eng[rev(order(country_counts.all.plot.eng$n)),]
country_counts.all.plot.order.non.eng <- country_counts.all.plot.non.eng[rev(order(country_counts.all.plot.non.eng$n)),]

# plotting!
ampbird <-ggplot(country_counts.all.plot.order.eng[country_counts.all.plot.order.eng$language=="English"&(country_counts.all.plot.order.eng$taxa=="Amphibians"|country_counts.all.plot.order.eng$taxa=="Birds"),]) + 
    geom_col(aes(x=n, y=reorder_within(country,n,taxa,fun=sum), fill=threat)) +
    scale_x_continuous(expand = c(0,0),breaks=seq(0,160,20)) + 
    theme_classic() + 
    scale_fill_manual(name="Threat",values=viridis::viridis_pal(option = "B")(10)[c(9,2,6,5,8,3)])+
    theme(axis.ticks = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=30),
        legend.position = "none",
        strip.text.x = element_text(face = "bold", color = "white", hjust = 1, size = 25),
        strip.text.y = element_text(face = "bold", color = "white", hjust = 0, size = 25),
        strip.background.x = element_rect(fill = "seagreen", linetype = "solid",
                                          color = "black", linewidth = 1),
        strip.background.y = element_rect(fill = "black", linetype = "solid",
                                          color = "gray30", linewidth = 1)) + 
    labs(x = "Number of studies", y = "Country")+
    facet_wrap(~taxa,scales="free",ncol=1)+
    scale_y_reordered(labels=function(x) gsub("___.+$", "", x))

mam <-ggplot(country_counts.all.plot.order.eng[country_counts.all.plot.order.eng$language=="English"&country_counts.all.plot.order.eng$taxa=="Mammals",]) + 
  geom_col(aes(x=n, y=reorder_within(country,n,taxa,fun=sum), fill=threat)) +
  scale_x_continuous(expand = c(0,0),breaks=seq(0,200,20), limits=c(0,210)) + 
  theme_classic() + 
  scale_fill_manual(name="Threat",values=viridis::viridis_pal(option = "B")(10)[c(9,2,6,5,8,3)])+
  theme(axis.ticks = element_blank(),
        axis.text=element_text(size=20),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size = 30),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 30),
        legend.position = "inside",
        strip.text.x = element_text(face = "bold", color = "white", hjust = 1, size = 25),
        strip.text.y = element_text(face = "bold", color = "white", hjust = 0, size = 25),
        strip.background.x = element_rect(fill = "seagreen", linetype = "solid",
                                          color = "black", linewidth = 1),
        strip.background.y = element_rect(fill = "black", linetype = "solid",
                                          color = "gray30", linewidth = 1)) + 
  labs(x = "Number of studies", y = "Country")+
  facet_wrap(language~taxa,scales="free_y",ncol=1)+
  scale_y_reordered(labels=function(x) gsub("___.+$", "", x))

ggarrange(ampbird,mam,ncol=2)

ggsave("barplot_english.svg", height = 30, width = 20,dpi=600,device="svg")


library(tidytext)
ggplot(country_counts.all.plot.order.non.eng[country_counts.all.plot.order.non.eng$language=="Non-English",]) + 
  geom_col(aes(x=n, y=reorder_within(country,n,taxa,fun=sum), fill=threat)) +
  scale_x_continuous(expand = c(0,0),breaks=seq(0,20,2),limits=c(0,21)) + 
  theme_classic() + 
  scale_fill_manual(name="Threat",values=viridis::viridis_pal(option = "B")(10)[c(9,2,6,5,8,3)])+
  theme(axis.ticks = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=30),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 30),
        legend.position = "inside",
        strip.placement = "outside",
        strip.text.x = element_text(face = "bold", color = "white", hjust = 1, size = 30),
        strip.text.y = element_text(face = "bold", color = "white", hjust = 0, size = 30),
        strip.background.x = element_rect(fill = "cornflowerblue", linetype = "solid",
                                          color = "black", linewidth = 1),
        strip.background.y = element_rect(fill = "black", linetype = "solid",
                                          color = "gray30", linewidth = 1)) +  
  labs(x = "Number of studies", y = "Country")+
  facet_wrap(~language+taxa,scales="free_y",ncol=1)+
  scale_y_reordered(labels=function(x) gsub("___.+$", "", x))

ggsave("barplot_non-english.svg", height = 30, width = 20,dpi=600,device="svg")





length(unique(country_counts.all.plot.order$country))
country_counts.all.plot.order %>% group_by(country) %>% summarise(sum.studies=sum(n)) %>% arrange(desc(sum.studies))
country_counts.all.plot.order %>% group_by(country) %>% summarise(sum.studies=sum(n)/1165*100) %>% arrange(desc(sum.studies))

country_counts.all.plot.order %>% filter(language=="English") %>% group_by(country) %>% summarise(sum.studies=sum(n)/1025*100) %>% arrange(desc(sum.studies))
country_counts.all.plot.order %>% filter(language=="Non-English") %>% group_by(country) %>% summarise(sum.studies=sum(n)/147*100) %>% arrange(desc(sum.studies))


country_counts.all.plot.order %>% group_by(country) %>% summarise(sum.studies=sum(n)) %>% arrange(desc(sum.studies)) %>% filter(sum.studies<6) %>% print(n=100)
51/92

####################

nrow(unique(studies[,c("pageid","country")])) # 1074 unique study-location combinations 
length(unique(studies[,c("pageid")])) #1025 unique studies, so 49 studies in more than one country
nrow(unique(ne.studies.all[,c("PaperID","country")])) # 147 unique study-location combinations 
length(unique(ne.studies.all[,c("PaperID")])) #147 unique studies, so no studies in more than one country


length(unique(studies$pageid)) #1025 studies English language
length(unique(ne.studies.all$PaperID)) #147 studies non-English language


## Studies per threat type
eng.agri <- length(unique(studies$pageid[studies$threattype == "Agriculture & aquaculture"])) 
non.eng.agri <- length(unique(ne.studies.all$PaperID[ne.studies.all$threattype == "Agriculture & aquaculture"])) 


eng.poll <- length(unique(studies$pageid[studies$threattype == "Pollution"])) 
non.eng.poll <- length(unique(ne.studies.all$PaperID[ne.studies.all$threattype == "Pollution"])) 


eng.cli <- length(unique(studies$pageid[studies$threattype == "Climate change & severe weather"])) 
non.eng.cli <- length(unique(ne.studies.all$PaperID[ne.studies.all$threattype == "Climate change & severe weather"])) 


eng.inv <- length(unique(studies$pageid[studies$threatname == "Invasive non-native species"])) 
non.eng.inv <- length(unique(ne.studies.all$PaperID[ne.studies.all$threatname == "Invasive non-native species"])) 


eng.hun <- length(unique(studies$pageid[studies$threatname == "Hunting/trapping terrestrial animals"])) 
non.eng.hun <- length(unique(ne.studies.all$PaperID[ne.studies.all$threatname == "Hunting/trapping terrestrial animals"])) 


eng.log <- length(unique(studies$pageid[studies$threatname == "Logging/wood harvesting"])) 
non.eng.log <- length(unique(ne.studies.all$PaperID[ne.studies.all$threatname == "Logging/wood harvesting"])) 


sum.dat <- rbind(data.frame(n=c(eng.log,eng.poll,eng.cli,eng.inv,eng.hun,eng.agri),
           language="English",
           threat=c("Logging","Pollution","Climate Change","Invasives","Hunting","Agriculture")),
data.frame(n=c(non.eng.log,non.eng.poll,non.eng.cli,non.eng.inv,non.eng.hun,non.eng.agri),
           language="Non-English",
           threat=c("Logging","Pollution","Climate Change","Invasives","Hunting","Agriculture")))

sum(sum.dat$n)
sum(sum.dat$n)-1025-147
sum.dat %>% group_by(threat) %>% summarise(sum(n))
sum.dat %>% group_by(language) %>% summarise(sum(n))

sum.dat %>% filter(language=="English")%>% group_by(threat) %>% summarise(sum(n)/1025)
sum.dat %>% filter(language=="Non-English")%>% group_by(threat) %>% summarise(sum(n)/147)
sum.dat %>% group_by(threat) %>% summarise(sum(n)/1165)

## Studies per taxa 
amp.agri <- length(unique(studies$pageid[studies$taxa=="Amphibians"&studies$threattype == "Agriculture & aquaculture"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Amphibians"&ne.studies.all$threattype == "Agriculture & aquaculture"])) 
bird.agri <- length(unique(studies$pageid[studies$taxa=="Birds"&studies$threattype == "Agriculture & aquaculture"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Birds"&ne.studies.all$threattype == "Agriculture & aquaculture"])) 
mam.agri <- length(unique(studies$pageid[studies$taxa=="Mammals"&studies$threattype == "Agriculture & aquaculture"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Mammals"&ne.studies.all$threattype == "Agriculture & aquaculture"])) 

amp.poll <- length(unique(studies$pageid[studies$taxa=="Amphibians"&studies$threattype == "Pollution"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Amphibians"&ne.studies.all$threattype == "Pollution"])) 
bird.poll <- length(unique(studies$pageid[studies$taxa=="Birds"&studies$threattype == "Pollution"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Birds"&ne.studies.all$threattype == "Pollution"])) 
mam.poll <- length(unique(studies$pageid[studies$taxa=="Mammals"&studies$threattype == "Pollution"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Mammals"&ne.studies.all$threattype == "Pollution"])) 

amp.cli <- length(unique(studies$pageid[studies$taxa=="Amphibians"&studies$threattype == "Climate change & severe weather"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Amphibians"&ne.studies.all$threattype == "Climate change & severe weather"])) 
bird.cli <- length(unique(studies$pageid[studies$taxa=="Birds"&studies$threattype == "Climate change & severe weather"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Birds"&ne.studies.all$threattype == "Climate change & severe weather"])) 
mam.cli <- length(unique(studies$pageid[studies$taxa=="Mammals"&studies$threattype == "Climate change & severe weather"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Mammals"&ne.studies.all$threattype == "Climate change & severe weather"])) 

amp.inv <- length(unique(studies$pageid[studies$taxa=="Amphibians"&studies$threatname == "Invasive non-native species"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Amphibians"&ne.studies.all$threatname == "Invasive non-native species"])) 
bird.inv <- length(unique(studies$pageid[studies$taxa=="Birds"&studies$threatname == "Invasive non-native species"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Birds"&ne.studies.all$threatname == "Invasive non-native species"])) 
mam.inv <- length(unique(studies$pageid[studies$taxa=="Mammals"&studies$threatname == "Invasive non-native species"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Mammals"&ne.studies.all$threatname == "Invasive non-native species"])) 

amp.hun <- length(unique(studies$pageid[studies$taxa=="Amphibians"&studies$threatname == "Hunting/trapping terrestrial animals"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Amphibians"&ne.studies.all$threatname == "Hunting/trapping terrestrial animals"])) 
bird.hun <- length(unique(studies$pageid[studies$taxa=="Birds"&studies$threatname == "Hunting/trapping terrestrial animals"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Birds"&ne.studies.all$threatname == "Hunting/trapping terrestrial animals"])) 
mam.hun <- length(unique(studies$pageid[studies$taxa=="Mammals"&studies$threatname == "Hunting/trapping terrestrial animals"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Mammals"&ne.studies.all$threatname == "Hunting/trapping terrestrial animals"])) 

amp.log <- length(unique(studies$pageid[studies$taxa=="Amphibians"&studies$threatname == "Logging/wood harvesting"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Amphibians"&ne.studies.all$threatname == "Logging/wood harvesting"])) 
bird.log <- length(unique(studies$pageid[studies$taxa=="Birds"&studies$threatname == "Logging/wood harvesting"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Birds"&ne.studies.all$threatname == "Logging/wood harvesting"])) 
mam.log <- length(unique(studies$pageid[studies$taxa=="Mammals"&studies$threatname == "Logging/wood harvesting"])) + length(unique(ne.studies.all$PaperID[ne.studies.all$taxa=="Mammals"&ne.studies.all$threatname == "Logging/wood harvesting"])) 


sum.dat.taxa <- rbind(data.frame(n=c(amp.log,amp.poll,amp.cli,amp.inv,amp.hun,amp.agri),
           threat=c("Logging","Pollution","Climate Change","Invasives","Hunting","Agriculture"),
           taxa="Amphibians"),
           data.frame(n=c(bird.log,bird.poll,bird.cli,bird.inv,bird.hun,bird.agri),
           threat=c("Logging","Pollution","Climate Change","Invasives","Hunting","Agriculture"),
           taxa="Birds"),
           data.frame(n=c(mam.log,mam.poll,mam.cli,mam.inv,mam.hun,mam.agri),
           threat=c("Logging","Pollution","Climate Change","Invasives","Hunting","Agriculture"),
           taxa="Mammals"))

sum(sum.dat.taxa$n)
sum.dat.taxa %>% group_by(threat) %>% summarise(sum(n))
sum.dat.taxa %>% group_by(taxa) %>% summarise(sum(n))

sum.dat.taxa %>% filter(taxa=="Amphibians")%>% group_by(threat) %>% summarise(sum(n)/173)
sum.dat.taxa %>% filter(taxa=="Birds")%>% group_by(threat) %>% summarise(sum(n)/392)
sum.dat.taxa %>% filter(taxa=="Mammals")%>% group_by(threat) %>% summarise(sum(n)/630)

length(unique(studies$pageid[studies$taxa == "Amphibians"])) 
length(unique(ne.studies.all$PaperID[ne.studies.all$taxa == "Amphibians"])) 
147+26
# 173

length(unique(studies$pageid[studies$taxa == "Birds"])) 
length(unique(ne.studies.all$PaperID[ne.studies.all$taxa == "Birds"])) 
297+95
# 392

length(unique(studies$pageid[studies$taxa == "Mammals"])) 
length(unique(ne.studies.all$PaperID[ne.studies.all$taxa == "Mammals"])) 
589+41
# 630

173+392+630
# 1195 #does not add up to 1165 as some studies study multiple taxa

