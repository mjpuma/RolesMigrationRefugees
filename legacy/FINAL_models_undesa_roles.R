### load package libraries
library(tidyverse)
# devtools::install_github("leifeld/tnam")
library(tnam)
library(readxl)
library(sna)
library(network)
library(igraph)
library(texreg)
library(matrixStats)
library(mclogit)

# set custom colors
custom_colors <- c(
  # "#9e0142",
                   "#d53e4f",
                   # "#f46d43",
                   # "#fdae61",
                   "#fee08b",
                   # "#e6f598",
                   # "#abdda4",
                   "#66c2a5"
                   # ,
                   # "#3288bd",
                   # "#5e4fa2"
                   )

# function to normalize from 0 to 1 for plotting arguments
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

### Set working directory
setwd("D:/Documents/Work/University of Florida/Postdoc 2021-2022/Project Documents/Roles of Mobility Networks")

##### Data Management #####
### Country Roles
undesa_roles <- read.csv("results/UNDESA_full_roles_k3.csv", header = T)
undesa_roles$X <- sub("\\.19.*", "", undesa_roles$X)
undesa_roles$X <- sub("\\.20.*", "", undesa_roles$X)
undesa_roles$id <- sub(".*\\.", "", undesa_roles$id)
colnames(undesa_roles) <- c("name","time","role")

### country reference from bilat_mig
States <- read.csv("data/country_list.csv", header = T)

# Migrant Flow Estimates 1990-2020
edges.mig <- read.csv("data/bilat_mig.csv", header = T)[,c(1:3,9)]
edges.mig$state.origin.name <- States$name[match(edges.mig$orig, States$iso3c)]
edges.mig$state.destination.name <- States$name[match(edges.mig$dest, States$iso3c)]
edges.mig <- edges.mig[,c(1,5,6,4)]
edges.mig$year0 <- edges.mig$year0+5
colnames(edges.mig)[1] <- "year5"
edges.mig.wide <- reshape(data = edges.mig,
                          varying = NULL,
                          v.names = NULL,
                          timevar = "year5",
                          idvar = c("state.origin.name", "state.destination.name"),
                          direction = "wide",
                          new.row.names = NULL,
                          sep = "_")
colnames(edges.mig.wide)[3:length(colnames(edges.mig.wide))] <- c("migflows.1995",
                                                                  "migflows.2000",
                                                                  "migflows.2005",
                                                                  "migflows.2010",
                                                                  "migflows.2015",
                                                                  "migflows.2020")

# UN POP data
pop <- read.csv("./data/WPP2022_Demographic_Indicators_Medium.csv") # call in UN world pop and demographic data
pop <- pop[as.numeric(pop$Time)>=1990 & as.numeric(pop$Time)<=2020,] # retain data between 1990 and 2020
pop <- pop[pop$Time%in%c("1995","2000","2005","2010","2015","2020"),]
pop <- pop[pop$LocTypeName == "Country/Area" & !pop$Location %in% c("Guernsey","Jersey"),] # selecting countries not equal to either of the Channel Islands (matching naming convention of country_list)
pop <- pop[!pop$Location %in% c("China, Taiwan Province of China", # removing additional countries and selecting appropriate columns
                                "Holy See",
                                "Kosovo (under UNSC res. 1244)",
                                "Saint Barthélemy",
                                "Saint Martin (French part)"),
           colnames(pop) %in% c("Location","Time","TPopulation1July")]
# changing names to match country_list
pop[which(pop$Location == "Türkiye"),]$Location <- "Turkey"
pop[which(pop$Location == "Antigua and Barbuda"),]$Location <- "Antigua & Barbuda"
pop[which(pop$Location == "Bolivia (Plurinational State of)"),]$Location <- "Bolivia"
pop[which(pop$Location == "Bosnia and Herzegovina"),]$Location <- "Bosnia & Herzegovina"
pop[which(pop$Location == "Brunei Darussalam"),]$Location <- "Brunei"
pop[which(pop$Location == "Cabo Verde"),]$Location <- "Cape Verde"
pop[which(pop$Location == "China, Hong Kong SAR"),]$Location <- "Hong Kong SAR"
pop[which(pop$Location == "China, Macao SAR"),]$Location <- "Macao SAR"
pop[which(pop$Location == "Côte d'Ivoire"),]$Location <- "Côte d’Ivoire"
pop[which(pop$Location == "Dem. People's Republic of Korea"),]$Location <- "North Korea"
pop[which(pop$Location == "Democratic Republic of the Congo"),]$Location <- "DR Congo"
pop[which(pop$Location == "Falkland Islands (Malvinas)"),]$Location <- "Falkland Islands"
pop[which(pop$Location == "Iran (Islamic Republic of)"),]$Location <- "Iran"
# pop[which(pop$Location == "Kosovo (under UNSC res. 1244)"),]$Location <- "Kosovo"
pop[which(pop$Location == "Lao People's Democratic Republic"),]$Location <- "Laos"
pop[which(pop$Location == "Micronesia (Fed. States of)"),]$Location <- "FS Micronesia"
pop[which(pop$Location == "Republic of Korea"),]$Location <- "South Korea"
pop[which(pop$Location == "Republic of Moldova"),]$Location <- "Moldova"
pop[which(pop$Location == "Russian Federation"),]$Location <- "Russia"
pop[which(pop$Location == "Saint Helena"),]$Location <- "St. Helena"
pop[which(pop$Location == "Saint Kitts and Nevis"),]$Location <- "St. Kitts & Nevis"
pop[which(pop$Location == "Saint Lucia"),]$Location <- "St. Lucia"
pop[which(pop$Location == "Saint Pierre and Miquelon"),]$Location <- "St. Pierre & Miquelon"
pop[which(pop$Location == "Saint Vincent and the Grenadines"),]$Location <- "St. Vincent & Grenadines"
pop[which(pop$Location == "Sao Tome and Principe"),]$Location <- "São Tomé & Príncipe"
pop[which(pop$Location == "Sint Maarten (Dutch part)"),]$Location <- "Sint Maarten"
pop[which(pop$Location == "State of Palestine"),]$Location <- "Palestine"
pop[which(pop$Location == "Syrian Arab Republic"),]$Location <- "Syria"
pop[which(pop$Location == "Trinidad and Tobago"),]$Location <- "Trinidad & Tobago"
pop[which(pop$Location == "Turks and Caicos Islands"),]$Location <- "Turks & Caicos Islands"
pop[which(pop$Location == "United Republic of Tanzania"),]$Location <- "Tanzania"
pop[which(pop$Location == "United States of America"),]$Location <- "United States"
pop[which(pop$Location == "United States Virgin Islands"),]$Location <- "U.S. Virgin Islands"
pop[which(pop$Location == "Venezuela (Bolivarian Republic of)"),]$Location <- "Venezuela"
pop[which(pop$Location == "Viet Nam"),]$Location <- "Vietnam"
pop[which(pop$Location == "Wallis and Futuna Islands"),]$Location <- "Wallis & Futuna"
pop <- pop[order(pop$Location),] # ordering alphabetically
colnames(pop) <- c("name","year0","pop") # renaming columns
pop$pop <- pop$pop*1e03

pop.wide <- reshape(data = pop,
                    varying = NULL,
                    v.names = NULL,
                    timevar = "year0",
                    idvar = c("name"),
                    direction = "wide",
                    new.row.names = NULL,
                    sep = ".")

## gdp
gdp <- as.data.frame(read.csv("data/GDP_Global_1990_2021.csv"))[,-c(2:4,36)]
colnames(gdp) <- c("name",1990:2020)

unique(c(gdp[!c(gdp$name,gdp$name)%in%States$name,]$name,
         gdp[!c(gdp$name,gdp$name)%in%States$name,]$name))

unique(States[!States$name%in%gdp$name,]$name)

gdp[which(gdp$name == "Antigua and Barbuda"),]$name <- "Antigua & Barbuda"
gdp[which(gdp$name == "Bahamas, The"),]$name <- "Bahamas"
gdp[which(gdp$name == "Bosnia and Herzegovina"),]$name <- "Bosnia & Herzegovina"
gdp[which(gdp$name == "Brunei Darussalam"),]$name <- "Brunei"
gdp[which(gdp$name == "Cabo Verde"),]$name <- "Cape Verde"
gdp[which(gdp$name == "Congo, Dem. Rep."),]$name <- "DR Congo"
gdp[which(gdp$name == "Congo, Rep."),]$name <- "Congo"
gdp[which(gdp$name == "Cote d'Ivoire"),]$name <- "Côte d’Ivoire"
gdp[which(gdp$name == "Curacao"),]$name <- "Curaçao"
gdp[which(gdp$name == "Czech Republic"),]$name <- "Czechia"
gdp[which(gdp$name == "Egypt, Arab Rep."),]$name <- "Egypt"
gdp[which(gdp$name == "Gambia, The"),]$name <- "Gambia"
gdp[which(gdp$name == "Hong Kong SAR, China"),]$name <- "Hong Kong SAR"
gdp[which(gdp$name == "Iran, Islamic Rep."),]$name <- "Iran"
gdp[which(gdp$name == "Korea, Dem. People's Rep."),]$name <- "North Korea"
gdp[which(gdp$name == "Korea, Rep."),]$name <- "South Korea"
# gdp[which(gdp$name == "Kosovo"),]$name
gdp[which(gdp$name == "Kyrgyz Republic"),]$name <- "Kyrgyzstan"
gdp[which(gdp$name == "Lao PDR"),]$name <- "Laos"
gdp[which(gdp$name == "Macao SAR, China"),]$name <- "Macao SAR"
gdp[which(gdp$name == "Micronesia, Fed. Sts."),]$name <- "FS Micronesia"
gdp[which(gdp$name == "Russian Federation"),]$name <- "Russia"
gdp[which(gdp$name == "Sao Tome and Principe"),]$name <- "São Tomé & Príncipe"
gdp[which(gdp$name == "Sint Maarten (Dutch part)"),]$name <- "Sint Maarten"
gdp[which(gdp$name == "Slovak Republic"),]$name <- "Slovakia"
gdp[which(gdp$name == "St. Kitts and Nevis"),]$name <- "St. Kitts & Nevis"
# gdp[which(gdp$name == "St. Martin (French part)"),]$name
gdp[which(gdp$name == "St. Vincent and the Grenadines"),]$name <- "St. Vincent & Grenadines"
gdp[which(gdp$name == "Syrian Arab Republic"),]$name <- "Syria"
gdp[which(gdp$name == "Trinidad and Tobago"),]$name <- "Trinidad & Tobago"
gdp[which(gdp$name == "Turks and Caicos Islands"),]$name <- "Turks & Caicos Islands"
gdp[which(gdp$name == "Venezuela, RB"),]$name <- "Venezuela"
gdp[which(gdp$name == "Virgin Islands (U.S.)"),]$name <- "U.S. Virgin Islands"
gdp[which(gdp$name == "Yemen, Rep."),]$name <- "Yemen"

gdp[gdp=='..'] <- NA
gdp <- gdp[!gdp$name %in% unique(gdp[!gdp$name %in% States$name,]$name),]
gdp.long <- reshape(data = as.data.frame(gdp[,1:32]),
                    varying = 2:32,
                    v.names = "gdp",
                    idvar = c("name"),
                    direction = "long",
                    sep = "")
for (i in 1:31){
  gdp.long[gdp.long$time==i,]$time <- i+1989
}
gdp.long$gdp <- as.numeric(gdp.long$gdp)
# convert back to wide format
gdp.wide <- reshape(data = gdp.long,
                    varying = NULL,
                    v.names = NULL,
                    timevar = "time",
                    idvar = c("name"),
                    direction = "wide",
                    new.row.names = NULL,
                    sep = ".")

### Call in SPEI Data
# MAX
spei.max <- read.csv("data/spei01-annual-max.csv")[,-1]
colnames(spei.max)[1] <- "COUNTRY"
unique(spei.max[!spei.max$COUNTRY %in% States$name,]$COUNTRY)
spei.max[which(spei.max$COUNTRY == "Bosnia and Herzegovina"),]$COUNTRY <- "Bosnia & Herzegovina"
spei.max[which(spei.max$COUNTRY == "Czech Republic"),]$COUNTRY <- "Czechia"
spei.max[which(spei.max$COUNTRY == "Democratic Republic of the Congo"),]$COUNTRY <- "DR Congo"
spei.max[which(spei.max$COUNTRY == "East Timor"),]$COUNTRY <- "Timor-Leste"
spei.max[which(spei.max$COUNTRY == "Guinea Bissau"),]$COUNTRY <- "Guinea-Bissau"
spei.max[which(spei.max$COUNTRY == "Ivory Coast"),]$COUNTRY <- "Côte d’Ivoire"
spei.max[which(spei.max$COUNTRY == "Macedonia"),]$COUNTRY <- "North Macedonia"
# spei.max[which(spei.max$COUNTRY == "Northern Cyprus"),]$COUNTRY <- "Cyprus"
spei.max[which(spei.max$COUNTRY == "Republic of Serbia"),]$COUNTRY <- "Serbia"
spei.max[which(spei.max$COUNTRY == "Republic of the Congo"),]$COUNTRY <- "Congo"
spei.max[which(spei.max$COUNTRY == "Saint Vincent and the Grenadines"),]$COUNTRY <- "St. Vincent & Grenadines"
spei.max[which(spei.max$COUNTRY == "Sao Tome and Principe"),]$COUNTRY <- "São Tomé & Príncipe"
spei.max[which(spei.max$COUNTRY == "Swaziland"),]$COUNTRY <- "Eswatini"
spei.max[which(spei.max$COUNTRY == "The Bahamas"),]$COUNTRY <- "Bahamas"
spei.max[which(spei.max$COUNTRY == "Trinidad and Tobago"),]$COUNTRY <- "Trinidad & Tobago"
spei.max[which(spei.max$COUNTRY == "United Republic of Tanzania"),]$COUNTRY <- "Tanzania"
spei.max[which(spei.max$COUNTRY == "United States of America"),]$COUNTRY <- "United States"
spei.max[which(spei.max$COUNTRY == "United States Virgin Islands"),]$COUNTRY <- "U.S. Virgin Islands"

# MIN
spei.min <- read.csv("data/spei01-annual-min.csv")[,-1]
colnames(spei.min)[1] <- "COUNTRY"
unique(spei.min[!spei.min$COUNTRY %in% States$statenme,]$COUNTRY)
spei.min[which(spei.min$COUNTRY == "Bosnia and Herzegovina"),]$COUNTRY <- "Bosnia & Herzegovina"
spei.min[which(spei.min$COUNTRY == "Czech Republic"),]$COUNTRY <- "Czechia"
spei.min[which(spei.min$COUNTRY == "Democratic Republic of the Congo"),]$COUNTRY <- "DR Congo"
spei.min[which(spei.min$COUNTRY == "East Timor"),]$COUNTRY <- "Timor-Leste"
spei.min[which(spei.min$COUNTRY == "Guinea Bissau"),]$COUNTRY <- "Guinea-Bissau"
spei.min[which(spei.min$COUNTRY == "Ivory Coast"),]$COUNTRY <- "Côte d’Ivoire"
spei.min[which(spei.min$COUNTRY == "Macedonia"),]$COUNTRY <- "North Macedonia"
# spei.min[which(spei.min$COUNTRY == "Northern Cyprus"),]$COUNTRY <- "Cyprus"
spei.min[which(spei.min$COUNTRY == "Republic of Serbia"),]$COUNTRY <- "Serbia"
spei.min[which(spei.min$COUNTRY == "Republic of the Congo"),]$COUNTRY <- "Congo"
spei.min[which(spei.min$COUNTRY == "Saint Vincent and the Grenadines"),]$COUNTRY <- "St. Vincent & Grenadines"
spei.min[which(spei.min$COUNTRY == "Sao Tome and Principe"),]$COUNTRY <- "São Tomé & Príncipe"
spei.min[which(spei.min$COUNTRY == "Swaziland"),]$COUNTRY <- "Eswatini"
spei.min[which(spei.min$COUNTRY == "The Bahamas"),]$COUNTRY <- "Bahamas"
spei.min[which(spei.min$COUNTRY == "Trinidad and Tobago"),]$COUNTRY <- "Trinidad & Tobago"
spei.min[which(spei.min$COUNTRY == "United Republic of Tanzania"),]$COUNTRY <- "Tanzania"
spei.min[which(spei.min$COUNTRY == "United States of America"),]$COUNTRY <- "United States"
spei.min[which(spei.min$COUNTRY == "United States Virgin Islands"),]$COUNTRY <- "U.S. Virgin Islands"

# MEAN
spei.mean <- read.csv("data/spei01-annual-mean.csv")[,-1]
colnames(spei.mean)[1] <- "COUNTRY"
unique(spei.mean[!spei.mean$COUNTRY %in% States$statenme,]$COUNTRY)
spei.mean[which(spei.mean$COUNTRY == "Bosnia and Herzegovina"),]$COUNTRY <- "Bosnia & Herzegovina"
spei.mean[which(spei.mean$COUNTRY == "Czech Republic"),]$COUNTRY <- "Czechia"
spei.mean[which(spei.mean$COUNTRY == "Democratic Republic of the Congo"),]$COUNTRY <- "DR Congo"
spei.mean[which(spei.mean$COUNTRY == "East Timor"),]$COUNTRY <- "Timor-Leste"
spei.mean[which(spei.mean$COUNTRY == "Guinea Bissau"),]$COUNTRY <- "Guinea-Bissau"
spei.mean[which(spei.mean$COUNTRY == "Ivory Coast"),]$COUNTRY <- "Côte d’Ivoire"
spei.mean[which(spei.mean$COUNTRY == "Macedonia"),]$COUNTRY <- "North Macedonia"
# spei.mean[which(spei.mean$COUNTRY == "Northern Cyprus"),]$COUNTRY <- "Cyprus"
spei.mean[which(spei.mean$COUNTRY == "Republic of Serbia"),]$COUNTRY <- "Serbia"
spei.mean[which(spei.mean$COUNTRY == "Republic of the Congo"),]$COUNTRY <- "Congo"
spei.mean[which(spei.mean$COUNTRY == "Saint Vincent and the Grenadines"),]$COUNTRY <- "St. Vincent & Grenadines"
spei.mean[which(spei.mean$COUNTRY == "Sao Tome and Principe"),]$COUNTRY <- "São Tomé & Príncipe"
spei.mean[which(spei.mean$COUNTRY == "Swaziland"),]$COUNTRY <- "Eswatini"
spei.mean[which(spei.mean$COUNTRY == "The Bahamas"),]$COUNTRY <- "Bahamas"
spei.mean[which(spei.mean$COUNTRY == "Trinidad and Tobago"),]$COUNTRY <- "Trinidad & Tobago"
spei.mean[which(spei.mean$COUNTRY == "United Republic of Tanzania"),]$COUNTRY <- "Tanzania"
spei.mean[which(spei.mean$COUNTRY == "United States of America"),]$COUNTRY <- "United States"
spei.mean[which(spei.mean$COUNTRY == "United States Virgin Islands"),]$COUNTRY <- "U.S. Virgin Islands"

# MAX.cv
spei.max.cv <- read.csv("data/spei01-annual-max-cv.csv")[,-1]
colnames(spei.max.cv)[1] <- "COUNTRY"
unique(spei.max.cv[!spei.max.cv$COUNTRY %in% States$statenme,]$COUNTRY)
spei.max.cv[which(spei.max.cv$COUNTRY == "Bosnia and Herzegovina"),]$COUNTRY <- "Bosnia & Herzegovina"
spei.max.cv[which(spei.max.cv$COUNTRY == "Czech Republic"),]$COUNTRY <- "Czechia"
spei.max.cv[which(spei.max.cv$COUNTRY == "Democratic Republic of the Congo"),]$COUNTRY <- "DR Congo"
spei.max.cv[which(spei.max.cv$COUNTRY == "East Timor"),]$COUNTRY <- "Timor-Leste"
spei.max.cv[which(spei.max.cv$COUNTRY == "Guinea Bissau"),]$COUNTRY <- "Guinea-Bissau"
spei.max.cv[which(spei.max.cv$COUNTRY == "Ivory Coast"),]$COUNTRY <- "Côte d’Ivoire"
spei.max.cv[which(spei.max.cv$COUNTRY == "Macedonia"),]$COUNTRY <- "North Macedonia"
# spei.max.cv[which(spei.max.cv$COUNTRY == "Northern Cyprus"),]$COUNTRY <- "Cyprus"
spei.max.cv[which(spei.max.cv$COUNTRY == "Republic of Serbia"),]$COUNTRY <- "Serbia"
spei.max.cv[which(spei.max.cv$COUNTRY == "Republic of the Congo"),]$COUNTRY <- "Congo"
spei.max.cv[which(spei.max.cv$COUNTRY == "Saint Vincent and the Grenadines"),]$COUNTRY <- "St. Vincent & Grenadines"
spei.max.cv[which(spei.max.cv$COUNTRY == "Sao Tome and Principe"),]$COUNTRY <- "São Tomé & Príncipe"
spei.max.cv[which(spei.max.cv$COUNTRY == "Swaziland"),]$COUNTRY <- "Eswatini"
spei.max.cv[which(spei.max.cv$COUNTRY == "The Bahamas"),]$COUNTRY <- "Bahamas"
spei.max.cv[which(spei.max.cv$COUNTRY == "Trinidad and Tobago"),]$COUNTRY <- "Trinidad & Tobago"
spei.max.cv[which(spei.max.cv$COUNTRY == "United Republic of Tanzania"),]$COUNTRY <- "Tanzania"
spei.max.cv[which(spei.max.cv$COUNTRY == "United States of America"),]$COUNTRY <- "United States"
spei.max.cv[which(spei.max.cv$COUNTRY == "United States Virgin Islands"),]$COUNTRY <- "U.S. Virgin Islands"

# MIN.cv
spei.min.cv <- read.csv("data/spei01-annual-min-cv.csv")[,-1]
colnames(spei.min.cv)[1] <- "COUNTRY"
unique(spei.min.cv[!spei.min.cv$COUNTRY %in% States$statenme,]$COUNTRY)
spei.min.cv[which(spei.min.cv$COUNTRY == "Bosnia and Herzegovina"),]$COUNTRY <- "Bosnia & Herzegovina"
spei.min.cv[which(spei.min.cv$COUNTRY == "Czech Republic"),]$COUNTRY <- "Czechia"
spei.min.cv[which(spei.min.cv$COUNTRY == "Democratic Republic of the Congo"),]$COUNTRY <- "DR Congo"
spei.min.cv[which(spei.min.cv$COUNTRY == "East Timor"),]$COUNTRY <- "Timor-Leste"
spei.min.cv[which(spei.min.cv$COUNTRY == "Guinea Bissau"),]$COUNTRY <- "Guinea-Bissau"
spei.min.cv[which(spei.min.cv$COUNTRY == "Ivory Coast"),]$COUNTRY <- "Côte d’Ivoire"
spei.min.cv[which(spei.min.cv$COUNTRY == "Macedonia"),]$COUNTRY <- "North Macedonia"
# spei.min.cv[which(spei.min.cv$COUNTRY == "Northern Cyprus"),]$COUNTRY <- "Cyprus"
spei.min.cv[which(spei.min.cv$COUNTRY == "Republic of Serbia"),]$COUNTRY <- "Serbia"
spei.min.cv[which(spei.min.cv$COUNTRY == "Republic of the Congo"),]$COUNTRY <- "Congo"
spei.min.cv[which(spei.min.cv$COUNTRY == "Saint Vincent and the Grenadines"),]$COUNTRY <- "St. Vincent & Grenadines"
spei.min.cv[which(spei.min.cv$COUNTRY == "Sao Tome and Principe"),]$COUNTRY <- "São Tomé & Príncipe"
spei.min.cv[which(spei.min.cv$COUNTRY == "Swaziland"),]$COUNTRY <- "Eswatini"
spei.min.cv[which(spei.min.cv$COUNTRY == "The Bahamas"),]$COUNTRY <- "Bahamas"
spei.min.cv[which(spei.min.cv$COUNTRY == "Trinidad and Tobago"),]$COUNTRY <- "Trinidad & Tobago"
spei.min.cv[which(spei.min.cv$COUNTRY == "United Republic of Tanzania"),]$COUNTRY <- "Tanzania"
spei.min.cv[which(spei.min.cv$COUNTRY == "United States of America"),]$COUNTRY <- "United States"
spei.min.cv[which(spei.min.cv$COUNTRY == "United States Virgin Islands"),]$COUNTRY <- "U.S. Virgin Islands"

# MEAN.cv
spei.mean.cv <- read.csv("data/spei01-annual-mean-cv.csv")[,-1]
colnames(spei.mean.cv)[1] <- "COUNTRY"
unique(spei.mean.cv[!spei.mean.cv$COUNTRY %in% States$statenme,]$COUNTRY)
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Bosnia and Herzegovina"),]$COUNTRY <- "Bosnia & Herzegovina"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Czech Republic"),]$COUNTRY <- "Czechia"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Democratic Republic of the Congo"),]$COUNTRY <- "DR Congo"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "East Timor"),]$COUNTRY <- "Timor-Leste"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Guinea Bissau"),]$COUNTRY <- "Guinea-Bissau"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Ivory Coast"),]$COUNTRY <- "Côte d’Ivoire"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Macedonia"),]$COUNTRY <- "North Macedonia"
# spei.mean.cv[which(spei.mean.cv$COUNTRY == "Northern Cyprus"),]$COUNTRY <- "Cyprus"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Republic of Serbia"),]$COUNTRY <- "Serbia"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Republic of the Congo"),]$COUNTRY <- "Congo"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Saint Vincent and the Grenadines"),]$COUNTRY <- "St. Vincent & Grenadines"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Sao Tome and Principe"),]$COUNTRY <- "São Tomé & Príncipe"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Swaziland"),]$COUNTRY <- "Eswatini"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "The Bahamas"),]$COUNTRY <- "Bahamas"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "Trinidad and Tobago"),]$COUNTRY <- "Trinidad & Tobago"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "United Republic of Tanzania"),]$COUNTRY <- "Tanzania"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "United States of America"),]$COUNTRY <- "United States"
spei.mean.cv[which(spei.mean.cv$COUNTRY == "United States Virgin Islands"),]$COUNTRY <- "U.S. Virgin Islands"


### Call in Uppsala Conflict Data
conflict <- read.csv("data/GEDEvent_v22_1.csv")
# conflict <- read.csv("data/Uppsala_Georeferenced_Conflict_Dataset.csv")
conflict <- conflict[c("country","country_id","year","best")]
unique(conflict[!conflict$country %in% States$name,]$country)
conflict[which(conflict$country == "Bosnia-Herzegovina"),]$country <- "Bosnia & Herzegovina"
conflict[which(conflict$country == "Cambodia (Kampuchea)"),]$country <- "Cambodia"
conflict[which(conflict$country == "DR Congo (Zaire)"),]$country <- "DR Congo"
conflict[which(conflict$country == "Ivory Coast"),]$country <- "Côte d’Ivoire"
conflict[which(conflict$country == "Kingdom of eSwatini (Swaziland)"),]$country <- "Eswatini"
conflict[which(conflict$country == "Macedonia, FYR"),]$country <- "North Macedonia"
conflict[which(conflict$country == "Madagascar (Malagasy)"),]$country <- "Madagascar"
conflict[which(conflict$country == "Myanmar (Burma)"),]$country <- "Myanmar"
conflict[which(conflict$country == "Russia (Soviet Union)"),]$country <- "Russia"
conflict[which(conflict$country == "Serbia (Yugoslavia)"),]$country <- "Serbia"
conflict[which(conflict$country == "Trinidad and Tobago"),]$country <- "Trinidad & Tobago"
conflict[which(conflict$country == "United States of America"),]$country <- "United States"
conflict[which(conflict$country == "Yemen (North Yemen)"),]$country <- "Yemen"
conflict[which(conflict$country == "Zimbabwe (Rhodesia)"),]$country <- "Zimbabwe"

nConflicts <- as.data.frame(xtabs(~ country + year, data = conflict))
nConflicts$country <- as.character(nConflicts$country)
nConflicts$year <- as.character(nConflicts$year)
colnames(nConflicts)[3] <- "nConflicts"
conflicts.dead <- aggregate(conflict$best, by = list(country = conflict$country, year = conflict$year), sum)
colnames(conflicts.dead)[3] = "conflict.sum.dead"
conflict.dat <- merge(nConflicts, conflicts.dead, by.x = c("country","year"), by.y = c("country","year"), all = TRUE)
conflict.dat <- conflict.dat[conflict.dat$year>=1990 & conflict.dat$year<=2020,]

# convert to wide format
conflict.wide <- reshape(data = conflict.dat,
                         varying = NULL,
                         v.names = NULL,
                         timevar = "year",
                         idvar = c("country"),
                         direction = "wide",
                         new.row.names = NULL,
                         sep = ".")
conflict.wide.ordered <- cbind(conflict.wide$country, conflict.wide[,2:length(colnames(conflict.wide))][,order(colnames(conflict.wide[2:length(colnames(conflict.wide))]))])
colnames(conflict.wide.ordered)[1] <- "country"

### Em-Dat Disaster Data
disaster <- read_xlsx("data/emdat_public_1990-2020.xlsx")[-c(1:6),c(2,7,11,35,39)]
colnames(disaster) <- c("Year","Disaster.Type","Country","Total.Deaths","Total.Affected")
unique(disaster$Country[!disaster$Country%in%States$name])
unique(States$name[!States$name%in%disaster$Country])
disaster[which(disaster$Country == "Bolivia (Plurinational State of)"),]$Country <- "Bolivia"
disaster[which(disaster$Country == "Antigua and Barbuda"),]$Country <- "Antigua & Barbuda"
disaster[which(disaster$Country == "Bahamas (the)"),]$Country <- "Bahamas"
disaster[which(disaster$Country == "Congo (the)"),]$Country <- "Congo"
disaster[which(disaster$Country == "Cook Islands (the)"),]$Country <- "Cook Islands"
disaster[which(disaster$Country == "Germany Fed Rep"),]$Country <- "Germany"
disaster[which(disaster$Country == "Dominican Republic (the)"),]$Country <- "Dominican Republic"
disaster[which(disaster$Country == "Comoros (the)"),]$Country <- "Comoros"
disaster[which(disaster$Country == "Cabo Verde"),]$Country <- "Cape Verde"
disaster[which(disaster$Country == "United Arab Emirates (the)"),]$Country <- "United Arab Emirates"
disaster[which(disaster$Country == "Iran (Islamic Republic of)"),]$Country <- "Iran"
disaster[which(disaster$Country == "Philippines (the)"),]$Country <- "Philippines"
disaster[which(disaster$Country == "Niger (the)"),]$Country <- "Niger"
disaster[which(disaster$Country == "United Kingdom of Great Britain and Northern Ireland (the)"),]$Country <- "United Kingdom"
disaster[which(disaster$Country == "Saint Kitts and Nevis"),]$Country <- "St. Kitts & Nevis"
disaster[which(disaster$Country == "Korea (the Republic of)"),]$Country <- "South Korea"
disaster[which(disaster$Country == "Netherlands (the)"),]$Country <- "Netherlands"
disaster[which(disaster$Country == "Lao People's Democratic Republic (the)"),]$Country <- "Laos"
disaster[which(disaster$Country == "Marshall Islands (the)"),]$Country <- "Marshall Islands"
disaster[which(disaster$Country == "Korea (the Democratic People's Republic of)"),]$Country <- "North Korea"
disaster[which(disaster$Country == "Russian Federation (the)"),]$Country <- "Russia"
disaster[which(disaster$Country == "Gambia (the)"),]$Country <- "Gambia"
disaster[which(disaster$Country == "Macedonia (the former Yugoslav Republic of)"),]$Country <- "North Macedonia"
disaster[which(disaster$Country == "United States of America (the)"),]$Country <- "United States"
disaster[which(disaster$Country == "Sudan (the)"),]$Country <- "Sudan"
disaster[which(disaster$Country == "Soviet Union"),]$Country <- "Russia"
# disaster[which(disaster$Country == "Taiwan (Province of China)"),]$Country <- "Taiwan"
disaster[which(disaster$Country == "Tanzania, United Republic of"),]$Country <- "Tanzania"
disaster[which(disaster$Country == "Venezuela (Bolivarian Republic of)"),]$Country <- "Venezuela"
disaster[which(disaster$Country == "Viet Nam"),]$Country <- "Vietnam"
disaster[which(disaster$Country == "Yugoslavia"),]$Country <- "Serbia & Montenegro"
disaster[which(disaster$Country == "Congo (the Democratic Republic of the)"),]$Country <- "DR Congo"
disaster[which(disaster$Country == "Serbia Montenegro"),]$Country <- "Serbia & Montenegro"
disaster[which(disaster$Country == "Saint Vincent and the Grenadines"),]$Country <- "St. Vincent & Grenadines"
disaster[which(disaster$Country == "Brunei Darussalam"),]$Country <- "Brunei"
disaster[which(disaster$Country == "Moldova (the Republic of)"),]$Country <- "Moldova"
disaster[which(disaster$Country == "Saint Lucia"),]$Country <- "St. Lucia"
disaster[which(disaster$Country == "Czech Republic (the)"),]$Country <- "Czechia"
disaster[which(disaster$Country == "Syrian Arab Republic"),]$Country <- "Syria"
disaster[which(disaster$Country == "Micronesia (Federated States of)"),]$Country <- "FS Micronesia"
disaster[which(disaster$Country == "Hong Kong"),]$Country <- "Hong Kong SAR"
disaster[which(disaster$Country == "Macao"),]$Country <- "Macao SAR"
disaster[which(disaster$Country == "Turks and Caicos Islands (the)"),]$Country <- "Turks & Caicos Islands"
disaster[which(disaster$Country == "Trinidad and Tobago"),]$Country <- "Trinidad & Tobago"
disaster[which(disaster$Country == "Wallis and Futuna"),]$Country <- "Wallis & Futuna"
disaster[which(disaster$Country == "Bosnia and Herzegovina"),]$Country <- "Bosnia & Herzegovina"
disaster[which(disaster$Country == "Virgin Island (U.S.)"),]$Country <- "U.S. Virgin Islands"
disaster[which(disaster$Country == "Virgin Island (British)"),]$Country <- "British Virgin Islands"
disaster[which(disaster$Country == "Cayman Islands (the)"),]$Country <- "Cayman Islands"
disaster[which(disaster$Country == "Saint Helena, Ascension and Tristan da Cunha"),]$Country <- "St. Helena"
disaster[which(disaster$Country == "Sao Tome and Principe"),]$Country <- "São Tomé & Príncipe"
disaster[which(disaster$Country == "Sint Maarten (Dutch part)"),]$Country <- "Sint Maarten"
disaster[which(disaster$Country == "Palestine, State of"),]$Country <- "Palestine"
disaster$Total.Deaths <- as.numeric(disaster$Total.Deaths)
disaster$Total.Affected <- as.numeric(disaster$Total.Affected)
disaster$Year <- as.numeric(disaster$Year)

emdat_storm <- disaster[disaster$Disaster.Type=="Storm",c(1,3:5)]
colnames(emdat_storm)[3:4] <- paste("storm",colnames(emdat_storm)[3:4],sep = ".")
emdat_storm <- as.data.frame(emdat_storm) %>% 
  group_by(Country, Year) %>% 
  summarise(storm.Total.Deaths = sum(storm.Total.Deaths, na.rm = T),
            storm.Total.Affected = sum(storm.Total.Affected, na.rm = T))
emdat_storm_wide <- reshape(data = as.data.frame(emdat_storm),
                            varying = NULL,
                            v.names = NULL,
                            timevar = "Year",
                            idvar = "Country",
                            direction = "wide",
                            new.row.names = NULL,
                            sep = ".")
emdat_storm_wide <- cbind(emdat_storm_wide$Country, emdat_storm_wide[,2:length(colnames(emdat_storm_wide))][,order(colnames(emdat_storm_wide[2:length(colnames(emdat_storm_wide))]))])
colnames(emdat_storm_wide)[1] <- "Country"

emdat_flood <- disaster[disaster$Disaster.Type=="Flood",c(1,3:5)]
colnames(emdat_flood)[3:4] <- paste("flood",colnames(emdat_flood)[3:4],sep = ".")
emdat_flood <- as.data.frame(emdat_flood) %>% 
  group_by(Country, Year) %>% 
  summarise(flood.Total.Deaths = sum(flood.Total.Deaths, na.rm = T),
            flood.Total.Affected = sum(flood.Total.Affected, na.rm = T))
emdat_flood_wide <- reshape(data = as.data.frame(emdat_flood),
                            varying = NULL,
                            v.names = NULL,
                            timevar = "Year",
                            idvar = "Country",
                            direction = "wide",
                            new.row.names = NULL,
                            sep = ".")
emdat_flood_wide <- cbind(emdat_flood_wide$Country, emdat_flood_wide[,2:length(colnames(emdat_flood_wide))][,order(colnames(emdat_flood_wide[2:length(colnames(emdat_flood_wide))]))])
colnames(emdat_flood_wide)[1] <- "Country"


### Calculate 5-year maxes, minimums, and moving averages
# pop delta
pop.delta <- pop.wide
colnames(pop.delta) <- c("Name","1995","2000","2005","2010","2015","2020")
pop.delta.long <- reshape(data = as.data.frame(pop.delta),
                          varying = 2:7,
                          v.names = "pop",
                          idvar = c("Name"),
                          direction = "long",
                          new.row.names = NULL,
                          sep = "")
pop.delta.long[pop.delta.long$time==1,]$time <- "1995"
pop.delta.long[pop.delta.long$time==2,]$time <- "2000"
pop.delta.long[pop.delta.long$time==3,]$time <- "2005"
pop.delta.long[pop.delta.long$time==4,]$time <- "2010"
pop.delta.long[pop.delta.long$time==5,]$time <- "2015"
pop.delta.long[pop.delta.long$time==6,]$time <- "2020"

# gdp delta
gdp.delta <- data.frame(Name = gdp.wide$name,
                        "1995" = rowMeans(gdp.wide[,3:7], na.rm = T),
                        "2000" = rowMeans(gdp.wide[,8:12], na.rm = T),
                        "2005" = rowMeans(gdp.wide[,13:17], na.rm = T),
                        "2010" = rowMeans(gdp.wide[,18:22], na.rm = T),
                        "2015" = rowMeans(gdp.wide[,23:27], na.rm = T),
                        "2020" = rowMeans(gdp.wide[,28:32], na.rm = T))
gdp.delta[gdp.delta=="NaN"] <- NA
gdp.delta.long <- reshape(data = as.data.frame(gdp.delta),
                          varying = 2:7,
                          v.names = "gdp",
                          idvar = c("Name"),
                          direction = "long",
                          new.row.names = NULL,
                          sep = "")
gdp.delta.long[gdp.delta.long$time==1,]$time <- "1995"
gdp.delta.long[gdp.delta.long$time==2,]$time <- "2000"
gdp.delta.long[gdp.delta.long$time==3,]$time <- "2005"
gdp.delta.long[gdp.delta.long$time==4,]$time <- "2010"
gdp.delta.long[gdp.delta.long$time==5,]$time <- "2015"
gdp.delta.long[gdp.delta.long$time==6,]$time <- "2020"

# spei.max delta
spei.max.delta <- data.frame(Name = spei.max$COUNTRY,
                             "1995" = rowMaxs(as.matrix(spei.max[,3:7]), na.rm = T),
                             "2000" = rowMaxs(as.matrix(spei.max[,8:12]), na.rm = T),
                             "2005" = rowMaxs(as.matrix(spei.max[,13:17]), na.rm = T),
                             "2010" = rowMaxs(as.matrix(spei.max[,18:22]), na.rm = T),
                             "2015" = rowMaxs(as.matrix(spei.max[,23:27]), na.rm = T),
                             "2020" = rowMaxs(as.matrix(spei.max[,28:31]), na.rm = T))
spei.max.delta.long <- reshape(data = as.data.frame(spei.max.delta),
                               varying = 2:7,
                               v.names = "spei.max",
                               idvar = c("Name"),
                               direction = "long",
                               new.row.names = NULL,
                               sep = ".")
spei.max.delta.long[spei.max.delta.long$time==1,]$time <- "1995"
spei.max.delta.long[spei.max.delta.long$time==2,]$time <- "2000"
spei.max.delta.long[spei.max.delta.long$time==3,]$time <- "2005"
spei.max.delta.long[spei.max.delta.long$time==4,]$time <- "2010"
spei.max.delta.long[spei.max.delta.long$time==5,]$time <- "2015"
spei.max.delta.long[spei.max.delta.long$time==6,]$time <- "2020"
# spei.min delta
spei.min.delta <- data.frame(Name = spei.min$COUNTRY,
                             "1995" = rowMins(as.matrix(spei.min[,3:7]), na.rm = T),
                             "2000" = rowMins(as.matrix(spei.min[,8:12]), na.rm = T),
                             "2005" = rowMins(as.matrix(spei.min[,13:17]), na.rm = T),
                             "2010" = rowMins(as.matrix(spei.min[,18:22]), na.rm = T),
                             "2015" = rowMins(as.matrix(spei.min[,23:27]), na.rm = T),
                             "2020" = rowMins(as.matrix(spei.min[,28:32]), na.rm = T))
spei.min.delta.long <- reshape(data = as.data.frame(spei.min.delta),
                               varying = 2:7,
                               v.names = "spei.min",
                               idvar = c("Name"),
                               direction = "long",
                               new.row.names = NULL,
                               sep = "")
spei.min.delta.long[spei.min.delta.long$time==1,]$time <- "1995"
spei.min.delta.long[spei.min.delta.long$time==2,]$time <- "2000"
spei.min.delta.long[spei.min.delta.long$time==3,]$time <- "2005"
spei.min.delta.long[spei.min.delta.long$time==4,]$time <- "2010"
spei.min.delta.long[spei.min.delta.long$time==5,]$time <- "2015"
spei.min.delta.long[spei.min.delta.long$time==6,]$time <- "2020"
# spei.mean delta
spei.mean.delta <- data.frame(Name = spei.mean$COUNTRY,
                              "1995" = rowMeans(spei.mean[,3:7], na.rm = T),
                              "2000" = rowMeans(spei.mean[,8:12], na.rm = T), 
                              "2005" = rowMeans(spei.mean[,13:17], na.rm = T), 
                              "2010" = rowMeans(spei.mean[,18:22], na.rm = T), 
                              "2015" = rowMeans(spei.mean[,23:27], na.rm = T), 
                              "2020" = rowMeans(spei.mean[,28:32], na.rm = T))
spei.mean.delta.long <- reshape(data = as.data.frame(spei.mean.delta),
                                varying = 2:7,
                                v.names = "spei.mean",
                                idvar = c("Name"),
                                direction = "long",
                                new.row.names = NULL,
                                sep = "")
spei.mean.delta.long[spei.mean.delta.long$time==1,]$time <- "1995"
spei.mean.delta.long[spei.mean.delta.long$time==2,]$time <- "2000"
spei.mean.delta.long[spei.mean.delta.long$time==3,]$time <- "2005"
spei.mean.delta.long[spei.mean.delta.long$time==4,]$time <- "2010"
spei.mean.delta.long[spei.mean.delta.long$time==5,]$time <- "2015"
spei.mean.delta.long[spei.mean.delta.long$time==6,]$time <- "2020"
# spei.max.cv delta
spei.max.cv.delta <- data.frame(Name = spei.max.cv$COUNTRY,
                                "1995" = rowMaxs(as.matrix(spei.max.cv[,3:7], na.rm = T)),
                                "2000" = rowMaxs(as.matrix(spei.max.cv[,8:12], na.rm = T)),
                                "2005" = rowMaxs(as.matrix(spei.max.cv[,13:17], na.rm = T)),
                                "2010" = rowMaxs(as.matrix(spei.max.cv[,18:22], na.rm = T)),
                                "2015" = rowMaxs(as.matrix(spei.max.cv[,23:27], na.rm = T)),
                                "2020" = rowMaxs(as.matrix(spei.max.cv[,28:32], na.rm = T)))
spei.max.cv.delta.long <- reshape(data = as.data.frame(spei.max.cv.delta),
                                  varying = 2:7,
                                  v.names = "spei.max.cv",
                                  idvar = c("Name"),
                                  direction = "long",
                                  new.row.names = NULL,
                                  sep = "")
spei.max.cv.delta.long[spei.max.cv.delta.long$time==1,]$time <- "1995"
spei.max.cv.delta.long[spei.max.cv.delta.long$time==2,]$time <- "2000"
spei.max.cv.delta.long[spei.max.cv.delta.long$time==3,]$time <- "2005"
spei.max.cv.delta.long[spei.max.cv.delta.long$time==4,]$time <- "2010"
spei.max.cv.delta.long[spei.max.cv.delta.long$time==5,]$time <- "2015"
spei.max.cv.delta.long[spei.max.cv.delta.long$time==6,]$time <- "2020"
# spei.min.cv delta
spei.min.cv.delta <- data.frame(Name = spei.min.cv$COUNTRY,
                                "1995" = rowMins(as.matrix(spei.min.cv[,3:7], na.rm = T)),
                                "2000" = rowMins(as.matrix(spei.min.cv[,8:12], na.rm = T)),
                                "2005" = rowMins(as.matrix(spei.min.cv[,13:17], na.rm = T)),
                                "2010" = rowMins(as.matrix(spei.min.cv[,18:22], na.rm = T)),
                                "2015" = rowMins(as.matrix(spei.min.cv[,23:27], na.rm = T)),
                                "2020" = rowMins(as.matrix(spei.min.cv[,28:32], na.rm = T)))
spei.min.cv.delta.long <- reshape(data = as.data.frame(spei.min.cv.delta),
                                  varying = 2:7,
                                  v.names = "spei.min.cv",
                                  idvar = c("Name"),
                                  direction = "long",
                                  new.row.names = NULL,
                                  sep = "")
spei.min.cv.delta.long[spei.min.cv.delta.long$time==1,]$time <- "1995"
spei.min.cv.delta.long[spei.min.cv.delta.long$time==2,]$time <- "2000"
spei.min.cv.delta.long[spei.min.cv.delta.long$time==3,]$time <- "2005"
spei.min.cv.delta.long[spei.min.cv.delta.long$time==4,]$time <- "2010"
spei.min.cv.delta.long[spei.min.cv.delta.long$time==5,]$time <- "2015"
spei.min.cv.delta.long[spei.min.cv.delta.long$time==6,]$time <- "2020"
# spei.mean.cv delta
spei.mean.cv.delta <- data.frame(Name = spei.mean.cv$COUNTRY,
                                 "1995" = rowMeans(spei.mean.cv[,3:7], na.rm = T),
                                 "2000" = rowMeans(spei.mean.cv[,8:12], na.rm = T), 
                                 "2005" = rowMeans(spei.mean.cv[,13:17], na.rm = T), 
                                 "2010" = rowMeans(spei.mean.cv[,18:22], na.rm = T), 
                                 "2015" = rowMeans(spei.mean.cv[,23:27], na.rm = T), 
                                 "2020" = rowMeans(spei.mean.cv[,28:32], na.rm = T))
spei.mean.cv.delta.long <- reshape(data = as.data.frame(spei.mean.cv.delta),
                                   varying = 2:7,
                                   v.names = "spei.mean.cv",
                                   idvar = c("Name"),
                                   direction = "long",
                                   new.row.names = NULL,
                                   sep = "")
spei.mean.cv.delta.long[spei.mean.cv.delta.long$time==1,]$time <- "1995"
spei.mean.cv.delta.long[spei.mean.cv.delta.long$time==2,]$time <- "2000"
spei.mean.cv.delta.long[spei.mean.cv.delta.long$time==3,]$time <- "2005"
spei.mean.cv.delta.long[spei.mean.cv.delta.long$time==4,]$time <- "2010"
spei.mean.cv.delta.long[spei.mean.cv.delta.long$time==5,]$time <- "2015"
spei.mean.cv.delta.long[spei.mean.cv.delta.long$time==6,]$time <- "2020"

# conflict.nConflicts delta
conflict.nConflicts.delta <- data.frame(Name = conflict.wide.ordered$country,
                                        "1995" = rowSums(as.matrix(conflict.wide.ordered[,c(34:38)]), na.rm = T),
                                        "2000" = rowSums(as.matrix(conflict.wide.ordered[,c(39:43)]), na.rm = T),
                                        "2005" = rowSums(as.matrix(conflict.wide.ordered[,c(44:48)]), na.rm = T),
                                        "2010" = rowSums(as.matrix(conflict.wide.ordered[,c(49:53)]), na.rm = T),
                                        "2015" = rowSums(as.matrix(conflict.wide.ordered[,c(54:58)]), na.rm = T),
                                        "2020" = rowSums(as.matrix(conflict.wide.ordered[,c(59:63)]), na.rm = T))
conflict.nConflicts.delta.long <- reshape(data = as.data.frame(conflict.nConflicts.delta),
                                          varying = 2:7,
                                          v.names = "nConflicts",
                                          idvar = c("Name"),
                                          direction = "long",
                                          new.row.names = NULL,
                                          sep = "")
conflict.nConflicts.delta.long[conflict.nConflicts.delta.long$time==1,]$time <- "1995"
conflict.nConflicts.delta.long[conflict.nConflicts.delta.long$time==2,]$time <- "2000"
conflict.nConflicts.delta.long[conflict.nConflicts.delta.long$time==3,]$time <- "2005"
conflict.nConflicts.delta.long[conflict.nConflicts.delta.long$time==4,]$time <- "2010"
conflict.nConflicts.delta.long[conflict.nConflicts.delta.long$time==5,]$time <- "2015"
conflict.nConflicts.delta.long[conflict.nConflicts.delta.long$time==6,]$time <- "2020"
# conflict.sum.dead delta
conflict.sum.dead.delta <- data.frame(Name = conflict.wide.ordered$country,
                                      "1995" = rowSums(as.matrix(conflict.wide.ordered[,c(3:7)]), na.rm = T),
                                      "2000" = rowSums(as.matrix(conflict.wide.ordered[,c(8:12)]), na.rm = T),
                                      "2005" = rowSums(as.matrix(conflict.wide.ordered[,c(13:17)]), na.rm = T),
                                      "2010" = rowSums(as.matrix(conflict.wide.ordered[,c(18:22)]), na.rm = T),
                                      "2015" = rowSums(as.matrix(conflict.wide.ordered[,c(23:27)]), na.rm = T),
                                      "2020" = rowSums(as.matrix(conflict.wide.ordered[,c(28:32)]), na.rm = T))
conflict.sum.dead.delta.long <- reshape(data = as.data.frame(conflict.sum.dead.delta),
                                        varying = 2:7,
                                        v.names = "conflict.sum.dead",
                                        idvar = c("Name"),
                                        direction = "long",
                                        new.row.names = NULL,
                                        sep = "")
conflict.sum.dead.delta.long[conflict.sum.dead.delta.long$time==1,]$time <- "1995"
conflict.sum.dead.delta.long[conflict.sum.dead.delta.long$time==2,]$time <- "2000"
conflict.sum.dead.delta.long[conflict.sum.dead.delta.long$time==3,]$time <- "2005"
conflict.sum.dead.delta.long[conflict.sum.dead.delta.long$time==4,]$time <- "2010"
conflict.sum.dead.delta.long[conflict.sum.dead.delta.long$time==5,]$time <- "2015"
conflict.sum.dead.delta.long[conflict.sum.dead.delta.long$time==6,]$time <- "2020"

# emdat.flood.affected.delta
emdat.flood.affected.delta <- data.frame(Name = emdat_flood_wide$Country,
                                         "1995" = rowSums(as.matrix(emdat_flood_wide[,c(3:7)]), na.rm = T),
                                         "2000" = rowSums(as.matrix(emdat_flood_wide[,c(8:12)]), na.rm = T),
                                         "2005" = rowSums(as.matrix(emdat_flood_wide[,c(13:17)]), na.rm = T),
                                         "2010" = rowSums(as.matrix(emdat_flood_wide[,c(18:22)]), na.rm = T),
                                         "2015" = rowSums(as.matrix(emdat_flood_wide[,c(23:27)]), na.rm = T),
                                         "2020" = rowSums(as.matrix(emdat_flood_wide[,c(28:31)]), na.rm = T))
emdat.flood.affected.delta.long <- reshape(data = as.data.frame(emdat.flood.affected.delta),
                                           varying = 2:7,
                                           v.names = "emdat.flood.affected",
                                           idvar = c("Name"),
                                           direction = "long",
                                           new.row.names = NULL,
                                           sep = "")
emdat.flood.affected.delta.long[emdat.flood.affected.delta.long$time==1,]$time <- "1995"
emdat.flood.affected.delta.long[emdat.flood.affected.delta.long$time==2,]$time <- "2000"
emdat.flood.affected.delta.long[emdat.flood.affected.delta.long$time==3,]$time <- "2005"
emdat.flood.affected.delta.long[emdat.flood.affected.delta.long$time==4,]$time <- "2010"
emdat.flood.affected.delta.long[emdat.flood.affected.delta.long$time==5,]$time <- "2015"
emdat.flood.affected.delta.long[emdat.flood.affected.delta.long$time==6,]$time <- "2020"

# emdat.flood.deaths.delta
emdat.flood.deaths.delta <- data.frame(Name = emdat_flood_wide$Country,
                                       "1995" = rowSums(as.matrix(emdat_flood_wide[,c(34:38)]), na.rm = T),
                                       "2000" = rowSums(as.matrix(emdat_flood_wide[,c(39:43)]), na.rm = T),
                                       "2005" = rowSums(as.matrix(emdat_flood_wide[,c(44:48)]), na.rm = T),
                                       "2010" = rowSums(as.matrix(emdat_flood_wide[,c(49:53)]), na.rm = T),
                                       "2015" = rowSums(as.matrix(emdat_flood_wide[,c(54:58)]), na.rm = T),
                                       "2020" = rowSums(as.matrix(emdat_flood_wide[,c(59:63)]), na.rm = T))
emdat.flood.deaths.delta.long <- reshape(data = as.data.frame(emdat.flood.deaths.delta),
                                         varying = 2:7,
                                         v.names = "emdat.flood.deaths",
                                         idvar = c("Name"),
                                         direction = "long",
                                         new.row.names = NULL,
                                         sep = "")
emdat.flood.deaths.delta.long[emdat.flood.deaths.delta.long$time==1,]$time <- "1995"
emdat.flood.deaths.delta.long[emdat.flood.deaths.delta.long$time==2,]$time <- "2000"
emdat.flood.deaths.delta.long[emdat.flood.deaths.delta.long$time==3,]$time <- "2005"
emdat.flood.deaths.delta.long[emdat.flood.deaths.delta.long$time==4,]$time <- "2010"
emdat.flood.deaths.delta.long[emdat.flood.deaths.delta.long$time==5,]$time <- "2015"
emdat.flood.deaths.delta.long[emdat.flood.deaths.delta.long$time==6,]$time <- "2020"

# emdat.storm.affected.delta
emdat.storm.affected.delta <- data.frame(Name = emdat_storm_wide$Country,
                                         "1995" = rowSums(as.matrix(emdat_storm_wide[,c(3:7)]), na.rm = T),
                                         "2000" = rowSums(as.matrix(emdat_storm_wide[,c(8:12)]), na.rm = T),
                                         "2005" = rowSums(as.matrix(emdat_storm_wide[,c(13:17)]), na.rm = T),
                                         "2010" = rowSums(as.matrix(emdat_storm_wide[,c(18:22)]), na.rm = T),
                                         "2015" = rowSums(as.matrix(emdat_storm_wide[,c(23:27)]), na.rm = T),
                                         "2020" = rowSums(as.matrix(emdat_storm_wide[,c(28:31)]), na.rm = T))
emdat.storm.affected.delta.long <- reshape(data = as.data.frame(emdat.storm.affected.delta),
                                           varying = 2:7,
                                           v.names = "emdat.storm.affected",
                                           idvar = c("Name"),
                                           direction = "long",
                                           new.row.names = NULL,
                                           sep = "")
emdat.storm.affected.delta.long[emdat.storm.affected.delta.long$time==1,]$time <- "1995"
emdat.storm.affected.delta.long[emdat.storm.affected.delta.long$time==2,]$time <- "2000"
emdat.storm.affected.delta.long[emdat.storm.affected.delta.long$time==3,]$time <- "2005"
emdat.storm.affected.delta.long[emdat.storm.affected.delta.long$time==4,]$time <- "2010"
emdat.storm.affected.delta.long[emdat.storm.affected.delta.long$time==5,]$time <- "2015"
emdat.storm.affected.delta.long[emdat.storm.affected.delta.long$time==6,]$time <- "2020"

# emdat.storm.deaths.delta
emdat.storm.deaths.delta <- data.frame(Name = emdat_storm_wide$Country,
                                       "1995" = rowSums(as.matrix(emdat_storm_wide[,c(34:38)]), na.rm = T),
                                       "2000" = rowSums(as.matrix(emdat_storm_wide[,c(39:43)]), na.rm = T),
                                       "2005" = rowSums(as.matrix(emdat_storm_wide[,c(44:48)]), na.rm = T),
                                       "2010" = rowSums(as.matrix(emdat_storm_wide[,c(49:53)]), na.rm = T),
                                       "2015" = rowSums(as.matrix(emdat_storm_wide[,c(54:58)]), na.rm = T),
                                       "2020" = rowSums(as.matrix(emdat_storm_wide[,c(59:63)]), na.rm = T))
emdat.storm.deaths.delta.long <- reshape(data = as.data.frame(emdat.storm.deaths.delta),
                                         varying = 2:7,
                                         v.names = "emdat.storm.deaths",
                                         idvar = c("Name"),
                                         direction = "long",
                                         new.row.names = NULL,
                                         sep = "")
emdat.storm.deaths.delta.long[emdat.storm.deaths.delta.long$time==1,]$time <- "1995"
emdat.storm.deaths.delta.long[emdat.storm.deaths.delta.long$time==2,]$time <- "2000"
emdat.storm.deaths.delta.long[emdat.storm.deaths.delta.long$time==3,]$time <- "2005"
emdat.storm.deaths.delta.long[emdat.storm.deaths.delta.long$time==4,]$time <- "2010"
emdat.storm.deaths.delta.long[emdat.storm.deaths.delta.long$time==5,]$time <- "2015"
emdat.storm.deaths.delta.long[emdat.storm.deaths.delta.long$time==6,]$time <- "2020"


## Merging additional environmental and conflict data
mig_stacked_score <- merge(pop.delta.long, gdp.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, spei.max.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, spei.min.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, spei.mean.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, spei.max.cv.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, spei.min.cv.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, spei.mean.cv.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, conflict.nConflicts.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, conflict.sum.dead.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, emdat.flood.affected.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, emdat.flood.deaths.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, emdat.storm.affected.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score <- merge(mig_stacked_score, emdat.storm.deaths.delta.long, by.x = c("Name", "time"), by.y = c("Name","time"), all = TRUE)
mig_stacked_score[mig_stacked_score=="-Inf"] <- NA

mig_stacked_score <- mig_stacked_score[mig_stacked_score$time>=1990 & mig_stacked_score$time<=2020,]
mig_stacked <- mig_stacked_score[mig_stacked_score$Name %in% States$name,]

## Averaging variables by POP
mig_stacked$conflict.sum.dead.percap <- mig_stacked$conflict.sum.dead/mig_stacked$pop
mig_stacked$emdat.flood.affected.percap <- mig_stacked$emdat.flood.affected/mig_stacked$pop
mig_stacked$emdat.flood.deaths.percap <- mig_stacked$emdat.flood.deaths/mig_stacked$pop
mig_stacked$emdat.storm.affected.percap <- mig_stacked$emdat.storm.affected/mig_stacked$pop
mig_stacked$emdat.storm.deaths.percap <- mig_stacked$emdat.storm.deaths/mig_stacked$pop

mig_stacked <- mig_stacked[,!colnames(mig_stacked) %in% c("conflict.sum.dead","emdat.flood.affected","emdat.flood.deaths","emdat.storm.affected","emdat.storm.deaths")]


## Creating variable dataset with no lag. 
mig_stacked_nolag <- merge(undesa_roles, mig_stacked[,!colnames(mig_stacked) %in% "pop"], by.x = c("name", "time"), by.y = c("Name","time"), all = TRUE)

# mig_stacked_nolag <- mig_stacked_nolag[,c(1:12,14,13,16,15)]
colnames(mig_stacked_nolag) <- c("NAME","YEAR","role","gdp","spei.max","spei.min","spei.mean","spei.max.cv","spei.min.cv","spei.mean.cv","nConflicts","conflict.sum.dead","emdat.flood.affected","emdat.flood.deaths","emdat.storm.affected","emdat.storm.deaths")

mig_stacked <- mig_stacked_nolag[mig_stacked_nolag$NAME %in% States$name, colnames(mig_stacked_nolag)%in% c("NAME",
                                                                                                            "YEAR",
                                                                                                            "role",
                                                                                                            "gdp",
                                                                                                            "spei.max",
                                                                                                            "spei.min",
                                                                                                            "spei.mean",
                                                                                                            "spei.max.cv",
                                                                                                            "spei.min.cv",
                                                                                                            "spei.mean.cv",
                                                                                                            "nConflicts",
                                                                                                            "conflict.sum.dead"
                                                                                                            ,
                                                                                                            "emdat.flood.affected",
                                                                                                            "emdat.flood.deaths"
                                                                                                            ,
                                                                                                            "emdat.storm.affected",
                                                                                                            "emdat.storm.deaths"
)]
mig_stacked$role <- as.numeric(mig_stacked$role)

### Scaling Attributes [NOT SCALING FOR MBLOGIT MODELS]
# mig_stacked[,4:ncol(mig_stacked)] <- scale(mig_stacked[,4:ncol(mig_stacked)], )

### Create binary role variable for desired role category:
node_attr <- mig_stacked
node_attr$role <- ifelse(node_attr$role=="1", # set desired role for model here
                         1,
                         0)

networks <- vector("list", length = length(seq(from=1995, to=2020, 5)))
l.role <- 
  l.gdp <- 
  l.spei.max <- 
  l.spei.min <- 
  l.spei.mean <- 
  l.spei.max.cv <- 
  l.spei.min.cv <- 
  l.spei.mean.cv <- 
  l.nConflicts <- 
  l.conflict.sum.dead <- 
  l.emdat.flood.affected <-
  l.emdat.flood.deaths <- 
  l.emdat.storm.affected <- 
  l.emdat.storm.deaths <- list()

for (i in seq(from=1995, to=2020, 5)){
  j <- (i/5)-398
  
  l.role[[j]] <- node_attr[node_attr$YEAR==i, c("role")]
  l.gdp[[j]] <- node_attr[node_attr$YEAR==i, c("gdp")]
  l.spei.max[[j]] <- node_attr[node_attr$YEAR==i, c("spei.max")]
  l.spei.min[[j]] <- node_attr[node_attr$YEAR==i, c("spei.min")]
  l.spei.mean[[j]] <- node_attr[node_attr$YEAR==i, c("spei.mean")]
  l.spei.max.cv[[j]] <- node_attr[node_attr$YEAR==i, c("spei.max.cv")]
  l.spei.min.cv[[j]] <- node_attr[node_attr$YEAR==i, c("spei.min.cv")]
  l.spei.mean.cv[[j]] <- node_attr[node_attr$YEAR==i, c("spei.mean.cv")]
  l.nConflicts[[j]] <- node_attr[node_attr$YEAR==i, c("nConflicts")]
  l.conflict.sum.dead[[j]] <- node_attr[node_attr$YEAR==i, c("conflict.sum.dead")]
  l.emdat.flood.affected[[j]] <- node_attr[node_attr$YEAR==i, c("emdat.flood.affected")]
  l.emdat.flood.deaths[[j]] <- node_attr[node_attr$YEAR==i, c("emdat.flood.deaths")]
  l.emdat.storm.affected[[j]] <- node_attr[node_attr$YEAR==i, c("emdat.storm.affected")]
  l.emdat.storm.deaths[[j]] <- node_attr[node_attr$YEAR==i, c("emdat.storm.deaths")]
  
  names(l.role[[j]]) <- 
    names(l.gdp[[j]]) <-
    names(l.spei.max[[j]]) <-
    names(l.spei.min[[j]]) <-
    names(l.spei.mean[[j]]) <-
    names(l.spei.max.cv[[j]]) <-
    names(l.spei.min.cv[[j]]) <-
    names(l.spei.mean.cv[[j]]) <-
    names(l.nConflicts[[j]]) <-
    names(l.conflict.sum.dead[[j]]) <-
    names(l.emdat.flood.affected[[j]]) <-
    names(l.emdat.flood.deaths[[j]]) <-
    names(l.emdat.storm.affected[[j]]) <-
    names(l.emdat.storm.deaths[[j]]) <-
    node_attr[node_attr$YEAR==i, c("NAME")]
  
  network_attributes <- edges.mig.wide[is.na(edges.mig.wide$state.origin.name)==F & is.na(edges.mig.wide$state.destination.name)==F,
                                       c(colnames(edges.mig.wide)[c(1,2)], colnames(edges.mig.wide)[grepl(i, colnames(edges.mig.wide))])]
  colnames(network_attributes) <- c("from","to","flow")
  network_attributes[,c("flow")] <- as.numeric(network_attributes[,c("flow")])
  network_attributes <- network_attributes[order(network_attributes[,c("from","to")]),]
  network_attributes[is.na(network_attributes)] <- 0
  networks[[j]] <- as.matrix(as_adjacency_matrix(graph_from_data_frame(network_attributes), attr = "flow", names = T))
}

names(l.role) <-
  names(l.gdp) <-
  names(l.spei.max) <-
  names(l.spei.min) <-
  names(l.spei.mean) <-
  names(l.spei.max.cv) <-
  names(l.spei.min.cv) <-
  names(l.spei.mean.cv) <-
  names(l.nConflicts) <-
  names(l.conflict.sum.dead) <-
  names(l.emdat.flood.affected) <-
  names(l.emdat.flood.deaths) <-
  names(l.emdat.storm.affected) <-
  names(l.emdat.storm.deaths) <-
  names(networks) <- paste0("t",1:length(networks))

### Model 1: Role = 1 #note{SIGNIFICANT: conflict.sum.dead, storm.deaths (unscaled) | gdp, spei.min.cv, flood.affected (scaled)}
model1 <- tnam(
  l.role ~
  covariate(l.gdp, coefname = "gdp")
  +
  # covariate(l.spei.max, coefname = "spei.max")
  # +
  # covariate(l.spei.min, coefname = "spei.min")
  # +
  # covariate(l.spei.mean, coefname = "spei.mean")
  # +
  #   covariate(l.spei.max.cv, coefname = "spei.max.cv")
  # +
  # covariate(l.spei.min.cv, coefname = "spei.min.cv")
  # +
  # covariate(l.spei.mean.cv, coefname = "spei.mean.cv")
  # +
  # covariate(l.nConflicts, coefname = "nConflicts")
  # +
  # covariate(l.conflict.sum.dead, coefname = "conflict.sum.dead", exponent = 2)
  # +
  covariate(l.emdat.flood.affected, coefname = "flood.affected")
  # +
  # covariate(l.emdat.flood.deaths, coefname = "flood.deaths")
  # +
  # covariate(l.emdat.storm.affected, coefname = "storm.affected")
  # +
    # covariate(l.emdat.storm.deaths, coefname = "storm.deaths")
  ,
  re.node = TRUE, time.linear = TRUE, family = binomial, na.action = na.omit
)
summary(model1)
screenreg(model1)

### Model 2: Role = 2 #note{SIGNIFICANT: flood.affected (unscaled) | gdp, spei.min.cv, flood.affected (scaled)}
model2 <- tnam(
  l.role ~
    covariate(l.gdp, coefname = "gdp")
    +
      # covariate(l.spei.max, coefname = "spei.max")
    # +
      # covariate(l.spei.min, coefname = "spei.min")
    # +
    # covariate(l.spei.mean, coefname = "spei.mean")
  # +
    # covariate(l.spei.max.cv, coefname = "spei.max.cv")
  # +
  #   covariate(l.spei.min.cv, coefname = "spei.min.cv")
  # +
  # covariate(l.spei.mean.cv, coefname = "spei.mean.cv")
  # +
    # covariate(l.nConflicts, coefname = "nConflicts")
  # +
  # covariate(l.conflict.sum.dead, coefname = "conflict.sum.dead")
  # +
  covariate(l.emdat.flood.affected, coefname = "flood.affected")
  # +
  # covariate(l.emdat.flood.deaths, coefname = "flood.deaths")
  # +
  # covariate(l.emdat.storm.affected, coefname = "storm.affected")
  # +
  # covariate(l.emdat.storm.deaths, coefname = "storm.deaths")
  ,
  re.node = TRUE, time.linear = TRUE, family = binomial, na.action = na.omit
)
summary(model2)
screenreg(model2)

### Model 3: Role = 3 #note{SIGNIFICANT: NULL (unscaled) | NULL (scaled)}
model3 <- tnam(
  l.role ~
  # covariate(l.gdp, coefname = "gdp")
  # +
  # covariate(l.spei.max, coefname = "spei.max")
  # +
  # covariate(l.spei.min, coefname = "spei.min")
  # +
  # covariate(l.spei.mean, coefname = "spei.mean")
  # +
  # covariate(l.spei.max.cv, coefname = "spei.max.cv")
  # +
  # covariate(l.spei.min.cv, coefname = "spei.min.cv")
  # +
  # covariate(l.spei.mean.cv, coefname = "spei.mean.cv")
  # +
  # covariate(l.nConflicts, coefname = "nConflicts")
  # +
  # covariate(l.conflict.sum.dead, coefname = "conflict.sum.dead")
  # +
  # covariate(l.emdat.flood.affected, coefname = "flood.affected")
  # +
  # covariate(l.emdat.flood.deaths, coefname = "flood.deaths")
  # +
  # covariate(l.emdat.storm.affected, coefname = "storm.affected")
  # +
  # covariate(l.emdat.storm.deaths, coefname = "storm.deaths")
  ,
  re.node = TRUE, time.linear = TRUE, family = binomial, na.action = na.omit
)
summary(model3)
screenreg(model3)

### K = 3 Exogenous Model Outputs: Odds Ratios and Likelihood Percentages (of scaled covariates and with SPEI.MIN.CV) [THESE ARE JUST FOR COMPARISON]
## Model 1
exp(summary(model1)$coefficients[-c(1,2),1])
(exp(summary(model1)$coefficients[-c(1,2),1])-1)*100
# A 1 standard deviation increase in the 5-year mean of per-capita GDP is associated with a 99.99998% (eβ = 0.00000...) decrease in the odds of belonging to Role 1.
# A 1 standard deviation increase in the coefficient of variation of 5-year minimum SPEI values is associated with a 100.00000% (eβ = 0.00000...) decrease in the odds of belonging to Role 1.
# A 1 standard deviation increase in the 5-year total of per-capita Flood Affected is associated with a 737.23825% (eβ = 8.37238...) increase in the odds of belonging to Role 1.

# Model 2
exp(summary(model2)$coefficients[-c(1,2),1])
(exp(summary(model2)$coefficients[-c(1,2),1])-1)*100
# A 1 standard deviation increase in the 5-year mean of per-capita GDP is associated with a 278669662979.71472% (eβ = 2786696630.79714727) increase in the odds of belonging to Role 1.
# A 1 standard deviation increase in the coefficient of variation of 5-year minimum SPEI values is associated with a 3599292693109087922068646868262880244622826080420820028044444226882820286222826406004242428468442646420648266080222620806464442420242864428264628608628260648084664800268686262084860620868888486200224288246844080086482204080008066622626628264062246602844.00000% (eβ = 35992926931090877700224448664428028040646866066224444628664428288804482206062684288048888062066642040440840686084082686668086044668608620026468068882684026482866420444266206800088082282268004046426204626602220468448220824088642828264882402626642044004.00000000) increase in the odds of belonging to Role 1.
# A 1 standard deviation increase in the 5-year total of per-capita Flood Affected is associated with a 96.65832 % (eβ = 0.03341681) decrease in the odds of belonging to Role 1.

# Model 3
exp(summary(model3)$coefficients[-c(1,2),1])
(exp(summary(model3)$coefficients[-c(1,2),1])-1)*100
# NULL MODEL

### K = 3 Model Outputs: Odds Ratios and Likelihood Percentages (of scaled covariates and without SPEI.MIN.CV)
## Model 1
exp(summary(model1)$coefficients[-c(1,2),1])
(exp(summary(model1)$coefficients[-c(1,2),1])-1)*100
# A 1 standard deviation increase in the 5-year mean of per-capita GDP is associated with a 99.9951% (eβ = 0.00004899923) decrease in the odds of belonging to Role 1.
# A 1 standard deviation increase in the 5-year total of per-capita Flood Affected is associated with a 479.2460% (eβ = 5.79245996343) increase in the odds of belonging to Role 1.

# Model 2
exp(summary(model2)$coefficients[-c(1,2),1])
(exp(summary(model2)$coefficients[-c(1,2),1])-1)*100
# A 1 standard deviation increase in the 5-year mean of per-capita GDP is associated with a 1490028.08861% (eβ = 14901.28088612) increase in the odds of belonging to Role 1.
# A 1 standard deviation increase in the 5-year total of per-capita Flood Affected is associated with a 90.78986% (eβ = 0.09210139) decrease in the odds of belonging to Role 1.

# Model 3
exp(summary(model3)$coefficients[-c(1,2),1])
(exp(summary(model3)$coefficients[-c(1,2),1])-1)*100
# NULL MODEL



### Migrant Flow Network Plots
for (i in 1:ncol(edges.mig.wide[,grep("migflows",colnames(edges.mig.wide))])){
  a <- na.omit(edges.mig.wide[edges.mig.wide[,c(paste("migflows", i*5+1990 , sep = "."))]>=1,][,c("state.origin.name", "state.destination.name", paste("migflows", i*5+1990, sep = "."))])
  colnames(a)[3] <- "flow"  
  b <- igraph::graph_from_data_frame(a, directed = T)
  c <- mig_stacked[mig_stacked$NAME%in%V(b)$name & mig_stacked$YEAR==i*5+1990,c(1,3)]
  set.seed(123)
  plot(b, 
       layout = layout.kamada.kawai(b),
       vertex.size = log(strength(b, mode = "total", weights = E(b)$flow))/2,
       vertex.color = custom_colors[c[match(V(b)$name,c$NAME),]$role],
       vertex.frame.color = custom_colors[c[match(V(b)$name,c$NAME),]$role],
       vertex.label = NA,
       edge.arrow.size = 0.1,
       edge.width = log(E(b)$flow)/2,
       edge.color = rgb(0,0,0,((range01(E(b)$flow))*0.19)+0.01),
       # edge.label = E(b)$flow,
       # edge.label.color = "#d53e4f",
       # edge.label.cex = 0.7,
       sub = paste0(i*5+1990, sep=""))
}



## mblogit Models
# Model 1
node_attr <- mig_stacked
node_attr$role <- ifelse(node_attr$role=="1", # set desired role for model here
                         1,
                         0)
node_attr$role <- as.factor(node_attr$role)
node_attr$role <- relevel(node_attr$role, ref = "0")

mblogit_model1 <- mblogit(formula = role ~ 
                            gdp +
                            # spei.max +
                            # spei.min +
                            # spei.mean +
                            # spei.max.cv +
                            # spei.min.cv +
                            # spei.mean.cv +
                            # nConflicts +
                            # conflict.sum.dead +
                          #   emdat.flood.affected
                          # +
                            emdat.flood.deaths
                          ,
                          random = ~1|NAME,
                          data = node_attr,
                          estimator = "ML",
                          method = "PQL",
                          maxit = 100
)
summary(mblogit_model1)

# Model 2
node_attr <- mig_stacked
node_attr$role <- ifelse(node_attr$role=="2", # set desired role for model here
                         1,
                         0)
node_attr$role <- as.factor(node_attr$role)
node_attr$role <- relevel(node_attr$role, ref = "0")


mblogit_model2 <- mblogit(formula = role ~ 
                            gdp +
                            # spei.max +
                          #   spei.min
                          # +
                          # spei.mean +
                          # spei.max.cv +
                          # spei.min.cv
                          # +
                          # spei.mean.cv +
                          # nConflicts +
                          # conflict.sum.dead
                          # +
                          #   emdat.flood.affected
                          # +
                          emdat.flood.deaths
                          ,
                          random = ~1|NAME,
                          data = node_attr,
                          estimator = "ML",
                          method = "PQL",
                          maxit = 100
)
summary(mblogit_model2)

# Model 3
node_attr <- mig_stacked
node_attr$role <- ifelse(node_attr$role=="3", # set desired role for model here
                         1,
                         0)
node_attr$role <- as.factor(node_attr$role)
node_attr$role <- relevel(node_attr$role, ref = "0")

mblogit_model3 <- mblogit(formula = role ~ 
                            gdp +
                            # spei.max +
                              spei.min
                            +
                            # spei.mean +
                            # spei.max.cv +
                            # spei.min.cv
                            # +
                            # spei.mean.cv +
                          #   nConflicts
                          # +
                          # conflict.sum.dead
                          # +
                          #   emdat.flood.affected
                          # +
                          emdat.flood.deaths
                          ,
                          random = ~1|NAME,
                          data = node_attr,
                          estimator = "ML",
                          method = "PQL",
                          maxit = 100
)
summary(mblogit_model3)

### Mblogit Model Outputs
sink(file = "results/mblogit_UNDESA_roles_model_results.txt")
print("Model: 1")
summary(mblogit_model1)
print("Model: 2")
summary(mblogit_model2)
print("Model: 3")
summary(mblogit_model3)
sink(file = NULL)
