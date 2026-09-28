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
unhcr_roles <- read.csv("results/UNHCR_full_roles_k3.csv", header = T)
unhcr_roles$X <- sub("\\.19.*", "", unhcr_roles$X)
unhcr_roles$X <- sub("\\.20.*", "", unhcr_roles$X)
unhcr_roles$id <- sub(".*\\.", "", unhcr_roles$id)
colnames(unhcr_roles) <- c("name","time","role")

### country reference from bilat_mig
States <- read.csv("data/country_list.csv", header = T)

### UNHCR Refugee Flows 1990-2020
edges.ref <- read.csv("./data/UNHCR_FLOW_Data_REF.csv", header = T)
edges.ref$state.origin.name <- States[match(edges.ref$orig.iso3, States$iso3c),]$name
edges.ref$state.destination.name <- States[match(edges.ref$dest.iso3, States$iso3c),]$name
edges.ref <- edges.ref[,-c(2,4)]

# UN POP data
pop <- read.csv("./data/WPP2022_Demographic_Indicators_Medium.csv") # call in UN world pop and demographic data
pop <- pop[as.numeric(pop$Time)>=1990 & as.numeric(pop$Time)<=2020,] # retain data between 1990 and 2020
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

spei.max.long <- reshape(data = as.data.frame(spei.max),
                         varying = 2:32,
                         v.names = "spei.max",
                         idvar = c("COUNTRY"),
                         direction = "long",
                         new.row.names = NULL,
                         sep = "")
spei.max.long$time <- spei.max.long$time + 1989

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

spei.min.long <- reshape(data = as.data.frame(spei.min),
                         varying = 2:32,
                         v.names = "spei.min",
                         idvar = c("COUNTRY"),
                         direction = "long",
                         new.row.names = NULL,
                         sep = "")
spei.min.long$time <- spei.min.long$time + 1989

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

spei.mean.long <- reshape(data = as.data.frame(spei.mean),
                          varying = 2:32,
                          v.names = "spei.mean",
                          idvar = c("COUNTRY"),
                          direction = "long",
                          new.row.names = NULL,
                          sep = "")
spei.mean.long$time <- spei.mean.long$time + 1989

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

spei.max.cv.long <- reshape(data = as.data.frame(spei.max.cv),
                            varying = 2:32,
                            v.names = "spei.max.cv",
                            idvar = c("COUNTRY"),
                            direction = "long",
                            new.row.names = NULL,
                            sep = "")
spei.max.cv.long$time <- spei.max.cv.long$time + 1989

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

spei.min.cv.long <- reshape(data = as.data.frame(spei.min.cv),
                            varying = 2:32,
                            v.names = "spei.min.cv",
                            idvar = c("COUNTRY"),
                            direction = "long",
                            new.row.names = NULL,
                            sep = "")
spei.min.cv.long$time <- spei.min.cv.long$time + 1989

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

spei.mean.cv.long <- reshape(data = as.data.frame(spei.mean.cv),
                             varying = 2:32,
                             v.names = "spei.mean.cv",
                             idvar = c("COUNTRY"),
                             direction = "long",
                             new.row.names = NULL,
                             sep = "")
spei.mean.cv.long$time <- spei.mean.cv.long$time + 1989

### Call in Uppsala Conflict Data
conflict <- read.csv("data/GEDEvent_v22_1.csv")
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

# calculate per capita 
conflict.dat <- merge(conflict.dat, pop, by.x = c("country","year"), by.y = c("name","year0"), all.x = TRUE)
conflict.dat$nConflicts <- conflict.dat$nConflicts/conflict.dat$pop
conflict.dat$conflict.sum.dead <- conflict.dat$conflict.sum.dead/conflict.dat$pop
conflict.dat <- conflict.dat[,!colnames(conflict.dat)%in%c("pop")]

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

# EMDAT Storm
emdat_storm <- disaster[disaster$Disaster.Type=="Storm",c(1,3:5)]
colnames(emdat_storm)[3:4] <- paste("storm",colnames(emdat_storm)[3:4],sep = ".")
emdat_storm <- as.data.frame(emdat_storm) %>% 
  group_by(Country, Year) %>% 
  summarise(storm.Total.Deaths = sum(storm.Total.Deaths, na.rm = T),
            storm.Total.Affected = sum(storm.Total.Affected, na.rm = T))

# calulate per capita
emdat_storm <- merge(emdat_storm, pop, by.x = c("Country","Year"), by.y = c("name","year0"), all.x = TRUE)
emdat_storm$storm.Total.Deaths <- emdat_storm$storm.Total.Deaths/emdat_storm$pop
emdat_storm$storm.Total.Affected <- emdat_storm$storm.Total.Affected/emdat_storm$pop
emdat_storm <- emdat_storm[,!colnames(emdat_storm)%in%c("pop")]

# EMDAT Flood
emdat_flood <- disaster[disaster$Disaster.Type=="Flood",c(1,3:5)]
colnames(emdat_flood)[3:4] <- paste("flood",colnames(emdat_flood)[3:4],sep = ".")
emdat_flood <- as.data.frame(emdat_flood) %>% 
  group_by(Country, Year) %>% 
  summarise(flood.Total.Deaths = sum(flood.Total.Deaths, na.rm = T),
            flood.Total.Affected = sum(flood.Total.Affected, na.rm = T))

# calulate per capita
emdat_flood <- merge(emdat_flood, pop, by.x = c("Country","Year"), by.y = c("name","year0"), all.x = TRUE)
emdat_flood$flood.Total.Deaths <- emdat_flood$flood.Total.Deaths/emdat_flood$pop
emdat_flood$flood.Total.Affected <- emdat_flood$flood.Total.Affected/emdat_flood$pop
emdat_flood <- emdat_flood[,!colnames(emdat_flood)%in%c("pop")]

## Creating variable dataset with no lag. 
ref_stacked_nolag <- merge(unhcr_roles, gdp.long, by.x = c("name", "time"), by.y = c("name","time"), all = TRUE)
ref_stacked_nolag <- merge(ref_stacked_nolag, spei.max.long, by.x = c("name", "time"), by.y = c("COUNTRY","time"), all = TRUE)
ref_stacked_nolag <- merge(ref_stacked_nolag, spei.min.long, by.x = c("name", "time"), by.y = c("COUNTRY","time"), all = TRUE)
ref_stacked_nolag <- merge(ref_stacked_nolag, spei.mean.long, by.x = c("name", "time"), by.y = c("COUNTRY","time"), all = TRUE)
ref_stacked_nolag <- merge(ref_stacked_nolag, spei.max.cv.long, by.x = c("name", "time"), by.y = c("COUNTRY","time"), all = TRUE)
ref_stacked_nolag <- merge(ref_stacked_nolag, spei.min.cv.long, by.x = c("name", "time"), by.y = c("COUNTRY","time"), all = TRUE)
ref_stacked_nolag <- merge(ref_stacked_nolag, spei.mean.cv.long, by.x = c("name", "time"), by.y = c("COUNTRY","time"), all = TRUE)
ref_stacked_nolag <- merge(ref_stacked_nolag, conflict.dat, by.x = c("name", "time"), by.y = c("country","year"), all = TRUE)
ref_stacked_nolag <- merge(ref_stacked_nolag, emdat_flood, by.x = c("name", "time"), by.y = c("Country","Year"), all = TRUE)
ref_stacked_nolag <- merge(ref_stacked_nolag, emdat_storm, by.x = c("name", "time"), by.y = c("Country","Year"), all = TRUE)


ref_stacked_nolag <- ref_stacked_nolag[,c(1:12,14,13,16,15)]
colnames(ref_stacked_nolag) <- c("NAME","YEAR","role","gdp","spei.max","spei.min","spei.mean","spei.max.cv","spei.min.cv","spei.mean.cv","nConflicts","conflict.sum.dead","emdat.flood.affected","emdat.flood.deaths","emdat.storm.affected","emdat.storm.deaths")

ref_stacked <- ref_stacked_nolag[ref_stacked_nolag$NAME %in% States$name, colnames(ref_stacked_nolag)%in% c("NAME",
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
ref_stacked$role <- as.numeric(ref_stacked$role)

### Scaling Attributes [NOT SCALING FOR MBLOGIT MODELS]
# ref_stacked[,4:ncol(ref_stacked)] <- scale(ref_stacked[,4:ncol(ref_stacked)])

### Create binary role variable for desired role category:
node_attr <- ref_stacked
node_attr$role <- ifelse(node_attr$role=="1", # set desired role for model here
1,
0)

networks <- vector("list", length = length(1990:2020))
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

for (i in 1990:2020){
  j <- i-1989
  
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
  # l.emdat.storm.affected[[j]] <- node_attr[node_attr$YEAR==i, c("emdat.storm.affected")]
  # l.emdat.storm.deaths[[j]] <- node_attr[node_attr$YEAR==i, c("emdat.storm.deaths")]
  
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
    # names(l.emdat.storm.affected[[j]]) <-
    # names(l.emdat.storm.deaths[[j]]) <-
    node_attr[node_attr$YEAR==i, c("NAME")]
  
  network_attributes <- edges.ref[is.na(edges.ref$state.origin.name)==F & is.na(edges.ref$state.destination.name)==F, c(colnames(edges.ref)[c(1,2)], colnames(edges.ref)[grepl(i, colnames(edges.ref))])]
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
  # names(l.emdat.storm.affected) <-
  # names(l.emdat.storm.deaths) <-
  names(networks) <- paste0("t",1:31)

### Model 1: Role = 1 #note{SIGNIFICANT: gdp, flood.affected}
model1 <- tnam(
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
  covariate(l.conflict.sum.dead, coefname = "conflict.sum.dead")
  # +
  # covariate(l.emdat.flood.affected, coefname = "flood.affected")
  # +
  # covariate(l.emdat.flood.deaths, coefname = "flood.deaths")
  # +
  # covariate(l.emdat.storm.affected, coefname = "storm.affected")
  # +
  # covariate(l.emdat.storm.deaths, coefname = "storm.deaths")
  ,
  re.node = TRUE,
  time.linear = TRUE,
  family = binomial,
  na.action = na.omit
)
summary(model1)
screenreg(model1)

### Model 2: Role = 2 #note{SIGNIFICANT: gdp, spei.mean}
model2 <- tnam(
  l.role ~
  covariate(l.gdp, coefname = "gdp")
  +
  # covariate(l.spei.max, coefname = "spei.max")
  # +
  # covariate(l.spei.min, coefname = "spei.min")
  # +
  covariate(l.spei.mean, coefname = "spei.mean")
  # +
  #   covariate(l.spei.max.cv, coefname = "spei.max.cv")
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
summary(model2)
screenreg(model2)

### Model 3: Role = 3 #note{SIGNIFICANT: nConflicts}
model3 <- tnam(
  l.role ~
  #   covariate(l.gdp, coefname = "gdp")
  # +
  #   covariate(l.spei.max, coefname = "spei.max")
  # +
  #   covariate(l.spei.min, coefname = "spei.min")
  # +
  #   covariate(l.spei.mean, coefname = "spei.mean")
  # +
  #   covariate(l.spei.max.cv, coefname = "spei.max.cv")
  # +
  #   covariate(l.spei.min.cv, coefname = "spei.min.cv")
  # +
  #   covariate(l.spei.mean.cv, coefname = "spei.mean.cv")
  # +
    covariate(l.nConflicts, coefname = "nConflicts")
  # +
  #   covariate(l.conflict.sum.dead, coefname = "conflict.sum.dead")
  # +
  #   covariate(l.emdat.flood.affected, coefname = "flood.affected")
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


### K = 3 Exogeneous Model Outputs: Odds Ratios and Likelihood Percentages (of scaled covariates) [THESE ARE JUST FOR COMPARISON]
## Model 1
exp(summary(model1)$coefficients[-c(1,2),1])
(exp(summary(model1)$coefficients[-c(1,2),1])-1)*100
# NULL MODEL IGNORE BELOW
# A 1 standard deviation increase in the annual per-capita GDP is associated with a 601.58240% (eβ = 7.0158240) increase in the odds of belonging to Role 1.
# A 1 standard deviation increase in the annual per-capita Flood Affected is associated with a 37.44448% (eβ = 0.6255552) decrease in the odds of belonging to Role 1.

# Model 2
exp(summary(model2)$coefficients[-c(1,2),1])
(exp(summary(model2)$coefficients[-c(1,2),1])-1)*100
# A 1 standard deviation increase in per-capita GDP is associated with a 98.37475% (eβ = 0.01625247) decrease in the odds of belonging to Role 2.
# A 1 standard deviation increase in mean SPEI values is associated with a 19.83296% (eβ = 0.80167036) decrease in the odds of belonging to Role 2.

# Model 3
exp(summary(model3)$coefficients[-c(1,2),1])
(exp(summary(model3)$coefficients[-c(1,2),1])-1)*100
# A 1 standard deviation increase in the number of annual conflicts per-capita is associated with a 99.2634% (eβ = 0.007366006) decrease in the odds of belonging to Role 3.


## Refugee Flow Networks
for (i in 1:ncol(edges.ref[,grep("refflows",colnames(edges.ref))])){
  a <- na.omit(edges.ref[edges.ref[,c(paste("refflows", i+1989 , sep = "."))]>=1,][,c("state.origin.name", "state.destination.name", paste("refflows", i+1989, sep = "."))])
  colnames(a)[3] <- "flow"  
  b <- igraph::graph_from_data_frame(a, directed = T)
  c <- ref_stacked[ref_stacked$NAME%in%V(b)$name & ref_stacked$YEAR==i+1989,c(1,3)]
  set.seed(123)
  plot(b, 
       layout = layout.kamada.kawai(b),
       vertex.size = log(strength(b, mode = "total", weights = E(b)$flow)+1)/2,
       vertex.color = custom_colors[c[match(V(b)$name,c$NAME),]$role],
       vertex.frame.color = custom_colors[c[match(V(b)$name,c$NAME),]$role],
       vertex.label = NA,
       edge.arrow.size = 0.1,
       edge.width = log(E(b)$flow)/2,
       edge.color = rgb(0,0,0,((range01(E(b)$flow))*0.19)+0.01),
       # edge.label = E(b)$flow,
       # edge.label.color = "#d53e4f",
       # edge.label.cex = 0.7,
       sub = paste0(i+1989, sep=""))
}

## mblogit Models
# Model 1
node_attr <- ref_stacked
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
                            nConflicts +
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
node_attr <- ref_stacked
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
                          nConflicts +
                          conflict.sum.dead
                          # +
                          #   emdat.flood.affected
                          # +
                          # emdat.flood.deaths
                          ,
                          random = ~1|NAME,
                          data = node_attr,
                          estimator = "ML",
                          method = "PQL",
                          maxit = 100
)
summary(mblogit_model2)

# Model 3
node_attr <- ref_stacked
node_attr$role <- ifelse(node_attr$role=="3", # set desired role for model here
                         1,
                         0)
node_attr$role <- as.factor(node_attr$role)
node_attr$role <- relevel(node_attr$role, ref = "0")

mblogit_model3 <- mblogit(formula = role ~ 
                          # gdp +
                          # spei.max +
                          #   spei.min
                          # +
                          # spei.mean +
                          # spei.max.cv +
                          # spei.min.cv
                          # +
                          # spei.mean.cv +
                          nConflicts
                          +
                          conflict.sum.dead
                          # +
                          #   emdat.flood.affected
                          # +
                          # emdat.flood.deaths
                          ,
                          random = ~1|NAME,
                          data = node_attr,
                          estimator = "ML",
                          method = "PQL",
                          maxit = 100
)
summary(mblogit_model3)

### Mblogit Model Outputs
sink(file = "results/mblogit_UNHCR_roles_model_results.txt")
print("Model: 1")
summary(mblogit_model1)
print("Model: 2")
summary(mblogit_model2)
print("Model: 3")
summary(mblogit_model3)
sink(file = NULL)
