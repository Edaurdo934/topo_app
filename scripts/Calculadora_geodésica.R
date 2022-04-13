
##################################################################################################
#-----------------------------------------UTM A GEODÉSICAS----------------------------------------

## IMPORTAMOS DATOS
#--------------------
rm(list = ls())
library(sf)
library(tidyverse)

Archivo <- file.choose()
utm <- read.csv(Archivo,
               header = TRUE,
               sep=";")
utm

#Convertimos a sf (simple feature)

#Ingresar zona:
zona = 18

codigo <- if(zona == 17) {
  32717
} else {
  if (zona == 18) {
    32718
  } else {
    32719
  }
}

utm <- st_as_sf(utm, coords=c(X = "ESTE",Y = "NORTE",Z = "h"), crs=codigo, remove=F)  

#Seleccionar solo la geometría

utm <- utm %>% select(geometry)

## CONVERTIR A GEODÉSICAS: 
#------------------------

latlong <- st_transform(utm, 4326)

latlong <- st_coordinates(latlong$geom) %>% as.data.frame()
options(digits = 11)
colnames(latlong) <- c("Longitud", "Latitud", "h") 
latlong <- latlong %>% select(Latitud, Longitud, h)

latlong

readr::write_csv(x = latlong, file = "data/PAF05_latlong.csv")



##################################################################################################
#-----------------------------------------GEODÉSICAS A UTM----------------------------------------

## IMPORTAMOS DATOS
#--------------------

library(sf)
library(tidyverse)

Archivo <- file.choose()
latlong <- read.csv(Archivo,
                header = TRUE,
                sep=";")
latlong

#Convertimos a sf (simple feature)

latlong <- st_as_sf(latlong, coords=c(X = "Longitud",Y = "Latitud",Z = "h"), crs=4326, remove=F)  

#Seleccionar solo la geometria

latlong <- latlong %>% select(geometry)


## CONVERTIR A UTM 
#--------------------

#Ingresar zona:
zona = 18

codigo <- if(zona == 17) {
  32717
} else {
  if (zona == 18) {
    32718
  } else {
    32719
  }
}

utm <- st_transform(latlong, codigo)

utm <- st_coordinates(utm$geom) %>% as.data.frame()
options(digits = 11)
colnames(utm) <- c("ESTE", "NORTE", "h") 
utm <- utm %>% select(NORTE, ESTE, h)

utm
readr::write_csv(x = utm, file = "data/PAF05_UTM.csv")


##################################################################################################
#------------------------------------------PSAD56 A WGS84-----------------------------------------

## IMPORTAMOS DATOS
#--------------------

library(sf)
library(tidyverse)

Archivo <- file.choose()
PSAD <- read.delim(Archivo, sep=";", header = T) 
PSAD

#Ingresar zona:
zona = 18

codigopsad <- if(zona == 17) {
  24877
} else {
  if (zona == 18) {
    24878
  } else {
    24879
  }
}

codigowgs <- if(zona == 17) {
  32717
} else {
  if (zona == 18) {
    32718
  } else {
    32719
  }
}


PSAD <- st_as_sf(PSAD, coords=c(X = "ESTE",Y = "NORTE",Z = "h"), crs=codigopsad, remove=F)  
PSAD <- PSAD %>% select(geometry)
WGS <- st_transform(PSAD, codigowgs)
options(digits = 11)
WGS <- st_coordinates(WGS$geom) %>% as.data.frame()
colnames(WGS) <- c("ESTE", "NORTE", "h") 
WGS <- WGS %>% select(NORTE, ESTE, h)
WGS

readr::write_csv(x = WGS, file = "data/WGSconv.csv")


##################################################################################################
#------------------------------------------WGS84 A PSAD56-----------------------------------------

## IMPORTAMOS DATOS
#--------------------

library(sf)
library(tidyverse)

Archivo <- file.choose()
WGS <- read.delim(Archivo, sep=";", header = T) 
WGS

#Ingresar zona:
zona = 18

codigowgs <- if(zona == 17) {
  32717
} else {
  if (zona == 18) {
    32718
  } else {
    32719
  }
}

codigopsad <- if(zona == 17) {
  24877
} else {
  if (zona == 18) {
    24878
  } else {
    24879
  }
}


WGS <- st_as_sf(WGS, coords=c(X = "ESTE",Y = "NORTE",Z = "h"), crs=codigowgs, remove=F)  
WGS <- WGS %>% select(geometry)
PSAD <- st_transform(WGS, codigopsad)
options(digits = 11)
PSAD <- st_coordinates(PSAD$geom) %>% as.data.frame()
colnames(PSAD) <- c("ESTE", "NORTE", "h") 
PSAD <- PSAD %>% select(NORTE, ESTE, h)
PSAD

readr::write_csv(x = PSAD, file = "data/PSAD56conv.csv")

