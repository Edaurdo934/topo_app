##Calculadora geodesica

calculadoraGeodesica <- function(datos, tipo, x, y, z, crs =NULL, crs2=NULL){
  if(tipo == "utmToGeodesic"){
    codigo <- crs
    
    utm <- st_as_sf(datos, coords=c(X = x,Y = y, Z = z), crs=codigo, remove=F)  
    
    #Seleccionar solo la geometría
    utm <- utm %>% select(geometry)
    
    ## CONVERTIR A GEODÉSICAS: 
    #------------------------
    
    latlong <- st_transform(utm, 4326)
    
    latlong <- st_coordinates(latlong$geom) %>% as.data.frame()
    
    colnames(latlong) <- c("Longitud", "Latitud", "h") 
    latlong <- latlong %>% select(Latitud, Longitud, h)
  
    crsFinal <- 4326
    resultados <- list(latlong, crsFinal)
    return(resultados)
    
  } else if(tipo == "geodesicToUtm"){
    #Convertimos a sf (simple feature)
    
    latlong <- st_as_sf(datos, coords=c(X = x,Y = y,Z = z), crs=4326, remove=F)  
    print(latlong)
    #Seleccionar solo la geometria
    
    latlong <- latlong %>% select(geometry)
  
    #Ingresar zona:
    codigo <- crs
    
    utm <- st_transform(latlong, codigo)
    
    utm <- st_coordinates(utm$geom) %>% as.data.frame()
    
    colnames(utm) <- c("ESTE", "NORTE", "h") 
    utm <- utm %>% select(NORTE, ESTE, h)
    
    crsFinal <- codigo
    resultados <- list(utm, crsFinal)
    return(resultados)
    
  } else if(tipo == "psad56ToWGS84"){
    codigopsad <- crs
    
    codigowgs <- crs2
    
    
    PSAD <- st_as_sf(datos, coords=c(X = x,Y = y,Z = z), crs=codigopsad, remove=F)  
    PSAD <- PSAD %>% select(geometry)
    WGS <- st_transform(PSAD, codigowgs)
    WGS <- st_coordinates(WGS$geom) %>% as.data.frame()
    colnames(WGS) <- c("ESTE", "NORTE", "h") 
    WGS <- WGS %>% select(NORTE, ESTE, h)
    
    crsFinal <- codigowgs
    resultados <- list(WGS, crsFinal) 
    return(resultados)
  } else if(tipo == "WGS84ToPsad56"){
    codigowgs <- crs
    
    codigopsad <- crs2
    
    WGS <- st_as_sf(WGS, coords=c(X = x,Y = y,Z = z), crs=codigowgs, remove=F)  
    WGS <- WGS %>% select(geometry)
    PSAD <- st_transform(WGS, codigopsad)
    PSAD <- st_coordinates(PSAD$geom) %>% as.data.frame()
    colnames(PSAD) <- c("ESTE", "NORTE", "h") 
    PSAD <- PSAD %>% select(NORTE, ESTE, h)
    
    crsFinal <- codigopsad
    resultados <- list(PSAD, crsFinal)
    return(resultados)
  }
}  
