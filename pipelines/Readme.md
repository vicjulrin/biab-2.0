Flujo de trabajo - Indicadores Biotablero
================
Instituto de Investigación de Recursos Biológicos Alexander von
Humboldt  
2022-11

``` r
## Nombrar librerias ####
packagesList<-list("dplyr", "plyr",  "raster", "terra", "rstac", "gdalcubes", "sf", "stars", "pbapply", "rasterDT", "rstudioapi")

## Cargar librerias ####
lapply(packagesList, library, character.only = TRUE)


## Establecer directorio ##
dir_file<- getSourceEditorContext()$path
dir_folder<- dirname(dir_file)

setwd(dir_folder)

## Cargar funciones
source("D:/IAvH 2022/STAC/funciones.R")

## tabla de metadatos
metadata<- read.csv("metadata covers stac.csv")


## Definir dirección para conexion con repositorio
repository_stac <- rep_connect_bridge(repository_url = "http://io.biodiversite-quebec.ca/stac/")
repository_stac

## Revisar colecciones disponibles en repositorio
repository_collections<- check_collections(dir_stac= repository_stac)
repository_collections



## Establecer conexión con colecciones
colections_stac<- load_collections(dir_stac= repository_stac,
                                   collection= c("gfw-treecover2000", "chelsa-clim" )  )
colections_stac

## Revisar colecciones disponibles por capa
repository_layers<- check_layers(collection_object= colections_stac)
repository_layers

## Generar colección de capas
stac_colections_layers<- load_layer_collection(collection_object= colections_stac,
                                               collection_layers = list(
                                                 "gfw-treecover2000"= c(""), "chelsa-clim"= c("bio9", "bio8"))
)

## Estimación por área de estudio
test1<- stac_layer(dir_polygon = "C:/WCS/contrafactual/defaultlayers/Departamento/Amazonas.shp",
                   colections_layers= stac_colections_layers,
                   res= 1000
                   )

```

