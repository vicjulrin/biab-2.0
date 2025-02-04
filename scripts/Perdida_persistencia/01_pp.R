### funci�n check_collections

# Cargar librerias
packages_list<-list("magrittr")
invisible(lapply(packages_list, library, character.only = TRUE))

# Organizar directorios
args <- commandArgs(trailingOnly=TRUE)
outputFolder <- args[1]

# Cargar archivos de entrada
input <- rjson::fromJSON(file=file.path(outputFolder, "input.json"))

dir_wkt<- input$dir_wkt_polygon
epsg_polygon<- input$epsg_polygon
resolution<- input$resolution
dir_colection<- input$dir_colection
folder_output<- input$folder_output

wkt_polygon<- readLines(dir_wkt)
layers <- list.files(dir_colection, "\\.tif$", recursive = TRUE, full.names = TRUE)
json_colleciton_file <- list.files(dir_colection, "\\.json$", recursive = TRUE, full.names = TRUE)
meadata_collecion_file <- list.files(dir_colection, "\\.csv$", recursive = TRUE, full.names = TRUE)
metadata<- read.csv(meadata_collecion_file)

# Especificar carpea donmde guardar resultados
folder_results<- paste0(dirname(dirname(dirname(outputFolder))), "/", folder_output)
dir.create(folder_results)


# Especificar info area de estudio
vector_polygon<- terra::vect(wkt_polygon, crs= sf::st_crs(epsg_polygon)$proj4string ) %>% terra::as.polygons()
crs_polygon<- terra::crs(vector_polygon)
box_polygon<-  sf::st_bbox(vector_polygon)


# Alinear con colecci�on
stac_collection<- gdalcubes::create_image_collection(files= layers, format= json_colleciton_file) 

# Cargar cubo
cube_collection<- gdalcubes::cube_view(srs = crs_polygon,  extent = list(t0 = gdalcubes::extent(stac_collection)$t0, t1 = gdalcubes::extent(stac_collection)$t1,
                                                left = box_polygon[1], right = box_polygon[3],
                                                top = box_polygon[4], bottom = box_polygon[2]),
              dx = resolution, dy = resolution, dt = "P1Y", aggregation = "first", resampling = "first", keep.asp= F)

cube <- gdalcubes::raster_cube(stac_collection, cube_collection)

# Cortar cubo por area de estudio
cube_mask<- gdalcubes::filter_geom(cube, geom= wkt_polygon, srs = crs_polygon )


# Convertir cubo a raster
cube_stars <- stars::st_as_stars(cube_mask) %>% terra::rast()
collection_rast<- lapply(cube_stars, function(x) { if(any( is.na(summary(raster::raster(x))) )){NULL}else{x} } ) %>%
  {Filter(function(x) !is.null(x), .)} %>% {setNames(., unlist(sapply(., function(x) names(x))) )} %>% terra::rast()

# estimar metricas de area
data_sum<- terra::freq(collection_rast, usenames=T) %>% dplyr::mutate(area_m2= count*resolution) %>% dplyr::select(-count) %>%
  dplyr::rename(collection= layer) %>% dplyr::mutate(layer= sapply(.$collection, function(x) stringr::str_split(x, "_B*time")[[1]][1]) ) %>%
  list(metadata) %>% plyr::join_all()

# guardar rasters
setwd(folder_results)
saveraster<- lapply(collection_rast, function(x) terra::writeRaster(x, paste0(names(x), ".tif"), overwrite=T) )

# guardar tabla resumen
setwd(folder_results)
write.csv(data_sum,"data_summ.csv")

# Imprimir resultado en logs
print(data_sum)

## Imprimir resultado - Formao json
output <- list("folder_output" = folder_results)
jsonData <- rjson::toJSON(output, indent=2)
write(jsonData, file.path(outputFolder,"output.json"))