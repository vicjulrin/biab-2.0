Flujo de trabajo - Indicadores Biotablero
================
Instituto de Investigación de Recursos Biológicos Alexander von
Humboldt  
2022-11

``` r
# Definir dirección de repositorio
dir_stac<- "http://io.biodiversite-quebec.ca/stac/"

# Establecer conección con repositorio
stac_link<- stac(dir_stac)

# Lista de colecciones disponibles
dir_data <- stac_link %>% collections() %>% get_request()
sapply(dir_data$collections, function(x) x$id)

```

