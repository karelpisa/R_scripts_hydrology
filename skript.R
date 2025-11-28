# extrahování shp filů z Datasetu HERA, nebo jiného rasteru
#paralelizace, běhá na rychlem koni

library(parallel)
library(terra)
library(sf)
library(exactextractr)
library(tidyverse)
library(data.table)
library(fst)
library(future)
library(future.apply)



setwd("/var")
rast_path <- "data/HERA/" 
shp_cr <- "../home/karelpisa/Dokumenty/sawewater/shp/CZ_povodi_stanice.shp"
shp_at <- "../home/karelpisa/Dokumenty/sawewater/shp/AT_catchments_stations_ETRS.shp"
# staty <-  grep("\\.shp$", list.files(shp_path), value = TRUE)[2] #zatím jen pov. 1. řádu v ČR...čekám na rakouskou stranu
folders <- list.dirs(rast_path, recursive = FALSE, full.names = FALSE)[2:6]
years <- seq(1990, 2020, 1)

#CESKO
data_extract <- function(year) {
  result_list <- list()
  shp <- st_read(shp_cr)
    for (folder in folders) {
      rast <- terra::rast(paste0(rast_path, folder, "/", folder, "_", year, ".nc"))
      shp_transformed <- st_transform(shp, st_crs(rast)) 
      vars <- exact_extract(rast, shp, fun = "mean", append_cols = "chp_4")
      
      vars %>%  
        as_tibble() %>% 
        pivot_longer(
          cols = -chp_4,
          names_to = "layer",
          values_to = paste0(folder)) %>% 
        group_by(chp_4) %>% 
        mutate(DTM = as.POSIXct(time(rast))) %>%
        ungroup() %>% 
        select(-layer) %>% 
        pivot_longer(cols = -c(DTM, chp_4), #pak upravit na rakouske SHP
                     names_to = "var", 
                      values_to = "value")-> df_year
      
      result_list[[paste0(folder, "_", year)]] <- df_year
      
    }
    
  df_FINAL <- bind_rows(result_list)
  return(df_FINAL)
}
  
plan(multisession, workers = 10)
list_results <- future_lapply(years, data_extract, future.seed = NULL)

FINAL <- bind_rows(list_results)
write_fst(FINAL, "../home/karelpisa/Dokumenty/sawewater/out/CZ_clima_HERA_stations_catchments_6hours.fst")

FINAL %>%
  mutate(DTMd = as.Date(DTM)) %>% 
  group_by(chp_4, DTMd, var) %>% 
  reframe(value = mean(value)) -> final_daily

write_fst(final_daily, "../home/karelpisa/Dokumenty/sawewater/out/CZ_clima_HERA_stations_catchments_daily.fst")

# RAKOUSKO
folder <- folders[1]
year <- 1990
data_extract <- function(year) {
  result_list <- list()
    shp <- st_read(shp_at)
    for (folder in folders) {
      rast <- terra::rast(paste0(rast_path, folder, "/", folder, "_", year, ".nc"))
      shp_transformed <- st_transform(shp, st_crs(rast)) 
      vars <- exact_extract(rast, shp_transformed, fun = "mean", append_cols = c("HZB")) 
      
      vars %>%  
        as_tibble() %>% 
        # distinct(c(DBCN, CHP_4), .keep_all = TRUE) %>% 
        pivot_longer(
          cols = -HZB,
          names_to = "layer",
          values_to = paste0(folder)) %>%
        group_by(HZB) %>% 
        mutate(DTM = as.POSIXct(time(rast))) %>% 
        select(-layer) %>% 
        pivot_longer(cols = -c(DTM, HZB), #pak upravit na rakouske SHP
                     names_to = "var", 
                     values_to = "value") -> df_year
      
      result_list[[paste0(folder, "_", year, "_")]] <- df_year

    }
    df_FINAL <- bind_rows(result_list)
    return(df_FINAL)
}


plan(multisession, workers = 10)
list_results <- future_lapply(years, data_extract, future.seed = NULL)

final <- bind_rows(list_results)
write_fst(final, "../home/karelpisa/Dokumenty/sawewater/out/AT_clima_6h.fst")

final %>%
  mutate(DTMd = as.Date(DTM)) %>% 
  group_by(HZB, DTMd, var) %>% 
  reframe(value = mean(value)) -> final_daily

write_fst(final_daily, "../home/karelpisa/Dokumenty/sawewater/out/AT_clima.fst")




