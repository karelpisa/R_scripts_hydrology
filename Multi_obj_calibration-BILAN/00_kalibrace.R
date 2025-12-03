library(bilan)
library(tidyverse)
library(fst)
library(hydroGOF)
library(data.table)
setwd("Multi_obj_calibration-BILAN/")

data <- read_csv("AT_final.csv") 



kalibrace <- function(data){  # kompletní kalibrační funkce
  s1 <- seq(1, 9, 1) #proměnné, které pak zastupují poměry v kalibračních strategiích
  s2 <- seq(9, 1, -1)
  strategies <- list()
  for (k in 1:9) {    #vytvoření kalibračních poměrů ET:R
    strategies[[k]] <- c(s1[k], s2[k])
  }
  for (strategy in strategies) {
  result_kalib <- list() #list výsledků
  params_kalib <- list() #list parametrů
  best_params <- list()
  for (i in seq_along(unique(data$IDiga))) {        
    id = unique(data$IDiga)[i]
    data |>         #vyfiltrování jedné stanice
      filter(IDiga == id) -> data_id
    
    data_id |>      #data do bilanu
      dplyr::select(DTM, R, P, T) -> kalibrace_id
  
    data_id |>      #DF evapotranspirace pro kalibraci
      filter(IDiga == id) |> 
      dplyr::select(DTM, ET_mmday) -> ET
  
    data_id |>      #DF odtoku pro kalibraci
      filter(IDiga == id) |> 
      dplyr::select(DTM, R) -> R
  
    b <- bil.new("d", data = kalibrace_id, modif = "critvars")
  
    bil.pet(b, pet_type = "latit", latitude = 48)
  
    bil.set.params.lower(b, list(Grd = 0.001))
    bil.set.params.lower(b, list(Spa = 30))
    bil.set.params.upper(b, list(Spa = 300))
    bil.set.params.upper(b, list(Alf = 0.65))
    
    r_optim <- function(var){     #funkce multi obj kalibrace (odtok)
      goff_r <- -1*hydroGOF::KGE(sim = var,
                            obs = R |> pull(R))
      goff_r
    }

    et_optim <- function(var){    #funkce multi obj kalibrace (ET)
      goff_et <- -1*hydroGOF::KGE(sim = var, #výpočet goff
                          obs = ET |> pull(ET_mmday))
      goff_et
      }

  bil.set.critvars(b,                        #nastavení kritérií kalibrace modelu
                   weights = c(strategy[1], strategy[2]),      #váhy kalibrujících proměnných
                   obs_vars = c("R", "ET"),  #proměnné z modelu, které jsou měřené (ale moc nechápu důvod tohodle parametru, ale je povinný)
                   obs_values = c(-1, -1),    #optimální hodnota výstupu funkcí definovaných nahoře
                   mod_vars = c("RM", "ET"),  #modelované veličiny (var) -> vstup do funkcí definovaných nahoře
                   crits = c("custom", "custom"),  #obě funkce jsou custom, tedy definované výše
                   funs = c(r_optim, et_optim))     #dkaz na ty konkrétní funkce
  
  bil.set.optimDE(b, method = "DE", ens_count = 20)     #metoda optimalizace
  result <- bil.optimize(b)           #optimalizace
  result |>                           #manipulace s výstupem kalibrace
    mutate(IDiga = id,
           KGE = hydroGOF::KGE(sim = result$RM,
                               obs = result$R),
           NSE = hydroGOF::NSE(sim = result$RM,
                               obs = result$R)) |> 
    left_join(ET, by = "DTM") -> result
  
  result_kalib[[paste0(id)]] <- result      #přidání do listu
    
  params <- bil.get.params(b)               #hodnota parametrů
  params_kalib[[paste0(id)]] <- params |>   #přidání do litu podle ID
    mutate(IDiga = id) |> 
    select(name, current, IDiga)
  best <- bil.get.ens.resul(b)
  best_params[[id]] <- best |> 
    mutate(IDiga = id)
  print(paste("Povodi", id, "Kalibrace Dokončena"))
  
  
  b2 <- bil.new(type = "d", data = kalibrace_id, modif = "critvars")
  bil.pet(b2, pet_type = "latit", latitude = 48)
  
  ensemble_list <- list()
  dir.create(paste0("strategie/", "R", strategy[1], "_ET", strategy[2]))
  dir.create(paste0("strategie/", "R", strategy[1], "_ET", strategy[2], "/", id))
  write_csv(best_params[[id]] |> mutate(IDiga = id), file = paste0("strategie/", "R", strategy[1], "_ET", strategy[2], "/", id, "/best_params_", id, ".csv"))
  for (j in 1:20) {
  # j = 20
  bil.set.params.curr(b2, params = list(Spa = best_params[[id]]$Spa[j],
                                    Alf = best_params[[id]]$Alf[j],
                                    Dgm = best_params[[id]]$Dgm[j],
                                    Soc = best_params[[id]]$Soc[j],
                                    Mec = best_params[[id]]$Mec[j],
                                    Grd = best_params[[id]]$Grd[j]))
  ensemble_run <- bil.run(b2)
  ensemble_run |> 
    left_join(ET, by = "DTM") |> 
    mutate(IDiga = id,
           KGE_R = hydroGOF::KGE(sim = RM,
                               obs = R),
           NSE_R = hydroGOF::NSE(sim = RM,
                               obs = R),
           KGE_ET = hydroGOF::KGE(sim = ET,
                                  obs = ET_mmday),
           NSE_ET = hydroGOF::NSE(sim = ET,
                                  obs = ET_mmday),
           ens_no = j) -> ensemble_result
  ensemble_list[[j]] <- ensemble_result
  
  print(paste("ensemble", j, "povodi", id, "OK"))
  }
  ensemble_id <- bind_rows(ensemble_list)
  write_csv(ensemble_id, file = paste0("strategie/", "R", strategy[1], "_ET", strategy[2], "/", id, "/", id, ".csv" ))
  print(paste("povodi", id, "OK"))
  }
  # params_kalib_final <- bind_rows(params_kalib)   #spojení do long formátu
  # write.fst(params_kalib_final, "00_SaveWater_domena/kalibrace/01_kalibrace/parametry_kalibrace.fst") #uložení
  
  # result_kalib_final <- bind_rows(result_kalib)   #spojení do long formátu
  # write.fst(result_kalib_final, "00_SaveWater_domena/kalibrace/01_kalibrace/kalibrace.fst")   #uložení
  
  # assign("params", params_kalib_final, envir = .GlobalEnv) #uložení i do globálního envrionmentu v R
  # assign("results", result_kalib_final, envir = .GlobalEnv) #uložení i do globálního envrionmentu v R
  # assign("best5params", best_parms_final, envir = .GlobalEnv) #uložení i do globálního envrionmentu v R
  }
  print(paste("Strategie", strategy, "OK"))
}




kalibrace(data)



  






