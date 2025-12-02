library(bilan)
library(tidyverse)
library(fst)
library(hydroGOF)
library(data.table)

#CZUmac
setwd("/Users/karelpisa/Library/CloudStorage/OneDrive-CZUvPraze/IGA/project/")
cesta <- "/Users/karelpisa/Library/CloudStorage/OneDrive-CZUvPraze/IGA/project/data"

#PC VUV
setwd("C:/Users/karel.pisa/OneDrive - CZU v Praze/IGA/project/")
cesta <- "C:/Users/karel.pisa/OneDrive - CZU v Praze/IGA/project"



  
data <- read_csv("00_SaveWater_domena/inputs_model/AT_final.csv")


# i = 1 #nastaveni workflow pro jetí mimo for cyklus


kalibrace <- function(data){  # kompletní kalibrační funkce
  # setwd(cesta)
  result_kalib <- list() #list výsledků (podle DBCN)
  params_kalib <- list() #list parametrů (podle DBCN)
  # for (i in seq_along(unique(data$IDiga))) {          ########řádek pro vše
  for (i in seq_along(unique(data$IDiga))[1:3]) {            ########řádek pro nastavení
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
  
    bil.pet(b, )
  
    bil.set.params.lower(b, list(Grd = 0.001))
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
                   weights = c(1, 2),      #váhy kalibrujících proměnných
                   obs_vars = c("R", "ET"),  #proměnné z modelu, které jsou měřené (ale moc nechápu důvod tohodle parametru, ale je povinný)
                   obs_values = c(-1, -1),    #optimální hodnota výstupu funkcí definovaných nahoře
                   mod_vars = c("RM", "ET"),  #modelované veličiny (var) -> vstup do funkcí definovaných nahoře
                   crits = c("custom", "custom"),  #obě funkce jsou custom, tedy definované výše
                   funs = c(r_optim, et_optim))     #dkaz na ty konkrétní funkce
  
  bil.set.optim(b, method = "DE")     #metoda optimalizace
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
    mutate(ID = id) |> 
    select(name, current, ID)

  }
  params_kalib_final <- bind_rows(params_kalib)   #spojení do long formátu
  write.fst(params_kalib_final, "00_SaveWater_domena/out/parametry_kalibrace.fst") #uložení
  
  result_kalib_final <- bind_rows(result_kalib)   #spojení do long formátu
  write.fst(result_kalib_final, "00_SaveWater_domena/out/kalibrace.fst")   #uložení
  
  assign("params", params_kalib_final, envir = .GlobalEnv) #uložení i do globálního envrionmentu v R
  assign("results", result_kalib_final, envir = .GlobalEnv) #uložení i do globálního envrionmentu v R
  }


kalibrace(data)


unique(results$NSE)
unique(results$KGE)


results |> 
  filter(KGE > 0.6) -> results_good


results|>
  select(DTM, R, RM, IDiga) |> 
  pivot_longer(cols = -c(DTM, IDiga),
  names_to = "variable",
  values_to = "value") -> plot_data_r
  






