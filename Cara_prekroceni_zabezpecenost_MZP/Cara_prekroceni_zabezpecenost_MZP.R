library(tidyverse)
MZP <-  function(data) {
  data <- data %>%
    mutate(year = year(DTM),
           month = month(DTM)) %>%
    filter(month == 8)
  
  MZP_df <- tibble(
    ID = numeric(),
    typ_Q = character(),
    MZP_rak = numeric()
  )
  
  for (id in unique(data$ID)) {
    data_perc <- data %>% filter(ID == id)
    
    MZP_bil <- quantile(data_perc$Q_bil_souc, probs = 0.05, na.rm = TRUE)
    MZP_gr4j <- quantile(data_perc$Q_GR4J_souc, probs = 0.05, na.rm = TRUE)
    
    MZP_df <- MZP_df %>%
      add_row(ID = id, typ_Q = "Q_bilan", MZP_rak = MZP_bil) %>%
      add_row(ID = id, typ_Q = "Q_GR4J", MZP_rak = MZP_gr4j)
  }
  
  return(MZP_df)
}

mzp_df <- MZP(dataMZP)  

zabezpec <- list()
for (typemod in unique(data_dif$type)) {
  data_dif %>% 
    mutate(PER = factor(PER, levels = c("NEAR", "MID", "FAR"))) %>% 
    filter(type == typemod) -> data_cmip
  
  for(id in unique(data_cmip$ID)){
    # dir.create(paste0("GRAFY/cara_prekroceni", typemod, "/", id), recursive = TRUE, showWarnings = FALSE)
    data_cmip %>% 
      filter(ID == id) -> data_id
    
    mzp_df %>% 
      filter(ID == id) -> mzp_id
    
    for (hydromodel in c("Q_bilan", "Q_GR4J") ) {
      data_id %>% 
        pivot_longer(cols = c(Q_bilan, Q_GR4J, Qobs),
                     names_to = "model_obs",
                     values_to = "Q_all") %>% 
        select(DTM, ID, Q_all, model_obs, SID, SCE, type, PER) %>% 
        filter(model_obs == hydromodel) -> data_long
      
      mzp_id %>% 
        filter(typ_Q == hydromodel) -> mzp      
      
      data_long %>% 
        group_by(ID, model_obs, SID, SCE, PER, type) %>% 
        arrange(desc(Q_all)) %>% 
        mutate(poradi = row_number(),
               n = n(),
               pravdep = poradi / (n+1)) %>% 
        ungroup() -> data_plot
      key <- paste(id, hydromodel, typemod, sep = "-")
      zabezpec[[key]] <- data_plot
      
      
      mzp_value <- mzp_id %>% 
        filter(typ_Q == hydromodel) %>% 
        pull(MZP_rak) %>% 
        as.numeric()
      
      data_plot %>%
        group_by(SID, SCE, type, PER) %>% 
        mutate(mzp = ifelse(Q_all < mzp_value, 1, 0)) %>% 
        filter(mzp == 0) %>%  
        filter(pravdep == max(pravdep)) -> prekroc
      
      
      p <- ggplot(data = data_plot)+
        geom_line(aes(x = pravdep, y = Q_all, color = interaction(SID, model_obs), group = interaction(SID, model_obs)),
                  alpha = 0.6)+
        geom_hline(yintercept = mzp_value, color = "red")+
        geom_vline(data = prekroc, aes(xintercept = pravdep, color = interaction(SID, model_obs)))+
        annotate("text", x = min(data_plot$pravdep), y = mzp_value, 
                 label = round(mzp_value, 3), 
                 vjust = -0.5, hjust = 0, color = "red")+
        geom_text(
          data = prekroc,
          aes(x = pravdep, y = 2000,  
              label = round(pravdep, 3),
              color = interaction(SID, model_obs)),
          vjust = -0.5,  
          angle = 90,    
          size = 2
        )+
        labs(x = "Pravděpodobnost",              
             y = "Průtok [m³/s]")+
        labs(color = NULL)+
        ggtitle(paste("Čára překročení pro minimální zůstatkový průtok. Povodí: ", id, "  Generace: ", typemod))+
        scale_y_log10()+
        facet_grid(SCE~PER)+
        theme_light()
      ggsave(paste0("GRAFY/zabezpecenost_vse/", id, "_", typemod, "_", hydromodel, ".png"), plot = p, width = 15, height = 10, dpi = 300)
      
    }
  }
}


write_rds(zabezpec, "zabezpecenost.rds")

write_csv(data_dif, "FINAL_dif.csv")

data_dif <- read_csv("FINAL_dif.csv")
#cara prekroceni - mesice
zabezpec_month <- list()
for (typemod in unique(data_dif$type)) {
  data_dif %>% 
    mutate(PER = factor(PER, levels = c("NEAR", "MID", "FAR"))) %>% 
    filter(type == typemod) -> data_cmip
  
  for(id in unique(data_cmip$ID)){
    # dir.create(paste0("GRAFY/cara_prekroceni", typemod, "/", id), recursive = TRUE, showWarnings = FALSE)
    data_cmip %>% 
      filter(ID == id) -> data_id
    
    mzp_df %>% 
      filter(ID == id) -> mzp_id
    
    for (hydromodel in c("Q_bilan", "Q_GR4J") ) {
      data_id %>% 
        pivot_longer(cols = c(Q_bilan, Q_GR4J, Qobs),
                     names_to = "model_obs",
                     values_to = "Q_all") %>% 
        select(DTM, ID, Q_all, model_obs, SID, SCE, type, PER) %>% 
        filter(model_obs == hydromodel) -> data_long
      
      mzp_id %>% 
        filter(typ_Q == hydromodel) -> mzp      
      
      data_long %>% 
        mutate(month = month(DTM)) -> data_long_month
      
      for (mon in unique(data_long_month$month)) {
        dir.create(paste0("GRAFY/zabezpecenost_mesice", "/", mon), recursive = TRUE, showWarnings = FALSE) 
        data_long_month %>% 
          filter(month == mon) %>% 
          group_by(ID, model_obs, SID, SCE, PER) %>% 
          arrange(desc(Q_all)) %>% 
          mutate(poradi = row_number(),
                 n = n(),
                 pravdep = poradi / (n+1)) %>% 
          ungroup() -> data_plot
        key <- paste(id, hydromodel, typemod, mon, sep = "-")
        zabezpec_month[[key]] <- data_plot
        
        data_plot %>%
          group_by(SID, SCE, type, PER) %>%
          mutate(mzp = ifelse(Q_all < mzp$MZP_rak, 1, 0)) %>%
          filter(mzp == 0) %>%
          filter(pravdep == max(pravdep)) -> prekroc
        
        
        # mzp_value <- as.numeric(mzp$MZP_rak[1])
        p <- ggplot(data = data_plot)+
          geom_line(aes(x = pravdep, y = Q_all, color = interaction(SID, model_obs), group = interaction(SID, model_obs)),
                    alpha = 0.6)+
          geom_hline(yintercept = mzp_value, color = "red")+
          geom_vline(data = prekroc, aes(xintercept = prekroc$pravdep, color = interaction(SID, model_obs)))+
          annotate("text", x = min(data_plot$pravdep), y = mzp_value,
                   label = round(mzp_value, 3),
                   vjust = -0.5, hjust = 0, color = "red")+
          geom_text(
            data = prekroc,
            aes(x = pravdep, y = 2000,  # nahoře v plotu
                label = round(pravdep, 3),
                color = interaction(SID, model_obs)),
            vjust = -0.5,  # posune text mírně nad vrchol
            angle = 90,    # otočí text, aby stál svisle
            size = 2
          )+
          labs(x = "Pravděpodobnost",
               y = "Průtok [m³/s]")+
          labs(color = NULL)+
          ggtitle(paste("Čára překročení pro minimální zůstatkový průtok. Povodí: ", id, "  Generace: ", typemod, "Měsíc: ", mon))+
          scale_y_log10()+
          facet_grid(SCE~PER)+
          theme_light()
        
        
        
        ggsave(paste0("GRAFY/zabezpecenost_mesice/", mon, "/", id, "_", typemod, "_", hydromodel, "_", mon, ".png"), plot = p, width = 15, height = 10, dpi = 300)
      }
    }
  }
}

write_rds(zabezpec_month, "zabezpecenost_mesice.rds")



#zabezpecenost - tabulka
zabezpec <- read_rds("zabezpecenost.rds")
df_zabezpec <- bind_rows(zabezpec)

mzp_df %>% 
  rename(model_obs = typ_Q) -> mzp_df

df_zabezpec %>% 
  left_join(mzp_df, by = c("ID", "model_obs")) %>% 
  group_by(ID, model_obs, SID, SCE, type, PER) %>% 
  filter(Q_all > MZP_rak) %>% 
  slice_min(Q_all, with_ties = FALSE)%>% 
  mutate(SID = factor(SID),
         SCE = factor(SCE))-> result


draw_heatmap <- function(data, per_val, model_val) {
  data %>%
    filter(PER == per_val, model_obs == model_val) %>%
    mutate(combo = interaction(factor(SID), factor(SCE), factor(type), sep = "_")) %>%
    ggplot(aes(x = factor(ID),
               y = factor(combo, levels = rev(levels(factor(combo)))),
               fill = pravdep,
               label = round(pravdep, 2))) + 
    geom_tile(color = "white")+
    geom_text(
      data = data %>%
        filter(PER == per_val, model_obs == model_val) %>%
        mutate(combo = interaction(factor(SID), factor(SCE), factor(type), sep = "_")) %>%
        filter(ID != "1007"),
      aes(label = round(pravdep, 2)),
      size = 3,
      color = "black"
    )+
    geom_text(
      data = data %>%
        filter(PER == per_val, model_obs == model_val) %>%
        mutate(combo = interaction(factor(SID), factor(SCE), factor(type), sep = "_")) %>%
        filter(ID == "1007"),
      aes(label = round(pravdep, 2)),
      size = 3,
      color = "red"
    )+
    scale_fill_gradient(
      low = "red", high = "green", na.value = "grey90"
    ) +
    labs(
      title = paste("PER:", per_val, "| Model:", model_val),
      x = element_blank(), y = "model", fill = "Pravdep"
    ) +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(size = 10),
          axis.text.y = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
    )
}


PER_levels <- c("NEAR", "MID", "FAR")
model_levels <- unique(result$model_obs)

plots <- list()

for (per in PER_levels) {
  for (model in model_levels) {
    key <- paste(per, model, sep = "_")
    plots[[key]] <- draw_heatmap(result, per, model)
  }
}

plots

for (i in seq_along(plots)) {
  ggsave(
    filename = paste0("GRAFY/zabezpecenost_tabulka/", names(plots)[i], ".png"),
    plot = plots[[i]],
    width = 7, height = 5, dpi = 300,
    units = "in"
  )
}




# zabezpecenost - mesice - tabulka
setwd("/var/data/temp/DuDy/vyhodnoceni_DUDY/")
zabezpec_month <- read_rds("zabezpecenost_mesice.rds")

mzp_df <- read_csv("")

df_zabezpec_mes <- bind_rows(zabezpec_month)

mzp_df %>% 
  rename(model_obs = typ_Q) -> mzp_df

df_zabezpec_mes %>% 
  left_join(mzp_df, by = c("ID", "model_obs")) %>% 
  group_by(ID, model_obs, SID, SCE, type, PER, month) %>% 
  filter(Q_all > MZP_rak) %>% 
  slice_min(Q_all, with_ties = FALSE)%>% 
  mutate(SID = factor(SID),
         SCE = factor(SCE))-> result_mes


draw_heatmap_mesic <- function(data, per_val, model_val, mesic) {
  data %>%
    filter(PER == per_val, model_obs == model_val, month == mesic) %>%
    mutate(combo = interaction(factor(SID), factor(SCE), factor(type), sep = "_")) %>%
    ggplot(aes(x = factor(ID),
               y = factor(combo, levels = rev(levels(factor(combo)))),
               fill = pravdep,
               label = round(pravdep, 2)))+  
    geom_tile(color = "white") +
    geom_text(
      data = data %>%
        filter(PER == per_val, model_obs == model_val, month == mesic) %>%
        mutate(combo = interaction(factor(SID), factor(SCE), factor(type), sep = "_")) %>%
        filter(ID != "1007"),
      size = 3,
      color = "black")+
    geom_text(data = data %>%
                filter(PER == per_val, model_obs == model_val, month == mesic) %>%
                mutate(combo = interaction(factor(SID), factor(SCE), factor(type), sep = "_")) %>%
                filter(ID == "1007"),
              aes(label = round(pravdep, 2)),
              size = 3,
              color = "red")+
    scale_fill_gradient(
      low = "red", high = "green", na.value = "grey90"
    ) +
    labs(
      title = paste("PER:", per_val, "| Hyd. Model:", model_val, "| Měsíc:", mesic),
      x = "ID povodí", y = "Model", fill = "Pravdep"
    ) +
    theme_light(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
}