check_daily_sequence <- function(dates) {
  dates <- as.Date(dates)
  
  dates <- sort(dates)
  
  diffs <- diff(dates)
  
  # Pokud jsou všechny rozdíly 1 den
  if (all(diffs == 1)) {
    message("Data jsou po dnech, nic nechybí.")
    return(TRUE)
  } else {
    # Najdeme chybějící dny
    missing_days <- unlist(lapply(which(diffs > 1), function(i) {
      seq(dates[i] + 1, dates[i + 1] - 1, by = "day")
    }))
    warning("Některé dny chybí!")
    return(as.Date(missing_days))
  }
}