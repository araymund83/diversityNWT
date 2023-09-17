# Load libraries --------------------------------------------
library(pacman)

pacman::p_load(fs, fst, glue, gtools, qs, terra, stringr, tidyverse)

rm(list = ls())

# Load data ---------------------------------------------------------------
path <- './tables/1km_comm'
files <- list.files(path, pattern = 'comm',full.names = TRUE)
names <- str_split(files, pattern = '_')
species<- sapply(names, `[[`, 3)
species <- unique(species)
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
#gcms <- sapply(names, `[[`, 4)
#gcms <- unique(gcms)

# Function ----------------------------------------------------------------
join_tble <- function(gcms){
  table <- map(.x = 1:length(gcms), function(i){
    cat('Start at gcm ', gcms[i], '\n')
  fls <- grep(gcms[i], files, value = TRUE)
  fle <- gtools::mixedsort(fls)
  fle <- as.character(fle)
  tbl <- lapply(fle, qs::qread)
  allTbl <- tbl %>% reduce(full_join, by = c('x', 'y'))  # reduce:comprime la funcion para poder aplicarla a una lista 
  colnames(allTbl) <- gsub(glue('_{gcms[i]}'), '', colnames(allTbl))
 
  out <- glue('./tables/1km_comm/allSpp')
  ifelse(!file.exists(out), dir_create(out), print("already exists"))
  qs::qsave(x = allTbl, file = glue('{out}/commAllSpp_113191_{gcms[i]}_1km.qs'))
  cat('Done!\n')
 })
}

# Apply this function -------------------------------------------------------
map(.x = gcms, .f = join_tble)

