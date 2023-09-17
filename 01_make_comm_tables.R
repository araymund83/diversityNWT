# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load( dplyr,fs, fst, glue, qs, terra, stringr, tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)


# Load data ---------------------------------------------------------------
path <- './outputs/1km'
files <- files <- list.files(path, pattern = '.tif', full.names = TRUE)
species <- str_split(files, pattern = '_')
species <- lapply(species, `[[`, 2)
species <- unique(species)
# species <- c("ALFL", "AMCR", "AMRE", "AMRO", "ATSP", "ATTW", "BAWW", "BBWA", 
#              "BBWO", "BCCH", "BHCO", "BHVI", "BLPW", "BOCH", "BRBL", "BRCR",
#              "CAWA", "CCSP", "CHSP", "CONW", "CORA", "COYE", "DEJU", "EAKI", 
#              "EVGR", "FOSP", "GCKI", "GCTH", "GRAJ", "HAFL", "HETH", "HOLA", 
#              "LCSP", "LEFL", "LEYE", "LISP", "MAWA", "NOWA", "OCWA", "OSFL",
#              "OVEN", "PAWA", "PHVI", "PISI", "PUFI", "RBGR", "RBNU", "RCKI",
#              "REVI", "RUBL", "RUGR", "RWBL", "SAVS", "SOSP", "SWSP", "SWTH", 
#              "TEWA" ,"TRES", "VATH", "WAVI", "WCSP", "WETA", "WEWP", "WIWA",
#              "WIWR", "WTSP", "WWCR", "YBFL", "YBSA", "YEWA", "YRWA")

# Make difference ---------------------------------------------------------
make_community <- function(spc){
  #spc <- species[1]
  cat('Start ', spc, '\n')
  fls <- grep(spc, files, value = TRUE)
  fls <- grep(paste0(c('2011', '2031', '2091'), collapse = '|'), fls, value = TRUE)
  names <- str_split(fls, pattern = '_')
  gcm <- c('CanESM2', 'CCSM4', 'INM-CM4')
  # gcm<- sapply(names, `[[`, 5)
  # gcm <- unique(gcm)
  yrs<- sapply(names, `[[`, 3)
  yrs <- unique(yrs)
  
  table <- map(.x = 1:length(gcm), function(i){
    cat('gcm ', gcm[i], '\n')
    fl <- grep(gcm[i], fls, value = TRUE)
    stk <- terra::rast(fl)
    tbl <- terra::as.data.frame(stk, xy = TRUE)
    names(tbl) <- c('x','y',glue('2011_{spc}_{gcm[i]}'),
      glue('2031_{spc}_{gcm[i]}'),
      glue('2091_{spc}_{gcm[i]}')
    )
    
    out <- glue('./tables/1km_comm/')
    ifelse(!file.exists(out), dir_create(out), print('Already exists'))
    qs::qsave(tbl, file = glue('{out}/comm_{spc}_{gcm[i]}.qs'))
    cat('Done!\n')
  })
  
  cat('----- Finish ----- \n')
  
}

# Apply the function ------------------------------------------------------
comm <- map(species, make_community)
