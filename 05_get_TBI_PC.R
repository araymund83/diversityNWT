# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(adespatial,crayon,dplyr,furrr, glue, future.apply, purrr, stringr,
               tictoc,tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
path <- './tables/1km_comm/allSpp/commMatrices_1km'
files <- list.files(path, pattern = 'commAll', full.names = TRUE)
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[2]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
set.seed(1234)

# # data ---------------------------------------------------------------
make_matrices <- function(gcms, yr1, yr2, files){
  lapply(gcms, function(gc){
    message(crayon::green('Making matrices form GCM:', gc))
    fls <- grep(gc, files, value = TRUE)
    fl1 <- grep(yr1, fls, value = TRUE)
    fl2 <- grep(yr2, fls, value = TRUE)
    df1 <- qs::qread(fl1)
    df1 <- df1 %>%
      mutate(pixelID = row_number()) %>%
      mutate(pixelID = glue('S_{pixelID}')) %>%
      remove_rownames() %>%
      column_to_rownames(var = 'pixelID')

    df2 <- qs::qread(fl2)
    df2 <- df2 %>%
      mutate(pixelID = row_number()) %>%
      mutate(pixelID = glue('S_{pixelID}')) %>%
      remove_rownames() %>%
      column_to_rownames(var = 'pixelID')
    names <- c(glue('{gc}_{yr1}'), glue('{gc}_{yr2}'))
    dfList<- setNames(list(df1, df2), names)
    return(dfList)
  })
}
# Apply the function ------------------------------------------------------
data <- make_matrices(gcms, yr1, yr2, files)
print(data)
qs::qsave(data, './outputs/listCommMatrices_1191_1km.qs')
data <- qs::qread('./outputs/listCommMatrices_1131_1km.qs')

# Get Temporal Beta diversity Index (TBI) ---------------------------------
#extract data frames for yr1 (before) and yr2 (after)
df1 <- map(1:length(data), function(i) data[[i]][[1]])
names(df1) <- glue('{gcms}_{yr1}')
df2 <- map(1:length(data), function(i) data[[i]][[2]])
names(df2) <- glue('{gcms}_{yr2}')


get_TBI<- function(gcms, df1, df2, yr1, yr2){
   tbi_gcm <- map(.x = 1:length(gcms), function(gc){
    message(crayon::green('Making TBI for GCM:', gc))
    tbi <- adespatial::TBI(df1, df2, method = '%diff', BCD = TRUE, nperm = 9999,
                           test.BC = TRUE, save.BC = TRUE, seed. = 1234,
                           test.t.perm = TRUE, 
                           clock = TRUE)
    
    out <- glue('./outputs/tbi_1km')
    ifelse(!file.exists(out), dir.create(out), print("Directory already exists"))
    qs::qsave(x = tbi, file = glue('{out}/tbi_{gc}_{yr1}_{yr2}2_1km.qs'))
    return(tbi)
  })
  return(tbi_gcm)  
}

# apply the function ------------------------------------------------------
#future::plan(future::cluster, workers = 30, gc = TRUE) # 'multicore' does not work with windows or Rstudio
future::plan(future::multisession, workers = 20, gc = TRUE) 
options(future.globals.maxSize= 4194304000)
tbi_gcms <- future_map(.x = gcms[[1]], df1[[1]], df2[[1]], yr1, yr2, .f = get_TBI, 
                       .options = furrr_options(seed = TRUE))
future:::ClusterRegistry('stop') 
