# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(dplyr, exactextractr, fasterize, fs, glue, gtools, qs, sf, stringr,
               terra, tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
path <- './tables/1km_comm/allSpp'
files <- list.files(path, pattern = '113191',full.names = TRUE)
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
yrs <- c('2011', '2031', '2091')

table_yrs <- function(gc){
  dfYrs <- map(.x = 1:length(gcms),.f = function(gc, yr){
  message(crayon::green('Making table form GCM:', gcms[gc]))
  fl <- grep(gcms[gc], files, value = TRUE)
  df <- qs::qread(fl)
  tbl <- map(.x = 1:length(yrs), function(yr) {
    cat('Start at yr ', yrs[yr], '\n')
    tblYr <- dplyr::select(df, starts_with(yrs[yr]))
    tblYr <- tblYr %>%  mutate(pixelID = row_number())
    tblYr <- tblYr %>%   mutate(pixelID = glue('S_{pixelID}'))
    colnames(tblYr) <- gsub(glue('{yrs[yr]}_'), '', colnames(tblYr))
    tblYr <- as_tibble(tblYr)
    tblYr <- tblYr %>% remove_rownames()
    tblYr<-  column_to_rownames(tblYr, var = 'pixelID')
    
    out <- glue('./tables/1km_comm/allSpp/commMatrices_1km')
    ifelse(!file.exists(out), dir_create(out), print("Directory already exists"))
    qs::qsave(x = tblYr, file = glue('{out}/commAllSpp_{yrs[yr]}_{gcms[gc]}_1km.qs'))
    cat('Done!\n')
    return(tbl)
   })
 })
}
# Apply the function ------------------------------------------------------
map(.x = gcms,  .f = table_yrs)

