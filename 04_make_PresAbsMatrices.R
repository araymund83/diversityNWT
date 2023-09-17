# Load libraries ----------------------------
require(pacman)
pacman::p_load(dplyr,glue,qs, sf, stringr,tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data -----------------------------------------
data <- qs::qread('./outputs/listCommMatrices1191_1km.qs')
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[3]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
thrs <- read_csv('./inputs/prevOcc3.csv') # sp threshold values 
region <- 1:3

##extract data frames for yr1 (before) and yr2 (after)
df1 <- map(1:length(data), function(i) data[[i]][[1]])
names(df1) <- glue('{gcms}_{yr1}')
df2 <- map(1:length(data), function(i) data[[i]][[2]])
names(df2) <- glue('{gcms}_{yr2}')

species <- colnames(df1[[1]])

get_presAbs<- function(species, thrs, gcms, df){
  presAbsDF <- lapply(species, function(sp) {
    message(crayon::green('Making presAbs matrix for:', sp))
    
    thr <- filter(thrs, spec == sp)
    thr_val <- unique(thr$meandens)
    df_sp <- df %>% select(all_of(sp))
    df_sp <- ifelse(df_sp >= thr_val, 1, 0)
    
    # Check if any NA values exist in the matrix
    if (any(is.na(df_sp))) {
      message(crayon::red('WARNING: NA values detected in the matrix for', sp))
    }
    return(df_sp)
  })
  
presAbsDFMatrix <- bind_cols(presAbsDF)
# Remove species with NA values
presAbsDFMatrix<- presAbsDFMatrix %>% select_if(~ !any(is.na(.)))
out <- './tables/pres_absMat'

 ifelse(!dir.exists(out), dir.create(out), print ('Directory already exists'))
 qs::qsave (presAbsDFMatrix, glue('{out}/presAbsmatrix_INM-CM4_{yr2}.qs'))
 return(presAbsDFMatrix)
}
result <- get_presAbs(species, thrs, gcms, df2[[3]])



                   