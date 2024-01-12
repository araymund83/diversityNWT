# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(dplyr,glue, purrr, qs, RColorBrewer, stringr, sf, terra, tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
#Load lat and long for each pixel
pixel_yx <- qs::qread('./inputs/xyPixelID.qs')
path <- './outputs/tbi_1km'
files <- list.files(path, pattern = '72spp', full.names = TRUE)
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[2]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")


limt <- sf::st_read('inputs/NT1_BCR6/NT1_BCR6_poly.shp') 
ecrg <- sf::st_read('inputs/ecoregions/NWT_ecoregions_dissolvec.shp')

# Extract by mask for the ecoregions ---------------------------------------
plot(st_geometry(ecrg))
limt <- sf::st_transform(x = limt, crs = targetCRS)
ecrg <- sf::st_transform(x = ecrg, crs = targetCRS)
plot(st_geometry(ecrg))
plot(st_geometry(limt))

##read the TBI outputs 
read_TBI <- function(gcms, yr2, files){
  map(gcms, ~ {
    gc <- .
    message(crayon::green('Loading TBI results for GCM:', gc))
    fls <- grep(gc, files, value = TRUE)
    fl1 <- grep(yr2, fls, value = TRUE)
    df1 <- qs::qread(fl1)
    return(df1)
  })
}

# Apply the function ------------------------------------------------------
data <- read_TBI(gcms, yr2, files)

#function to extract the BCD.mat element for each gcm
# get_elem <- function(x, elem){
#   if(!is.list(x[[1]])) x[[elem]]
#   else lapply(x, get_elem, elem)
# }

# I had to create the list because the function read_TBI only reads the first element of 
#each list. 


CanESM2<- qread('./outputs/tbi_1km/tbi72spp_CanESM2_2011_20312_1km.qs')
CCSM4<- qread('./outputs/tbi_1km/tbi72spp_CCSM4_2011_20312_1km.qs')
INMCM4<- qread('./outputs/tbi_1km/tbi72spp_INM-CM4_2011_20312_1km.qs')


data <- list(CanESM2_91, CCSM4_91, INMCM4_91)
#data <- list(CanESM2, CCSM4, INMCM4)
make_TBI_Ras <- function(data, gcms, yr = yr2, var){
  #extract the BCD matrix from the TBI results to work with ggplot
  p_TBI <- lapply(data,function(x) x[[1]][[1]][["p.TBI"]]) 
  valsTBI <-  lapply(data,function(x) x[[1]][[1]][["TBI"]]) 
  names(p_TBI) <- gcms
  names(valsTBI) <- gcms
  
#convert lists to df and merge with pixels xy
  p_TBIDF <- data.frame(do.call(cbind, p_TBI))
  tbiDF <- data.frame(do.call(cbind, valsTBI))
  
  p_TBIDF <- p_TBIDF %>%  mutate(pixelID = glue('S_{1:nrow(.)}'))
  tbiDF <- tbiDF %>% mutate(pixelID = glue('S_{1:nrow(.)}'))
 
  #join xy coordinates to tbi and pvalue dataframe
  p_TBIDF <- left_join(p_TBIDF, pixel_yx, by = 'pixelID')
  tbiDF <- left_join(tbiDF, pixel_yx, by = 'pixelID')
  
  tbiDF <- tbiDF %>% pivot_longer(cols = 1:3, names_to = 'gcm', values_to = 'tbi')
  p_TBIDF <- p_TBIDF %>% pivot_longer(cols = 1:3, names_to = 'gcm', values_to = 'p_val')
  
  rasDF <- left_join(tbiDF, p_TBIDF, by = c("pixelID", 'x', 'y', 'gcm'))
  sig_cutOff <- 0.05
  rasDF <- rasDF %>% mutate(yr = yr, pValue2 = case_when(p_val >= sig_cutOff ~ 'no_sig',
                        p_val < sig_cutOff ~ 'significant') ) %>%
    rename(Longitude = x, Latitude = y)
  
  color_2vals<- c('#b3b3b3', '#018571')
  var = var

 ggVar <-  ggplot() +
  geom_tile(data = rasDF, aes(x = Longitude, y = Latitude, fill = pValue2)) +
    scale_fill_manual(values = color_2vals, 
                     labels = c('No significant', 'Significant')) +
    geom_sf(data = limt, fill = NA, col = '#36454F') +
    geom_sf(data = ecrg,fill = NA,col = '#36454F') +
    ggtitle(label = glue('TBI significance'),
            subtitle = glue('{yr1}-{yr2}')) +
    theme_bw() +
    theme(legend.position = 'bottom',legend.key.width = unit(2, 'line'),
          plot.title = element_text(size = 16, face = 'bold', hjust = 0, vjust = 0.7),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12, face = 'bold'),
          strip.text = element_text(size = 14)) +
    labs(x = 'Longitude', y = 'Latitude', fill = var) +
    coord_sf() +
    facet_wrap(~ gcm)
  
  out <-  glue('./figs/BC_Ras_72spp/1km')
  ifelse(!dir.exists(out), dir.create(out, recursive = TRUE),
         print('Folder already exists'))
  ggsave(plot = ggVar, filename = glue('{out}/{var}_{yr1}_{yr2}.png'), 
         units = 'in', width = 10, height = 7, dpi = 300)
}

# Apply the function ------------------------------------------------------
 make_TBI_Ras(data = data, gcms = gcms, var = 'pValue2' )

## This code generates the plot for discrete variables 
 
 
ggTBI <-  ggplot() +
   geom_tile(data = rasDF, aes(x = Longitude, y = Latitude, fill =tbi)) +
   scale_fill_gradientn(colours = brewer.pal(n = 10, name = 'BrBG'), limit = c(0,1)) +
   geom_sf(data = limt, fill = NA, col = '#36454F') +
   geom_sf(data = ecrg,fill = NA,col = '#36454F') +
   ggtitle(label = glue('TBI index'),
           subtitle = glue('{yr1}-{yr2}')) +
   theme_bw() +
   theme(legend.position = 'bottom',legend.key.width = unit(2, 'line'),
         plot.title = element_text(size = 16, face = 'bold', hjust = 0, vjust = 0.7),
         plot.subtitle = element_text(size = 14),
         axis.title = element_text(size = 14),
         axis.text.x = element_text(size = 12),
         axis.text.y = element_text(size = 12),
         legend.text = element_text(size = 11),
         legend.title = element_text(size = 12, face = 'bold'),
         strip.text = element_text(size = 14)) +
   labs(x = 'Longitude', y = 'Latitude', fill = 'TBI') +
   coord_sf(lims_method = 'geometry_bbox') +
   facet_wrap(~ gcm)
 
 out <-  glue('./figs/BC_Ras_72spp/1km')
 ifelse(!dir.exists(out), dir.create(out, recursive = TRUE),
        print('Folder already exists'))
 ggsave(plot = ggTBI, filename = glue('{out}/TBI_{yr1}_{yr2}_.png'), 
        units = 'in', width = 10, height = 7, dpi = 300)
 

 
                   