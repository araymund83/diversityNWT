## Loading R packages
library(pacman)
p_load( ade4, ape, BAT, gawdis, ggplot2, ggrepel, glue, hypervolume, tidyverse,vegan)

##############################

# Load data ---------------------------------------------------------------
path <- './inputs'
traits <- read.csv(glue('{path}/birdSpecies_traits2023.csv'))
data <- qs::qread('./outputs/listCommMatrices_1131_1km.qs')
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[2]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')

splu <- traits %>% select(species, Habitat) %>% 
  mutate(spID = as.character(1:nrow(traits)))

#convert character traits to categorical traits
traits <- traits %>% mutate_at(c(7:14), ~scale(.) %>% as.vector) %>% # normalize numerical traits 
  mutate(Masslog10 = log10(Mass)) %>% 
  mutate(across(c(2:6,16), as.factor)) %>% 
  column_to_rownames('species')
traits2 <- traits %>% select(2:6,7:13,16)
write.csv(traits2, glue('{path}/birdTraitsStandard_2023.csv'))

traitsFG<- traits2 %>%
           mutate(FG = glue('{Feeding}-{Forage_Breed1}')) %>% 
           select(FG, Masslog10) 
qs::qsave(traits2, './inputs/traitsMassFG_2023.qs')
