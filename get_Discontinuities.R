source('./R/GapRarityIndex.R')

traits <- qs::qread(glue('./inputs/traits17FGStandard_2023.qs'))
#create a  vector with names that will be used to calculate the Gap Rarity Index (GRI)
bodyMass <- traits$Masslog
names(bodyMass)<- rownames(traits)

#call
hnull<-Neutral.Null(bodyMass)
Bootstrap.gaps<-DD(bodyMass,hnull)

out <- './outputs/DD_gaps'
ifelse(!file.exists(out), dir_create(out), print('Already exists'))

write.csv(table(Bootstrap.gaps[,2]), './outputs/DD_gaps/sppGaps.csv')
#it detected 16 gaps 

