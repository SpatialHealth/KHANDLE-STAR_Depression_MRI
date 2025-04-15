
###########################
## Prepare Analysis Data ##
###########################

#written by: Emma Gause 
#Date: 11/15/23
#Last updated: "

#load libraries:
library("tidyverse")
library("dplyr")

#Create path to directory
datadir <- "[directory path]"

#read in long data
long <- readRDS(paste0(datadir, "Analysis/MRI_Long_110823.rds"))

##------------------------------------------------------------------------##

str(long)

#keep only wave 1 
wave1 <- long %>% filter(waveseq==1)

#remove those without depression information
keep <- wave1 %>% filter(!is.na(NIHTLBX_depr_theta))

#create variable for time from baseline to MRI
keep <- keep %>% mutate(time2mri = Age_at_Date - age_bl_x)
summary(keep$time2mri)


#remove those missing any outcome:
  #Total_wmh
  #Cerebrum_gray
  #Cerebrum_white

dat <- keep %>% filter(!is.na(Total_wmh)&!is.na(Cerebrum_gray)&!is.na(Cerebrum_white))

#look into complete cases
#adjust for:
  #female
  #age_bl_x
  #Cerebrum_tcv
  #time2mri
  #college OR inc_gt55k

#income and education have same function in DAG, which to use?
table(dat$college, dat$inc_gt55k, useNA = "ifany", deparse.level = 2) 
#college is not missing -- use this

#get complete cases
cc <- dat %>% filter(!is.na(female)&!is.na(age_bl_x)&!is.na(Cerebrum_tcv)&!is.na(time2mri)&!is.na(college))
#No additional missingness!!!! yay!


saveRDS(dat, paste0(datadir, "Analysis/Analysis_Set_111523.rds"))












