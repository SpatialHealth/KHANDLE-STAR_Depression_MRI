
####################################################################
## Add MRI Measures to the Cohort Data Cleaned for SENAS Analysis ##
####################################################################

#written by: Emma Gause 
#Date: 07/19/23
#Last updated: 11/08/23

#load libraries:
library("tidyverse")
library("dplyr")

#Create path to directory
datadir1 <- "[directory path]"
datadir2 <- "[directory path]"
dataexp <- "[directory path]"

#read in cohort data with covariates
dat <- readRDS(paste0(datadir1, "Analysis_Data_081823.rds"))

#read in mri data
kmri <- readRDS(paste0(datadir2, "KHANDLE/Data/Imaging/khandle_T1_analysis_052322_age.rds"))
smri <- readRDS(paste0(datadir2, "STAR/Data/Imaging/k-star_T1_analysis_052322_age.rds"))

##------------------------------------------------------------------------##

str(dat)
str(kmri)
str(smri)

#create the study IDs to be compatible
head(dat$id)
head(dat$STUDYID)
head(kmri$StudyID)
tail(smri$StudyID)

#pad study IDs for cohort data to match MRI

#[redacted to preserve ID anonymity]

##------------------------------------------------------------------------##

#rbind mri measures and then merge to cohort data
colnames(kmri)
colnames(smri)
mri <- rbind(kmri, smri)

#prepare cohort data to merge to MRI set
#We want to compare m:1 so we can assess the difference in time between depression and MRI scan
str(dat)
str(mri)
  # XXXXX and XXXXX exist in MRI but not in cohort data...
  # these are Native American or missing race exclusions - they will be removed in inner_join

data <- inner_join(dat, mri, by = "StudyID", relationship = "many-to-one")

#see if they all merged --> race_fact should have no missingness 
table(data$race_fact, useNA = "ifany")
#looks good!

#how many unique IDs do we have?
ids <- data %>% select(StudyID) %>% unique()
  #560 - this is what we expect

saveRDS(data, paste0(dataexp, "MRI_Long_110823.rds"))

##------------------------------------------------------------------------##





