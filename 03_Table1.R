
####################
## Create Table 1 ##
####################

#written by: Emma Gause 
#Date: 11/15/23
#Last updated: "

#load libraries:
library("tidyverse")
library("dplyr")
library("tableone")
library("ggplot2")
library("ggdist")
library("gghalves")
library("ggthemes")

#Create path to directory
datadir <- "[directory path]"

#read in long data
data <- readRDS(paste0(datadir, "Analysis/Analysis_Set_111523.rds"))

##------------------------------------------------------------------------##

#create log of WMH var
data <- data %>% mutate(log_wmh = log(Total_wmh))
summary(data$log_wmh)

##------------------------------------------------------------------------##

#Create Table 1 from complete cases
str(data)

tab1r <- CreateCatTable(vars = c("study","female", "college", "inc_gt55k", 
                                 "married_partner", "social", "exc_health",
                                 "ever_smoke", "heavy_drink_fact", "daily_physical"),
                        data = data, strata = "race_fact", test = FALSE, includeNA = TRUE,
                        addOverall = TRUE)
tab1r

tab2r <- CreateContTable(vars = c("NIHTLBX_depr_theta", "age_bl_x", "Cerebrum_tcv", 
                                  "log_wmh", "Cerebrum_gray", "Cerebrum_white"),
                         data = data, strata = "race_fact",
                         funcNames = c("n", "miss", "p.miss", "mean", "sd", "median", "p25", "p75"),
                         test = FALSE, addOverall = TRUE)
summary(tab2r)


##------------------------------------------------------------------------##

#create plot showing boxplots of depressive symptom distribution by race/ethnicity

colnames(data)


tiff(filename = "[directory path]/Figures/Baseline_Depressive_Symptoms.tiff",
     res = 300, width = 2000, height = 2000, units = "px")

expP <- data %>% ggplot(aes(race_fact, NIHTLBX_depr_theta, fill = race_fact)) + 
  geom_boxplot(width = .18, outlier.shape = NA, alpha = 0.8) + 
  gghalves::geom_half_point(side = "l", range_scale = 0, shape = 95, size = 8, alpha = .3) + 
  theme_gdocs() + 
  labs(x = "Race/Ethnicity Groups",
       y = "Baseline Depressive Symptoms")

expP + theme(legend.position = "none")

dev.off()


