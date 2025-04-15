
##########################
## Frequentist Analysis ##
##########################

#written by: Emma Gause 
#Date: 11/15/23
#Last updated: "

#load libraries:
library("tidyverse")
library("dplyr")

#Create path to directory
datadir <- "[directory path]"

#read in long data
dat <- readRDS(paste0(datadir, "Analysis/Analysis_Set_111523.rds"))

##------------------------------------------------------------------------##

#take the log of WMH because highly skewed
dat <- dat %>% mutate(log_wmh = log(Total_wmh))

#CENTER CONTINUOUS VARIABLES [dep already done]
#this gives better performance for MCMC
dat <- dat %>% mutate(age_cen = age_bl_x-mean(age_bl_x),
                      c_offset_cen = Cerebrum_tcv-mean(Cerebrum_tcv),
                      log_wmh_cen = log_wmh-mean(log_wmh),
                      c_grey_cen = Cerebrum_gray-mean(Cerebrum_gray),
                      c_white_cen = Cerebrum_white-mean(Cerebrum_white),
                      time2mri_cen = time2mri-mean(time2mri))

#break into race/eth groups
asian <- dat %>% filter(race_fact=="Asian")
black <- dat %>% filter(race_fact=="Black")
latinx <- dat %>% filter(race_fact=="LatinX")
white <- dat %>% filter(race_fact=="White")

##------------------------------------------------------------------------##

#############
## OVERALL ##

## Unadjusted
wmh <- lm(data = dat, log_wmh_cen~NIHTLBX_depr_theta + c_offset_cen)

cgrey <- lm(data = dat, c_grey_cen~NIHTLBX_depr_theta + c_offset_cen)

cwhite <- lm(data = dat, c_white_cen~NIHTLBX_depr_theta + c_offset_cen)

cbind(round(coef(wmh), 2), round(confint(wmh), 2))
cbind(round(coef(cgrey), 2), round(confint(cgrey), 2))
cbind(round(coef(cwhite), 2), round(confint(cwhite), 2))


## Adjusted
wmhx <- lm(data = dat, log_wmh_cen~NIHTLBX_depr_theta + c_offset_cen +
            female+age_cen+time2mri_cen+college)

cgreyx <- lm(data = dat, c_grey_cen~NIHTLBX_depr_theta + c_offset_cen +
              female+age_cen+time2mri_cen+college)

cwhitex <- lm(data = dat, c_white_cen~NIHTLBX_depr_theta + c_offset_cen +
               female+age_cen+time2mri_cen+college)


cbind(round(coef(wmhx), 2), round(confint(wmhx), 2))
cbind(round(coef(cgreyx), 2), round(confint(cgreyx), 2))
cbind(round(coef(cwhitex), 2), round(confint(cwhitex), 2))


##------------------------------------------------------------------------##

#############
## Asian ##

## Unadjusted
awmh <- lm(data = asian, log_wmh_cen~NIHTLBX_depr_theta + c_offset_cen)

acgrey <- lm(data = asian, c_grey_cen~NIHTLBX_depr_theta + c_offset_cen)

acwhite <- lm(data = asian, c_white_cen~NIHTLBX_depr_theta + c_offset_cen)

cbind(round(coef(awmh), 2), round(confint(awmh), 2))
cbind(round(coef(acgrey), 2), round(confint(acgrey), 2))
cbind(round(coef(acwhite), 2), round(confint(acwhite), 2))



## Adjusted
awmhx <- lm(data = asian, log_wmh_cen~NIHTLBX_depr_theta + c_offset_cen +
            female+age_cen+time2mri_cen+college)

acgreyx <- lm(data = asian, c_grey_cen~NIHTLBX_depr_theta + c_offset_cen +
              female+age_cen+time2mri_cen+college)

acwhitex <- lm(data = asian, c_white_cen~NIHTLBX_depr_theta + c_offset_cen +
               female+age_cen+time2mri_cen+college)

cbind(round(coef(awmhx), 2), round(confint(awmhx), 2))
cbind(round(coef(acgreyx), 2), round(confint(acgreyx), 2))
cbind(round(coef(acwhitex), 2), round(confint(acwhitex), 2))



##------------------------------------------------------------------------##

#############
## Black ##

## Unadjusted
bwmh <- lm(data = black, log_wmh_cen~NIHTLBX_depr_theta + c_offset_cen)

bcgrey <- lm(data = black, c_grey_cen~NIHTLBX_depr_theta + c_offset_cen)

bcwhite <- lm(data = black, c_white_cen~NIHTLBX_depr_theta + c_offset_cen)

cbind(round(coef(bwmh), 2), round(confint(bwmh), 2))
cbind(round(coef(bcgrey), 2), round(confint(bcgrey), 2))
cbind(round(coef(bcwhite), 2), round(confint(bcwhite), 2))



## Adjusted
bwmhx <- lm(data = black, log_wmh_cen~NIHTLBX_depr_theta + c_offset_cen +
              female+age_cen+time2mri_cen+college)

bcgreyx <- lm(data = black, c_grey_cen~NIHTLBX_depr_theta + c_offset_cen +
                female+age_cen+time2mri_cen+college)

bcwhitex <- lm(data = black, c_white_cen~NIHTLBX_depr_theta + c_offset_cen +
                 female+age_cen+time2mri_cen+college)

cbind(round(coef(bwmhx), 2), round(confint(bwmhx), 2))
cbind(round(coef(bcgreyx), 2), round(confint(bcgreyx), 2))
cbind(round(coef(bcwhitex), 2), round(confint(bcwhitex), 2))


##------------------------------------------------------------------------##

#############
## LatinX ##

## Unadjusted
lwmh <- lm(data = latinx, log_wmh_cen~NIHTLBX_depr_theta + c_offset_cen)

lcgrey <- lm(data = latinx, c_grey_cen~NIHTLBX_depr_theta + c_offset_cen)

lcwhite <- lm(data = latinx, c_white_cen~NIHTLBX_depr_theta + c_offset_cen)

cbind(round(coef(lwmh), 2), round(confint(lwmh), 2))
cbind(round(coef(lcgrey), 2), round(confint(lcgrey), 2))
cbind(round(coef(lcwhite), 2), round(confint(lcwhite), 2))



## Adjusted
lwmhx <- lm(data = latinx, log_wmh_cen~NIHTLBX_depr_theta + c_offset_cen +
              female+age_cen+time2mri_cen+college)

lcgreyx <- lm(data = latinx, c_grey_cen~NIHTLBX_depr_theta + c_offset_cen +
                female+age_cen+time2mri_cen+college)

lcwhitex <- lm(data = latinx, c_white_cen~NIHTLBX_depr_theta + c_offset_cen +
                 female+age_cen+time2mri_cen+college)

cbind(round(coef(lwmhx), 2), round(confint(lwmhx), 2))
cbind(round(coef(lcgreyx), 2), round(confint(lcgreyx), 2))
cbind(round(coef(lcwhitex), 2), round(confint(lcwhitex), 2))



##------------------------------------------------------------------------##

#############
## White ##

## Unadjusted
wwmh <- lm(data = white, log_wmh_cen~NIHTLBX_depr_theta + c_offset_cen)

wcgrey <- lm(data = white, c_grey_cen~NIHTLBX_depr_theta + c_offset_cen)

wcwhite <- lm(data = white, c_white_cen~NIHTLBX_depr_theta + c_offset_cen)

cbind(round(coef(wwmh), 2), round(confint(wwmh), 2))
cbind(round(coef(wcgrey), 2), round(confint(wcgrey), 2))
cbind(round(coef(wcwhite), 2), round(confint(wcwhite), 2))



## Adjusted
wwmhx <- lm(data = white, log_wmh_cen~NIHTLBX_depr_theta + c_offset_cen +
              female+age_cen+time2mri_cen+college)

wcgreyx <- lm(data = white, c_grey_cen~NIHTLBX_depr_theta + c_offset_cen +
                female+age_cen+time2mri_cen+college)

wcwhitex <- lm(data = white, c_white_cen~NIHTLBX_depr_theta + c_offset_cen +
                 female+age_cen+time2mri_cen+college)

cbind(round(coef(wwmhx), 2), round(confint(wwmhx), 2))
cbind(round(coef(wcgreyx), 2), round(confint(wcgreyx), 2))
cbind(round(coef(wcwhitex), 2), round(confint(wcwhitex), 2))
