
##################################
## Bayesian Analysis - Adjusted ##
## Grey matter volume           ##
##################################

#written by: Emma Gause 
#Date: 11/15/23
#Last updated: 03/13/24

#load libraries:
library("tidyverse")
library("dplyr")
library("nimble")
library("posterior")
library("bayesplot")
library("ggplot2")
library("ggthemes")


#Create path to directory
datadir <- "[directory path]"

#read in long data
dat <- readRDS(paste0(datadir, "Analysis/Analysis_Set_111523.rds"))


tic <- Sys.time()

##------------------------------------------------------------------------##

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


#Constants
#Set initial values for MCMC samples (very important for convergence)
inits <- list(alpha = 0, beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0, sigma = 1)

################################
## OVERALL ##

no <- nrow(dat) # number of observations

####
#IMPRECISE

ocode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-1.5, sd = 4) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components ... edited due to non-convergence
  
  
  # regression formula
  for (i in 1:no) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = no)
odata <- list(y = dat$c_grey_cen, x1 = dat$NIHTLBX_depr_theta, x2 = dat$c_offset_cen,
              x3 = dat$female, x4 = dat$age_cen, x5 = dat$time2mri_cen, x6 = dat$college)


#model inputs
omodel <- nimbleModel(ocode, constants = constants, data = odata, inits = inits)
mcmcConf_o <- configureMCMC(omodel)

#Run the MCMC simulations.
overall_mcmc <- nimbleMCMC(
  code = ocode,
  data = odata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)

####
#SHIFTED LOWER

ocode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-2.5, sd = 2) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:no) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = no)
odata <- list(y = dat$c_grey_cen, x1 = dat$NIHTLBX_depr_theta, x2 = dat$c_offset_cen,
              x3 = dat$female, x4 = dat$age_cen, x5 = dat$time2mri_cen, x6 = dat$college)


#model inputs
omodel <- nimbleModel(ocode, constants = constants, data = odata, inits = inits)
mcmcConf_o <- configureMCMC(omodel)

#Run the MCMC simulations.
overall_low <- nimbleMCMC(
  code = ocode,
  data = odata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)

####
#SHIFTED HIGHER

ocode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-0.5, sd = 2) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:no) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = no)
odata <- list(y = dat$c_grey_cen, x1 = dat$NIHTLBX_depr_theta, x2 = dat$c_offset_cen,
              x3 = dat$female, x4 = dat$age_cen, x5 = dat$time2mri_cen, x6 = dat$college)


#model inputs
omodel <- nimbleModel(ocode, constants = constants, data = odata, inits = inits)
mcmcConf_o <- configureMCMC(omodel)

#Run the MCMC simulations.
overall_high <- nimbleMCMC(
  code = ocode,
  data = odata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)




################################
## Asian ##

na <- nrow(asian) # number of observations

####
#IMPRECISE

acode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-1.5, sd = 4) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:na) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = na)
adata <- list(y = asian$c_grey_cen, x1 = asian$NIHTLBX_depr_theta, x2 = asian$c_offset_cen,
              x3 = asian$female, x4 = asian$age_cen, x5 = asian$time2mri_cen, x6 = asian$college)

#model inputs
amodel <- nimbleModel(acode, constants = constants, data = adata, inits = inits)
mcmcConf_a <- configureMCMC(amodel)

#Run the MCMC simulations.
asian_mcmc <- nimbleMCMC(
  code = acode,
  data = adata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)

####
#ISHIFTED LOWER

acode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-2.5, sd = 2) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:na) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = na)
adata <- list(y = asian$c_grey_cen, x1 = asian$NIHTLBX_depr_theta, x2 = asian$c_offset_cen,
              x3 = asian$female, x4 = asian$age_cen, x5 = asian$time2mri_cen, x6 = asian$college)

#model inputs
amodel <- nimbleModel(acode, constants = constants, data = adata, inits = inits)
mcmcConf_a <- configureMCMC(amodel)

#Run the MCMC simulations.
asian_low <- nimbleMCMC(
  code = acode,
  data = adata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)

####
#SHIFTED HIGHER

acode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-0.5, sd = 2) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:na) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = na)
adata <- list(y = asian$c_grey_cen, x1 = asian$NIHTLBX_depr_theta, x2 = asian$c_offset_cen,
              x3 = asian$female, x4 = asian$age_cen, x5 = asian$time2mri_cen, x6 = asian$college)

#model inputs
amodel <- nimbleModel(acode, constants = constants, data = adata, inits = inits)
mcmcConf_a <- configureMCMC(amodel)

#Run the MCMC simulations.
asian_high <- nimbleMCMC(
  code = acode,
  data = adata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)



################################
## Black ##

nb <- nrow(black) # number of observations

####
#IMPRECISE

bcode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-1.5, sd = 4) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:nb) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = nb)
bdata <- list(y = black$c_grey_cen, x1 = black$NIHTLBX_depr_theta, x2 = black$c_offset_cen,
              x3 = black$female, x4 = black$age_cen, x5 = black$time2mri_cen, x6 = black$college)


#model inputs
bmodel <- nimbleModel(bcode, constants = constants, data = bdata, inits = inits)
mcmcConf_b <- configureMCMC(bmodel)

#Run the MCMC simulations.
black_mcmc <- nimbleMCMC(
  code = bcode,
  data = bdata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)

####
#SHIFTED LOWER

bcode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-2.5, sd = 2) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:nb) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = nb)
bdata <- list(y = black$c_grey_cen, x1 = black$NIHTLBX_depr_theta, x2 = black$c_offset_cen,
              x3 = black$female, x4 = black$age_cen, x5 = black$time2mri_cen, x6 = black$college)


#model inputs
bmodel <- nimbleModel(bcode, constants = constants, data = bdata, inits = inits)
mcmcConf_b <- configureMCMC(bmodel)

#Run the MCMC simulations.
black_low <- nimbleMCMC(
  code = bcode,
  data = bdata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)

####
#SHIFTED HIGHER

bcode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-0.5, sd = 2) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:nb) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = nb)
bdata <- list(y = black$c_grey_cen, x1 = black$NIHTLBX_depr_theta, x2 = black$c_offset_cen,
              x3 = black$female, x4 = black$age_cen, x5 = black$time2mri_cen, x6 = black$college)


#model inputs
bmodel <- nimbleModel(bcode, constants = constants, data = bdata, inits = inits)
mcmcConf_b <- configureMCMC(bmodel)

#Run the MCMC simulations.
black_high <- nimbleMCMC(
  code = bcode,
  data = bdata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)


################################
## LatinX ##

nl <- nrow(latinx) # number of observations

####
#IMPRECISE

lcode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-1.5, sd = 4) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:nl) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = n)
ldata <- list(y = latinx$c_grey_cen, x1 = latinx$NIHTLBX_depr_theta, x2 = latinx$c_offset_cen,
              x3 = latinx$female, x4 = latinx$age_cen, x5 = latinx$time2mri_cen, x6 = latinx$college)


#model inputs
lmodel <- nimbleModel(lcode, constants = constants, data = ldata, inits = inits)
mcmcConf_l <- configureMCMC(lmodel)

#Run the MCMC simulations.
latinx_mcmc <- nimbleMCMC(
  code = lcode,
  data = ldata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)

####
#SHIFTED LOWER

lcode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-2.5, sd = 2) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:nl) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = n)
ldata <- list(y = latinx$c_grey_cen, x1 = latinx$NIHTLBX_depr_theta, x2 = latinx$c_offset_cen,
              x3 = latinx$female, x4 = latinx$age_cen, x5 = latinx$time2mri_cen, x6 = latinx$college)


#model inputs
lmodel <- nimbleModel(lcode, constants = constants, data = ldata, inits = inits)
mcmcConf_l <- configureMCMC(lmodel)

#Run the MCMC simulations.
latinx_low <- nimbleMCMC(
  code = lcode,
  data = ldata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)


####
#SHIFTED HIGHER

lcode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-0.5, sd = 2) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:nl) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = n)
ldata <- list(y = latinx$c_grey_cen, x1 = latinx$NIHTLBX_depr_theta, x2 = latinx$c_offset_cen,
              x3 = latinx$female, x4 = latinx$age_cen, x5 = latinx$time2mri_cen, x6 = latinx$college)


#model inputs
lmodel <- nimbleModel(lcode, constants = constants, data = ldata, inits = inits)
mcmcConf_l <- configureMCMC(lmodel)

#Run the MCMC simulations.
latinx_high <- nimbleMCMC(
  code = lcode,
  data = ldata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)




################################
## White ##

nw <- nrow(white) # number of observations

####
#IMPRECISE

wcode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-1.5, sd = 4) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:nw) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = n)
wdata <- list(y = white$c_grey_cen, x1 = white$NIHTLBX_depr_theta, x2 = white$c_offset_cen,
              x3 = white$female, x4 = white$age_cen, x5 = white$time2mri_cen, x6 = white$college)


#model inputs
wmodel <- nimbleModel(wcode, constants = constants, data = wdata, inits = inits)
mcmcConf_w <- configureMCMC(wmodel)

#Run the MCMC simulations.
white_mcmc <- nimbleMCMC(
  code = wcode,
  data = wdata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)

####
#SHIFTED LOWER

wcode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-2.5, sd = 2) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:nw) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = n)
wdata <- list(y = white$c_grey_cen, x1 = white$NIHTLBX_depr_theta, x2 = white$c_offset_cen,
              x3 = white$female, x4 = white$age_cen, x5 = white$time2mri_cen, x6 = white$college)


#model inputs
wmodel <- nimbleModel(wcode, constants = constants, data = wdata, inits = inits)
mcmcConf_w <- configureMCMC(wmodel)

#Run the MCMC simulations.
white_low <- nimbleMCMC(
  code = wcode,
  data = wdata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)

####
#SHIFTED HIGHER

wcode <- nimbleCode({
  # priors for parameters
  alpha ~ dnorm(0, sd = 1) # prior for alpha
  beta1 ~ dnorm(-0.5, sd = 2) # prior for depression
  beta2 ~ dnorm(0.1, sd = 1) # prior for beta2 // ICV offset [slightly larger volumes]
  beta3 ~ dnorm(0, sd = 1) # prior for beta3 // gender [no hypothesis for female gender?]
  beta4 ~ dnorm(-1, sd = 2) # prior for beta4 // age [higher age, smaller volume]
  beta5 ~ dnorm(-0.1, sd = 1) # prior for beta5 // time2mri [higher time, smaller volume, but small]
  beta6 ~ dnorm(1, sd = 2) # prior for beta6 //college [more college, more volume]
  sigma ~ dnorm(1, 30) # prior for variance components  edited due to non-convergence
  
  
  # regression formula
  for (i in 1:nw) {
    # n is the number of observations we have in the data
    mu[i] <- alpha+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]+beta4*x4[i]+beta5*x5[i]+beta6*x6[i] # manual entry of linear predictors
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

#Final preparation of data into lists for input into `NIMBLE` MCMC running
constants <- list(n = n)
wdata <- list(y = white$c_grey_cen, x1 = white$NIHTLBX_depr_theta, x2 = white$c_offset_cen,
              x3 = white$female, x4 = white$age_cen, x5 = white$time2mri_cen, x6 = white$college)


#model inputs
wmodel <- nimbleModel(wcode, constants = constants, data = wdata, inits = inits)
mcmcConf_w <- configureMCMC(wmodel)

#Run the MCMC simulations.
white_high <- nimbleMCMC(
  code = wcode,
  data = wdata,
  constants = constants,
  inits = inits,
  niter = 10000, # run 10000 samples 
  nburnin = 1000, # burn in for 1000 iterations (10% is about right)
  setSeed = 111523,
  samplesAsCodaMCMC = TRUE
)



##------------------------------------------------------------------------##

#Plot Them! 

o <- as_draws(overall_mcmc)
oo <- subset_draws(o, c("beta1"))

ol <- as_draws(overall_low)
ool <- subset_draws(ol, c("beta1"))

oh <- as_draws(overall_high)
ooh <- subset_draws(oh, c("beta1"))


a <- as_draws(asian_mcmc)
aa <- subset_draws(a, c("beta1"))

al <- as_draws(asian_low)
aal <- subset_draws(al, c("beta1"))

ah <- as_draws(asian_high)
aah <- subset_draws(ah, c("beta1"))


b <- as_draws(black_mcmc)
bb <- subset_draws(b, c("beta1"))

bl <- as_draws(black_low)
bbl <- subset_draws(bl, c("beta1"))

bh <- as_draws(black_high)
bbh <- subset_draws(bh, c("beta1"))


l <- as_draws(latinx_mcmc)
ll <- subset_draws(l, c("beta1"))

l_l <- as_draws(latinx_low)
ll_l <- subset_draws(l_l, c("beta1"))

lh <- as_draws(latinx_high)
llh <- subset_draws(lh, c("beta1"))


w <- as_draws(white_mcmc)
ww <- subset_draws(w, c("beta1"))

wl <- as_draws(white_low)
wwl <- subset_draws(wl, c("beta1"))

wh <- as_draws(white_high)
wwh <- subset_draws(wh, c("beta1"))


osum <- c(mean(oo), sd(oo))
osuml <- c(mean(ool), sd(ool))
osumh <- c(mean(ooh), sd(ooh))

asum <- c(mean(aa), sd(aa))
asuml <- c(mean(aal), sd(aal))
asumh <- c(mean(aah), sd(aah))

bsum <- c(mean(bb), sd(bb))
bsuml <- c(mean(bbl), sd(bbl))
bsumh <- c(mean(bbh), sd(bbh))

lsum <- c(mean(ll), sd(ll))
lsuml <- c(mean(ll_l), sd(ll_l))
lsumh <- c(mean(llh), sd(llh))

wsum <- c(mean(ww), sd(ww))
wsuml <- c(mean(wwl), sd(wwl))
wsumh <- c(mean(wwh), sd(wwh))


sums <- rbind(osum, osuml, osumh,
              asum, asuml, asumh,
              bsum, bsuml, bsumh,
              lsum, lsuml, lsumh,
              wsum, wsuml, wsumh)

colnames(sums) <- c("mu", "sd")

strata <- c("overall", "overall", "overall",
                     "asian", "asian", "asian",
                     "black", "black", "black",
                     "latinX", "latinX", "latinX",
                     "white", "white", "white")
type <- c("imprecise", "shifted lower", "shifted higher",
          "imprecise", "shifted lower", "shifted higher",
          "imprecise", "shifted lower", "shifted higher",
          "imprecise", "shifted lower", "shifted higher",
          "imprecise", "shifted lower", "shifted higher")


sumsx <- cbind(strata, type, sums)
sumsx <- as.data.frame(sumsx)

write_csv(sumsx, paste0(datadir, "Results/Grey_Matter_Sensitivity_Posteriors.csv"))

##------------------------------------------------------------------------##

toc <- Sys.time()
toc - tic
