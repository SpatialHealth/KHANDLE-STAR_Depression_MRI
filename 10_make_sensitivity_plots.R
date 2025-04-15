##############################
## Create Sensitivity Plots ##
## All Outcomes             ##
##############################

#Edited by: Emma Gause 
#Adapted from code by: Sarah Ackley
#Date: 03/27/24
#Last updated: "

#load libraries:
library("tidyverse")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("ggpubr")


#Create path to directory
datadir <- "[directory path]"

#read in plot data
dat <- read_csv(paste0(datadir, "Results/Prior_Sensitivity.csv"))

##------------------------------------------------------------------------##

str(dat)

#we'll do each outcome separately since the x-axis limits may differ
table(dat$outcome)
wmh <- dat %>% filter(outcome == "Log(WMH)")
grey <- dat %>% filter(outcome == "Grey Matter")
white <- dat %>% filter(outcome == "White Matter")

#set up functions
trail.round <- function(x,round.digits=2){
  sprintf(paste("%.",round.digits,"f",sep=""), 
          round(x,round.digits))
}


#for loop to create plots for each sensitivity analysis and each strata

#start with wmh


wmh.list <- list()
grey.list <- list()
white.list <- list()

for(i in 1:nrow(wmh)){
  
  #first get all the correct mus and sigmas and the function to create their normal distributions
  prior.mu=wmh$prior_mu[[i]]
  prior.sigma=wmh$prior_sigma[[i]]
  prior.fun <- function(x)dnorm(x,prior.mu,prior.sigma)
  post.mu=wmh$post_mu[[i]]
  post.sigma=wmh$post_sigma[[i]]
  post.fun <- function(x)dnorm(x,post.mu,post.sigma)
  
  #get 95% credible interval values
  l95 <- qnorm(0.025, mean = post.mu, sd = post.sigma)
  u95 <- qnorm(0.975, mean = post.mu, sd = post.sigma)
  
  #first set up the plot titles
  sens.title <- bquote(atop(.("Prior:")~
                        mu == .(paste(trail.round(prior.mu,2),",",sep="")) ~ 
                        sigma == .(trail.round(prior.sigma,2)),
                        paste("Estimate = ", .(trail.round(post.mu,2)), 
                              " (", .(trail.round(l95,2)), ", ", .(trail.round(u95,2)), ")")))
  
                      
  #now create the plots
  plot <- ggplot(data.frame(x = c(-2, 3)), #specify the x-axis limits here
                 aes(x = x)) +
    stat_function(fun = prior.fun, 
                  geom = "area", 
                  fill = "#69b3a2", 
                  alpha = .5)+
    stat_function(fun = post.fun, 
                  geom = "area", 
                  fill = "#404080", 
                  alpha = .5)+
    stat_function(fun = prior.fun,color='lightgrey')+
    stat_function(fun = post.fun,color='darkgrey')+
    xlab("")+
    ylab("")+
    theme_tufte()+
    ylim(c(0,5.5))+
    ggtitle(sens.title)+
    theme(plot.title = element_text(size = 10, face = "bold", vjust = -20, hjust = 0.05))+
    theme(plot.margin = margin(c(-30,-5,-5,-5)))
  
  wmh.list[[i]] <- plot
  
  #save each plot individually 
  ggsave(plot,
         file=paste0(datadir, "Figures/Sensitivity/",wmh$strata[i], "_", wmh$outcome[i], "_", wmh$sensitivity[i], ".jpeg"),
         width=4,
         height=3)
  
}



#Next grey matter
for(i in 1:nrow(grey)){
  
  #first get all the correct mus and sigmas and the function to create their normal distributions
  prior.mu=grey$prior_mu[[i]]
  prior.sigma=grey$prior_sigma[[i]]
  prior.fun <- function(x)dnorm(x,prior.mu,prior.sigma)
  post.mu=grey$post_mu[[i]]
  post.sigma=grey$post_sigma[[i]]
  post.fun <- function(x)dnorm(x,post.mu,post.sigma)
  
  #get 95% credible interval values
  l95 <- qnorm(0.025, mean = post.mu, sd = post.sigma)
  u95 <- qnorm(0.975, mean = post.mu, sd = post.sigma)
  
  #first set up the plot titles
  sens.title <- bquote(atop(.("Prior:")~
                             mu == .(paste(trail.round(prior.mu,2),",",sep="")) ~ 
                             sigma == .(trail.round(prior.sigma,2)),
                           paste("Estimate = ", .(trail.round(post.mu,2)), 
                                 " (", .(trail.round(l95,2)), ", ", .(trail.round(u95,2)), ")")))
  
  
  #now create the plots
  plot <- ggplot(data.frame(x = c(-6, 5)), #specify the x-axis limits here
                 aes(x = x)) +
    stat_function(fun = prior.fun, 
                  geom = "area", 
                  fill = "#69b3a2", 
                  alpha = .5)+
    stat_function(fun = post.fun, 
                  geom = "area", 
                  fill = "#404080", 
                  alpha = .5)+
    stat_function(fun = prior.fun,color='lightgrey')+
    stat_function(fun = post.fun,color='darkgrey')+
    xlab("")+
    ylab("")+
    theme_tufte()+
    ylim(c(0,1.1))+
    ggtitle(sens.title)+
    theme(plot.title = element_text(size = 10, face = "bold", vjust = -20, hjust = 0.05))+
    theme(plot.margin = margin(c(-30,-5,-5,-5)))
  
  #save each plot individually
  ggsave(plot,
         file=paste0(datadir, "Figures/Sensitivity/",grey$strata[i], "_", grey$outcome[i], "_", grey$sensitivity[i], ".jpeg"),
         width=4,
         height=3)
  
  grey.list[[i]] <- plot
  
}


#Finally white matter
for(i in 1:nrow(white)){
  
  #first get all the correct mus and sigmas and the function to create their normal distributions
  prior.mu=white$prior_mu[[i]]
  prior.sigma=white$prior_sigma[[i]]
  prior.fun <- function(x)dnorm(x,prior.mu,prior.sigma)
  post.mu=white$post_mu[[i]]
  post.sigma=white$post_sigma[[i]]
  post.fun <- function(x)dnorm(x,post.mu,post.sigma)
  
  #get 95% credible interval values
  l95 <- qnorm(0.025, mean = post.mu, sd = post.sigma)
  u95 <- qnorm(0.975, mean = post.mu, sd = post.sigma)
  
  #first set up the plot titles
  sens.title <- bquote(atop(.("Prior:")~
                             mu == .(paste(trail.round(prior.mu,2),",",sep="")) ~ 
                             sigma == .(trail.round(prior.sigma,2)),
                           paste("Estimate = ", .(trail.round(post.mu,2)), 
                                 " (", .(trail.round(l95,2)), ", ", .(trail.round(u95,2)), ")")))
  
  
  #now create the plots
  plot <- ggplot(data.frame(x = c(-6, 5)), #specify the x-axis limits here
                 aes(x = x)) +
    stat_function(fun = prior.fun, 
                  geom = "area", 
                  fill = "#69b3a2", 
                  alpha = .5)+
    stat_function(fun = post.fun, 
                  geom = "area", 
                  fill = "#404080", 
                  alpha = .5)+
    stat_function(fun = prior.fun,color='lightgrey')+
    stat_function(fun = post.fun,color='darkgrey')+
    xlab("")+
    ylab("")+
    theme_tufte()+
    ylim(c(0,0.9))+
    ggtitle(sens.title)+
    theme(plot.title = element_text(size = 10, face = "bold", vjust = -20, hjust = 0.05))+
    theme(plot.margin = margin(c(-30,-5,-5,-5)))
  
  #save each plot individually 
  ggsave(plot,
         file=paste0(datadir, "Figures/Sensitivity/",white$strata[i], "_", white$outcome[i], "_", white$sensitivity[i], ".jpeg"),
         width=4,
         height=3)
  
  white.list[[i]] <- plot
  
}

#create one legend
dat$legend <- if_else(dat$prior_mu<1, 1, 0)
legp <- dat %>% ggplot(aes(x = post_mu, fill = as.factor(legend))) + 
  geom_histogram(alpha = 0.8) + 
  scale_fill_manual(values = c("#69b3a2", "#404080"), 
                    labels = c("Prior", "Posterior")) +
  labs(fill="") + 
  theme(legend.direction ="horizontal") + 
  theme(legend.position = "bottom") + 
  theme_hc()
legp
  
#ggsave(legp, file = paste0(datadir, "Figures/Sensitivity_LEGEND_ONLY.tiff"),
#       width=4, height=3)


  


##------------------------------------------------------------------------##


