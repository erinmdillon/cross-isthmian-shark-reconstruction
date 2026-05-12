
#######################################################################
############ AGE-UNCERTAINTY BOOTSTRAP ################################
#######################################################################

## Overview ##
#This script explores how uncertainty in the age-depth model outputs from BChron 
#affect our denticle accumulation model estimates in the Pacific and produces figure S6.


## 0. Set-up ##

# Libraries
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(parallel)
library(doSNOW)
library(progress)

#Source functions
source("functions/model_fit.R")

#Read in denticle count data
dd_pan_age_clean <- read.csv("clean/dd_pan_age_clean.csv",header=T,na.strings="NA") %>% 
  select(-c(X))

#Read in BChron posterior
pacific_posterior_raw <- read.csv("data/bchron_iter.csv",header=T)


## 1. Fitting the selected model ##

#Subsetting and renaming cleaned data for use here
pacific_data <- dd_pan_age_clean %>% 
  filter(basin=="Pacific") %>% 
  mutate(age_group = fct_relevel(age_group,c("pre_exp","modern")))

#Model
m_selected <- glmmTMB(total_dd_count ~ age_group + log(sed_weight_kg) + (1|sample_age) + offset(log(no_years)),
                      data=pacific_data,
                      family=nbinom2)
summary(m_selected)

#Just keep the columns that correspond with the denticle data
keep_samples <- pacific_data$unique
keep_samples_sanitized <- make.names(keep_samples)

#Now subset and reorder - keep just columns that correspond with a sample row (in order)
pacific_posterior <- pacific_posterior_raw[, match(keep_samples_sanitized, colnames(pacific_posterior_raw))]


## 2. Running time bootstrap in parallel. This will use your computer cores to speed up the process ##

nboot <- nrow(pacific_posterior)

cores <- detectCores() - 2 # Leaving one core for the general use of the PC.
cl <- makeCluster(cores)
registerDoSNOW(cl)

pb <- txtProgressBar(min = 1, max = nboot, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# This program fits the model structure used that the selected model uses, but it
# changes the time offset by the values in the posterior distribution produced by BChron.

time_fits <- foreach(i = 1:nboot, 
                     .options.snow = opts, 
                     .packages = c("glmmTMB")) %dopar% {
                       #Extracting offsets from posterior
                       boot_offsets <- as.numeric(pacific_posterior[i,]) + 0.001
                       
                       #Fitting model and simulating data from it
                       model_boot <- glmmTMB(total_dd_count ~ age_group + log(sed_weight_kg) + (1|sample_age) + offset(log(boot_offsets)),
                                             data=pacific_data,
                                             family=nbinom2)
                       
                       conv <- model_boot$fit$convergence
                       
                       sim_boot <- simulate_data(model_boot, nsim=1)
                       
                       #Extracting confidence intervals for the parameters
                       conf_int <- confint(model_boot)[1:3,]
                       rownames(conf_int) <- c("pre_exp", "modern", "b")
                       colnames(conf_int) <- c("ci_low", "ci_high", "estimate")
                       
                       return(list(model=model_boot, sim=sim_boot$sim, mu=sim_boot$mu, conv=conv, conf_int=conf_int))
                     }
stopCluster(cl)


## 3. Extracting model parameters from bootstrap samples ##

# This code extracts the parameters from each bootstrap model to create parameter bootstrap distributions
# and determine where do the values obtained in model_selected fall. If the fitted values
# obtained with model_selected fall inside of the bootstrap distributions, we can conclude
# that there is a level of robustness, with respect to time uncertainty
# in the parameter estimates. The plots are saved in a pdf file.

#Extracting convergence. This excludes fits that didn't converge
convergence <- sapply(time_fits, function(t) t$conv)
time_fits <- time_fits[convergence == 0]

#Extracting modern and pre_exp intercepts
modern_intercepts <- sapply(time_fits, function(t) t$conf_int["modern", "estimate"])
pre_exp_intercepts <- sapply(time_fits, function(t) t$conf_int["pre_exp", "estimate"])

#Extracting log(sed_weight_kg) slopes
b_slopes <- sapply(time_fits, function(t) t$conf_int["b", "estimate"])

#Extracting dispersion parameters
phis <- sapply(time_fits, function(t) sigma(t$model))

#Plotting distributions
pdf("figures/figS6.pdf",width=12,height=4)
par(mfrow=c(1,3))

plot(density(pre_exp_intercepts), lwd=2, main="A  Intercept (pre-exploitation)",xlab="Estimate",ylab="Density")
abline(v=fixef(m_selected)$cond[1], col="red", lwd=2, lty=2)

plot(density(modern_intercepts), lwd=2, main="B  Recent difference",xlab="Estimate",ylab="")
abline(v=fixef(m_selected)$cond[2], col="red", lwd=2, lty=2)

plot(density(b_slopes), lwd=2, main="C  Sample weight coefficient",xlab="Estimate",ylab="")
abline(v=fixef(m_selected)$cond[3], col="red", lwd=2, lty=2)

par(mfrow=c(1,1))
dev.off()


## 4. Difference significance visualization ##

# Every error bar that passes through zero is a non-significant difference from zero.
# This means that there is no difference between modern and pre_exp for that specific
# bootstrap iteration.

# This is a way of considering both time uncertainty and sampling uncertainty when deciding
# whether there is evidence for a difference between modern and pre_exp.

#Extracting low and high confidence intervals for the pre_exp difference
pre_exp_ci <- t(sapply(time_fits, 
                       function(t) t$conf_int["pre_exp", c("ci_low", "ci_high", "estimate")])) %>%
  data.frame() %>%
  mutate(pcolor = ifelse(sign(ci_low) == sign(ci_high), "significant", "non-significant"))

#Creating plot
ggplot(pre_exp_ci, aes(x=1:nrow(pre_exp_ci), y=estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high)) +
  geom_hline(yintercept=0, col="red", lty=2, linewidth=1) +
  theme_bw()

#Splitting the dataframe into smaller chunks
pre_exp_ci_split <- pre_exp_ci %>%
  mutate(n = 1:nrow(pre_exp_ci), chunk = (row_number() - 1) %/% 50 + 1) %>%
  group_split(chunk)

#Creating plots
plot_list <- list()

for(i in 1:length(pre_exp_ci_split)) {
  plot_list[[i]] <- ggplot(pre_exp_ci_split[[i]], aes(x=n, y=estimate, color=pcolor)) +
    geom_point(size=1.5) +
    geom_errorbar(aes(ymin=ci_low, ymax=ci_high), linewidth=0.8) +
    geom_hline(yintercept=0, col="red", lty=2, linewidth=1) +
    theme_bw() + 
    labs(x="Bootstrap Iteration", y="Pre-Exp Differences and 95% CI") +
    scale_color_manual(values=c("significant"="#4DAF4A", "non-significant"="#377EB8"))
}
