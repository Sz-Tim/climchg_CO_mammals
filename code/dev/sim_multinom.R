# Multinomial approach
# Simulated data
# Tim Szewcyzk


# This script simulates data in the form of the historical records. It assumes
# that pr(detect individual) is known for each species, calculated from the
# mark/recapture data.
#
# Then, it fits a JAGS model where the number of individuals of each species
# detected at each elevation is drawn from a multinomial distribution, where the
# probability for each species is weighted by: 
#   1) the overall number of detections
#   2) the individual-level detection probability
#   3) the number of elevational bins away from its nearest detection
#   4) the patchiness of its detections within its interpolated range.
#
# Currently, we're assuming that the species are within a taxon (e.g., voles),
# but I'm not sure that's strictly necessary.
#
# The variables are: 
# n.spp: number of species
# n.el: number of elevational bins
# delta: individual-level detection probability, length=n.spp
# b: intercept + slopes for regression of 3) & 4)
# spAbund: total detections along gradient, length=n.spp
# Z: true interpolated range, nrows=n.el, ncol=n.spp
# N: true abundances, nrows=n.el, ncol=n.spp
# Y: detections, nrows=n.el, ncol=n.spp
# Ytot: total detections at each elevation, length=n.el
# spp_pool: sampling pool = N * delta, nrows=n.el, ncol=n.spp
# binsAway: number of bins away detected range, nrows=n.el, ncol=n.spp
# interpPatchy: proportion non-detections in interpolated range, length=n.spp
# rng: true range boundaries, nrow=2, ncol=n.spp




# set up
# if needed: install.packages(c("tidyverse", "rjags", "ggmcmc"))
library(tidyverse); library(rjags); library(ggmcmc)


########
## Data generation
########

# Set parameters
n.spp <- 9
n.el <- 33
delta <- runif(n.spp, min=0, max=0.5)
b <- c(0, -1, 1)

# initialize storage objects
Z <- N <- Y <- spp_pool <- binsAway <- interpRng <- matrix(0, nrow=n.el, ncol=n.spp)
Ytot <- rep(0, n.el)
interpPatchy <- rep(0, n.spp)

# compute occupancy, abundances, detections
rng <- do.call("cbind", map(1:n.spp, ~sort(sample.int(n.el, 2))))
rng.obs <- matrix(0, ncol=n.spp, nrow=2)
for(i in 1:n.el) {
  Z[i,] <- i >= rng[1,] & i <= rng[2,]
  N[i,] <- rpois(n.spp, rpois(n.spp, runif(n.spp, 50, 5000)) *
                            runif(n.spp, 0, 1) * Z[i,])
  
  spp_pool[i,] <- N[i,]*delta
  Ytot[i] <- ceiling(rlnorm(1, log(20), log(4)))*any(Z[i,])
  if(Ytot[i] > 0) {
    Y[i,] <- rmultinom(1, Ytot[i], prob=spp_pool[i,])
  } else {
    Y[i,] <- 0
  }
}
spAbund <- colSums(Y)


# compute interpPatchy & binsAway
for(s in 1:n.spp) {
  # observed interpolated range
  rng.obs[,s] <- range(which(Y[,s]>0))
  interpRng.s <- rng.obs[1,s]:rng.obs[2,s]
  interpRng[,s] <- 1:n.el %in% interpRng.s
  # within interpolated range, proportion of bins with Y.s=0
  interpPatchy[s] <- sum(Y[interpRng.s,s]==0)/length(interpRng.s)
  for(i in (1:n.el)[-interpRng.s]) {
    # number of bins away from last detection (outside interpolated range only)
    binsAway[i,s] <- min(abs(i - rng.obs[,s]))
  }
}

# aggregate true distributions
true.df <- data.frame(el=rep(1:n.el, times=n.spp),
                      spp=paste("Species", rep(1:n.spp, each=n.el)),
                      Z=as.numeric(Z),
                      N=c(N),
                      Y=c(Y),
                      binsAway=c(binsAway),
                      Ytot=rep(Ytot, times=n.spp))

ggplot(true.df, aes(x=el)) + 
  geom_area(aes(y=Z), fill="blue", alpha=0.25) +
  geom_area(aes(y=as.numeric(Y>0)), fill="red", alpha=0.25) +
  facet_wrap(~spp) + labs(x="Elevational bin", y="True presence or absence")




########
## Fit JAGS model
########
# data to feed model
jags_d <- list(n.spp=n.spp, 
               n.el=n.el, 
               Y=Y, 
               Ytot=Ytot, 
               delta=delta, 
               spAbund=spAbund,
               interpPatchy=c(scale(interpPatchy)),
               distAway=matrix(scale(c(binsAway*100)), ncol=n.spp))

# parameters to store
pars <- c("Z", "N", "b")

# fit model
mod <- jags.model(file="code/multinom_b_global.txt", data=jags_d,
                  n.chains=3, n.adapt=1000, 
                  inits=list(Z=matrix(1, n.el, n.spp)))
out <- coda.samples(mod, variable.names=pars, n.iter=2000, thin=5)

# summarize output
gg.b <- ggs(out, "b")
gg.Z <- ggs(out, "Z") %>%
  mutate(el=str_split_fixed(Parameter, ",", 2)[,1] %>%
           str_remove("Z\\[") %>% as.numeric,
         spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
           str_remove("\\]") %>% as.numeric) %>%
  rename(Z=value)
gg.N <- ggs(out, "N") %>%
  mutate(el=str_split_fixed(Parameter, ",", 2)[,1] %>%
           str_remove("N\\[") %>% as.numeric,
         spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
           str_remove("\\]") %>% as.numeric) %>%
  rename(N=value)
gg.NZ <- full_join(gg.Z, gg.N, c("el", "spp", "Chain", "Iteration")) %>%
  mutate(NZ=N*Z,
         spp=paste("Species", spp))
gg.NZ.sum <- gg.NZ %>% group_by(el, spp) %>%
  summarise(mn_prPres=mean(NZ>0),
            mn_N=mean(N),
            mn_Z=mean(Z),
            prPres_gt05=as.numeric(mn_prPres>0.05)) %>%
  full_join(true.df, c("el", "spp"))





########
## Plots
########

# plot: slope probability densities
# ggs_density(gg.b) + facet_wrap(~Parameter) + geom_vline(xintercept=0)

# plot: true Z [blue area], Y [red area], mean(prob N>0) [line]
ggplot(gg.NZ.sum, aes(x=el, y=mn_prPres)) + 
  geom_area(data=true.df, aes(y=Z), fill="blue", alpha=0.25) +
  geom_area(data=true.df, aes(y=as.numeric(Y>0)), fill="red", alpha=0.25) +
  geom_point(aes(colour=log(Ytot+1))) + 
  geom_line() + facet_wrap(~spp) + ylim(0,1) + 
  labs(x="Elevational bin", y="Predicted probability of presence")

# plot: true Z (blue area), Y (red area), mean(probability of presence)>5%
# ggplot(gg.NZ.sum, aes(x=el, y=prPres_gt05)) + 
#   geom_area(data=true.df, aes(y=Z), fill="blue", alpha=0.25) +
#   geom_area(data=true.df, aes(y=as.numeric(Y>0)), fill="red", alpha=0.25) +
#   geom_line() + facet_wrap(~spp) + ylim(0,1) + 
#   labs(x="Elevational bin", y="Predicted presence\n(mean prob > 0.05)")
