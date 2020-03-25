# Simulated data for evaluating model performance and sensitivity
# Tim Szewczyk


library(tidyverse); library(rjags); library(ggmcmc); theme_set(theme_bw())
source("code/00_fn.R")
sp.i <- read_csv("data/by_trapType/Mamm_Summary_Data.csv") # true data
prDet.df <- readxl::read_xlsx("data/orig/Prob_Detection.xlsx", 1) %>%
  select(Species, 5:36) %>% 
  pivot_longer(2:33, names_to="Site", values_to="prDet", values_drop_na=T)

# assumptions
# - relAbundance (on avg) at each trapping location = relAbundance of bin
# - pr(detect individual) is constant within species and across years
# - elevational range = minimum and maximum during sampling time period

# For each species j in elevational bin i:
#   r.mn.ij = b0.j + b1.j * el.i + b2.j * el.i^2
#
# In each year: 
#   r.ijk ~ Norm(r.mn.ij*(1 - N[k-1,]/K), r.sd.j) 

n.sim <- 10

J <- 44  # number of species
tmax <- 200  # number of years to simulate
b <- 50  # bin size
yrs.obs <- (-29:0)+tmax  # years to draw samples from


beta <- c(0.2, -0.3, -0.3)  # mean slopes (avg b0, b1, b2 among species)
beta.sd <- c(0.2, 0.5, 0.1)  # sd in slopes (sd of b0, b1, b2 among species)

Ytot_real <- read_csv(paste0("data/by_trapType/Sampling_", b, "m.csv"))
Ytot <- Ytot_real$ShermanFR_H
els <- Ytot_real$Elev  # elevational bins
n.els <- length(els)



################################################################################
##-- simulation and modelling loop
##------------------------------------------------------------------------------
gg.NZ <- rmse.ls <- vector("list", n.sim)

for(s in 1:n.sim) {
  
  ##-- simulate communities
  comm.true <- simulate_communities(J, els, beta, beta.sd, tmax)
  comm.true$rng.large <- cbind(apply(comm.true$rng[,1,yrs.obs], 1, min, na.rm=T),
                               apply(comm.true$rng[,2,yrs.obs], 1, max, na.rm=T))
  comm.true$rng.med <- cbind(apply(comm.true$rng[,1,yrs.obs], 1, median, na.rm=T),
                             apply(comm.true$rng[,2,yrs.obs], 1, median, na.rm=T))
  
  
  ##-- sample from communities
  obs <- sample_community(J, els, yrs.obs, comm.true, Ytot)
  
  
  ##-- fit model
  jags_d <- list(J=J, 
                 n.el=n.els, 
                 y=obs$Y, 
                 Y=Ytot, 
                 delta=comm.true$pDet, 
                 LAMBDA=obs$spAbund,
                 interpPatchy=c(scale(obs$interpPatchy)),
                 distAway=matrix(scale(c(obs$binsAway*b)), ncol=J))
  pars <- c("lambda", "Z")
  mod <- jags.model(file="code/00_multinom_b_global.txt", data=jags_d,
                    n.chains=4, n.adapt=1000, 
                    inits=list(Z=matrix(1, n.els, J)))
  out <- coda.samples(mod, variable.names=pars, n.iter=2000, thin=5)
  
  
  ##-- aggregate output
  true.df <- data.frame(bin=rep(1:n.els, times=J),
                        Elevation=rep(els, times=J),
                        spp=rep(1:J, each=n.els),
                        N=c(t(apply(comm.true$N[yrs.obs,,], 2:3, mean))),
                        Z=c(t(apply(comm.true$Z[yrs.obs,,], 2:3, mean)))>0.05,
                        Y=c(obs$Y),
                        Ytot=rep(Ytot, times=J),
                        binsAway=c(obs$binsAway),
                        distAway=c(obs$binsAway)*b,
                        distAway_sc=c(scale(c(obs$binsAway)*b)),
                        interpPatchy=rep(obs$interpPatchy, each=n.els),
                        interp.rng=c(obs$interp.rng),
                        spAbund=rep(obs$spAbund, each=n.els),
                        delta=rep(comm.true$pDet, each=n.els))
  gg.NZ[[s]] <- full_join(ggs(out, "Z") %>%
                            mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
                                     str_remove("Z\\[") %>% as.numeric,
                                   spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
                                     str_remove("\\]") %>% as.numeric) %>%
                            rename(Z=value), 
                          ggs(out, "lambda") %>%
                            mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
                                     str_remove("lambda\\[") %>% as.numeric,
                                   spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
                                     str_remove("\\]") %>% as.numeric) %>%
                            rename(N=value), 
                          by=c("bin", "spp", "Chain", "Iteration")) %>%
    mutate(NZ=N*Z,
           spp=as.numeric(spp)) %>% 
    group_by(bin, spp) %>%
    summarise(mn_N=mean(N),
              mn_prPres=mean(NZ>0)) %>%
    full_join(true.df, c("bin", "spp")) %>%
    mutate(binSize=b, 
           sim=s)
  
  rmse.ls[[s]] <- data.frame(spp=1:J, 
                             true.lo=comm.true$rng.med[,1],
                             true.hi=comm.true$rng.med[,2],
                             obs.lo=obs$rng[1,],
                             obs.hi=obs$rng[2,]) %>%
    full_join(gg.NZ[[s]] %>% filter(mn_prPres >= 0.05) %>%
                arrange(spp, bin) %>% group_by(spp) %>%
                summarise(mod.lo=first(bin), 
                          mod.hi=last(bin)),
              by="spp") %>%
    mutate_at(2:7, ~.*b+min(els)) %>% 
    filter(!is.na(true.lo)) %>%
    summarise(obs.lo=sqrt(mean((obs.lo-true.lo)^2, na.rm=T)),
              obs.hi=sqrt(mean((obs.hi-true.hi)^2, na.rm=T)),
              mod.lo=sqrt(mean((mod.lo-true.lo)^2, na.rm=T)),
              mod.hi=sqrt(mean((mod.hi-true.hi)^2, na.rm=T))) %>%
    mutate(sim=s)
  cat("\n--------------\n", 
      "Finished", s, "of", n.sim,
      "\n--------------\n")
}


write.csv(do.call('rbind', rmse.ls), paste0("out/RMSE_", b, ".csv"))
do.call('rbind', rmse.ls) %>% 
  summarise(mnDiff.lo=mean(mod.lo-obs.lo, na.rm=T), 
            mnDiff.hi=mean(mod.hi-obs.hi, na.rm=T), 
            pctDiff.lo=mean((mod.lo-obs.lo)/obs.lo, na.rm=T)*100,
            pctDiff.hi=mean((mod.hi-obs.hi)/obs.hi, na.rm=T)*100)






