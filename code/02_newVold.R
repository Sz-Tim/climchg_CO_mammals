# Cross-validation for old vs. new models
# Tim Szewczyk

# This script splits data into testing/training sets so that the model can be
# fit with the training set and the ability to predict the testing set can be
# evaluated. This is not as straightforward as with a simpler model like a 
# linear regression, since the entire goal of the model is to estimate 
# uncertainty about unobserved data, and it relies on the observations of the
# entire community. Further, since each mtn x era is fit independently, fitting
# with one mountain range or era and predicting the other is not really 
# consistent with the model usage or structure.

# Instead, we can use the simulated data 

library(tidyverse); library(rjags); library(ggmcmc); library(dclone)
source("code/00_fn.R"); theme_set(theme_bw())
sp.i <- read_csv("data/by_trapType/Mamm_Summary_Data.csv") # true data
pDet.full <- read_csv("data/prDet_processed.csv")


# assumptions
# - relAbundance (on avg) at each trapping location = relAbundance of bin
# - pr(detect individual) is drawn from a species-specific distribution
# - elevational range = minimum and maximum during sampling time period

# For each species j in elevational bin i:
#   r.mn.ij = b0.j + b1.j * el.i + b2.j * el.i^2
#
# In each year: 
#   r.ijk ~ Norm(r.mn.ij*(1 - N[k-1,]/K), r.sd.j) 

sets <- expand.grid(b=50, 
                    mtn=c("FR", "SJ"), 
                    type="all",
                    effort=1,
                    era=c("H", "C")) 


n.sim <- 50
tmax <- 200  # number of years to simulate
yrs.obs <- (-29:0)+tmax  # years to draw samples from

beta <- c(0.2, -0.3, -0.3)  # mean slopes (avg b0, b1, b2 among species)
beta.sd <- c(0.2, 0.5, 0.15)  # sd in slopes (sd of b0, b1, b2 among species)


for(set in 1:nrow(sets)) {
  
  # parameters for 'set' details
  b <- sets$b[set]  
  mtn <- sets$mtn[set]
  type <- sets$type[set]
  effort <- sets$effort[set]
  era <- sets$era[set]
  
  Y.real.sample <- read_csv(paste0("data/sample_els_", era, ".csv")) %>%
    mutate(elBin=el %/% b * b,
           Sample=sp.i$Sample[match(sp, sp.i$Abbrev)]) 
  if(mtn=="FR") {
    Y.real.sample <- Y.real.sample %>% 
      filter(county %in% c("Boulder", "Larimer", "Big Thompson", "Cedar Park",
                           "Cow Camp", "Denver", "SylvanDale"))
  } else if(mtn=="SJ") {
    Y.real.sample <- Y.real.sample %>%
      filter(county %in% c("Dolores", "La Plata", "Montezuma", "San Juan",
                           "East San Juans", "West San Juans"))
  }
  Y.real <- Y.real.sample %>% group_by(elBin) %>% summarise(Y=round(n()*effort))
  els <- seq(min(Y.real$elBin), max(Y.real$elBin), by=b)
  n.els <- length(els)
  Ytot <- Y.real$Y[match(els, Y.real$elBin)] %>% replace_na(0)
  J <- n_distinct(Y.real.sample$sp)
  pDet.set <- pDet.full %>% 
    filter(Species %in% filter(sp.i, Abbrev %in% Y.real.sample$sp)$Species)
  
  
  
  
  ################################################################################
  ##-- simulation and modelling loop
  ##------------------------------------------------------------------------------
  gg.NZ.new <- gg.NZ.old <- rmse.ls <- vector("list", n.sim)
  
  ##-- simulate communities
  comm.true <- simulate_communities(J, els, pDet.set, beta, beta.sd, tmax)
  
  for(s in 1:n.sim) {
    
    ##-- sample from communities
    obs <- sample_community(J, els, yrs.obs, comm.true, Ytot)
    
    
    ##-- fit model
    jags_d <- list(J=J, 
                   n.el=n.els, 
                   y=obs$Y, 
                   Y=rowSums(obs$Y), 
                   delta_dat=comm.true$pDet.par[,1]/rowSums(comm.true$pDet.par),
                   delta_shp=as.matrix(comm.true$pDet.par), 
                   LAMBDA=obs$spAbund,
                   interpPatchy=c(scale(obs$interpPatchy)),
                   distAway=matrix(scale(c(obs$binsAway*b)), ncol=J))
    pars <- c("lambda", "Z")
    cl <- makeCluster(4)
    out.new <- jags.parfit(cl=cl, data=jags_d, params=pars,
                           model="code/00_multinom_b_global.txt", 
                           inits=list(Z=matrix(1, n.els, J)),
                           n.thin=5, n.chains=4, n.adapt=5000, n.update=1000)
    stopCluster(cl)
    
    cl <- makeCluster(4)
    out.old <- jags.parfit(cl=cl, data=jags_d, params=pars,
                           model="code/00_multinom_b_global_deltaAsData.txt", 
                           inits=list(Z=matrix(1, n.els, J)),
                           n.thin=5, n.chains=4, n.adapt=5000, n.update=1000)
    stopCluster(cl)
    
    
    
    ##-- aggregate output
    true.df <- data.frame(bin=rep(1:n.els, times=J),
                          Elevation=rep(els, times=J),
                          spp=rep(1:J, each=n.els),
                          N=c(t(apply(comm.true$N[yrs.obs,,], 2:3, mean))),
                          Z=c(t(apply(comm.true$Z[yrs.obs,,], 2:3, mean)))>0.05,
                          Y=c(obs$Y),
                          Ytot=rep(rowSums(obs$Y), times=J),
                          binsAway=c(obs$binsAway),
                          distAway=c(obs$binsAway)*b,
                          distAway_sc=c(scale(c(obs$binsAway)*b)),
                          interpPatchy=rep(obs$interpPatchy, each=n.els),
                          interp.rng.obs=c(obs$interp.rng),
                          interp.rng.true=c(comm.true$interp.rng.true),
                          spAbund=rep(obs$spAbund, each=n.els),
                          delta=c(t(apply(comm.true$pDet.ar[yrs.obs,,], 
                                          2:3, mean))))
    gg.NZ.old[[s]] <- full_join(ggs(out.old, "Z") %>%
                                  mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
                                           str_remove("Z\\[") %>% as.numeric,
                                         spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
                                           str_remove("\\]") %>% as.numeric) %>%
                                  rename(Z.mod=value), 
                                ggs(out.old, "lambda") %>%
                                  mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
                                           str_remove("lambda\\[") %>% as.numeric,
                                         spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
                                           str_remove("\\]") %>% as.numeric) %>%
                                  rename(N.mod=value), 
                                by=c("bin", "spp", "Chain", "Iteration")) %>%
      mutate(NZ.mod=N.mod*Z.mod,
             spp=as.numeric(spp)) %>% 
      group_by(bin, spp) %>%
      summarise(mn_N.mod=mean(N.mod),
                mn_prPres=mean(NZ.mod>0)) %>%
      full_join(., true.df, by=c("bin", "spp")) %>%
      mutate(sim=s, b=b, mtn=mtn, type=type, effort=effort, era=era)
    gg.NZ.new[[s]] <- full_join(ggs(out.new, "Z") %>%
                                  mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
                                           str_remove("Z\\[") %>% as.numeric,
                                         spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
                                           str_remove("\\]") %>% as.numeric) %>%
                                  rename(Z.mod=value), 
                                ggs(out.new, "lambda") %>%
                                  mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
                                           str_remove("lambda\\[") %>% as.numeric,
                                         spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
                                           str_remove("\\]") %>% as.numeric) %>%
                                  rename(N.mod=value), 
                                by=c("bin", "spp", "Chain", "Iteration")) %>%
      mutate(NZ.mod=N.mod*Z.mod,
             spp=as.numeric(spp)) %>% 
      group_by(bin, spp) %>%
      summarise(mn_N.mod=mean(N.mod),
                mn_prPres=mean(NZ.mod>0)) %>%
      full_join(., true.df, by=c("bin", "spp")) %>%
      mutate(sim=s, b=b, mtn=mtn, type=type, effort=effort, era=era)
    
    rmse.ls[[s]] <- data.frame(spp=1:J, 
                               true.lo=comm.true$rng.med[,1],
                               true.hi=comm.true$rng.med[,2],
                               obs.lo=na_if(obs$rng[1,], 0),
                               obs.hi=na_if(obs$rng[2,], 0)) %>%
      full_join(gg.NZ.old[[s]] %>% filter(mn_prPres >= 0.05) %>%
                  arrange(spp, bin) %>% group_by(spp) %>%
                  summarise(mod.old.lo=first(bin), 
                            mod.old.hi=last(bin)),
                by="spp") %>%
      full_join(gg.NZ.new[[s]] %>% filter(mn_prPres >= 0.05) %>%
                  arrange(spp, bin) %>% group_by(spp) %>%
                  summarise(mod.new.lo=first(bin), 
                            mod.new.hi=last(bin)),
                by="spp") %>%
      mutate_at(2:9, ~.*b+min(els)) %>% 
    mutate(sim=s, b=b, mtn=mtn, type=type, effort=effort, era=era)
    cat("--------------\n", 
        "Finished", s, "of", n.sim, "in set", set, "/", nrow(sets),
        "\n--------------\n")
  }
  
  set_pars <- paste(paste0(names(sets), c(as.character(b), 
                                          as.character(mtn), 
                                          as.character(type),
                                          as.character(effort*100),
                                          as.character(era))), 
                    collapse="_")
  
  write_csv(do.call('rbind', rmse.ls), 
            paste0("out/newVold_RMSE_", set_pars, ".csv"))
  # write_csv(do.call('rbind', gg.NZ.old),
  #           paste0("out/ggNZ_old_", set_pars, ".csv"))
  # write_csv(do.call('rbind', gg.NZ.new),
  #           paste0("out/ggNZ_new_", set_pars, ".csv"))
}


RMSE <- map_dfr(dir("out", "newVold_RMSE", full.names=T), read_csv)
write_csv(RMSE, "out/newVold_out.csv")


RMSE %>% filter(!is.na(obs.lo)) %>%
  group_by(mtn, era, sim) %>%
  summarise(obs.lo.rmse=sqrt(mean((obs.lo-true.lo)^2, na.rm=T)),
            mod.old.lo.rmse=sqrt(mean((mod.old.lo-true.lo)^2, na.rm=T)),
            mod.new.lo.rmse=sqrt(mean((mod.new.lo-true.lo)^2, na.rm=T)),
            obs.hi.rmse=sqrt(mean((obs.hi-true.hi)^2, na.rm=T)),
            mod.old.hi.rmse=sqrt(mean((mod.old.hi-true.hi)^2, na.rm=T)),
            mod.new.hi.rmse=sqrt(mean((mod.new.hi-true.hi)^2, na.rm=T))) %>%
  mutate(diff.old.lo=mod.old.lo.rmse-obs.lo.rmse,
         diff.new.lo=mod.new.lo.rmse-obs.lo.rmse,
         diff.old.hi=mod.old.hi.rmse-obs.hi.rmse,
         diff.new.hi=mod.new.hi.rmse-obs.hi.rmse) %>%
  select(1:3,10:13) %>% ungroup %>%
  pivot_longer(4:7, names_to="bound", values_to="diff") %>%
  mutate(boundary=factor(str_sub(bound, -2L, -1L), levels=c("lo", "hi")), 
         model=str_sub(bound, 6, 8)) %>%
  group_by(boundary, mtn, model, era) %>%
  summarise(mnDiff=mean(diff), seDiff=sd(diff)/sqrt(max(sim))) %>%
  ggplot(aes(x=mtn, y=mnDiff, ymin=mnDiff-2*seDiff, ymax=mnDiff+2*seDiff,
             colour=model)) +
  geom_hline(yintercept=0, linetype=2) + 
  geom_point(position=position_dodge(width=0.25)) + 
  geom_linerange(position=position_dodge(width=0.25)) + 
  facet_grid(era~boundary) + 
  scale_colour_brewer("Model\nversion", type="qual", palette=2) +
  labs(x="Mountain", y="Change in RMSE (m)")

RMSE %>% filter(!is.na(obs.lo)) %>%
  group_by(mtn, era, sim) %>%
  summarise(obs.lo.rmse=sqrt(mean((obs.lo-true.lo)^2, na.rm=T)),
            mod.old.lo.rmse=sqrt(mean((mod.old.lo-true.lo)^2, na.rm=T)),
            mod.new.lo.rmse=sqrt(mean((mod.new.lo-true.lo)^2, na.rm=T)),
            obs.hi.rmse=sqrt(mean((obs.hi-true.hi)^2, na.rm=T)),
            mod.old.hi.rmse=sqrt(mean((mod.old.hi-true.hi)^2, na.rm=T)),
            mod.new.hi.rmse=sqrt(mean((mod.new.hi-true.hi)^2, na.rm=T))) %>%
  mutate(diff.old.lo=(mod.old.lo.rmse-obs.lo.rmse)/obs.lo.rmse,
         diff.new.lo=(mod.new.lo.rmse-obs.lo.rmse)/obs.lo.rmse,
         diff.old.hi=(mod.old.hi.rmse-obs.hi.rmse)/obs.hi.rmse,
         diff.new.hi=(mod.new.hi.rmse-obs.hi.rmse)/obs.hi.rmse) %>%
  select(1:3,10:13) %>% ungroup %>%
  pivot_longer(4:7, names_to="bound", values_to="diff") %>%
  mutate(boundary=factor(str_sub(bound, -2L, -1L), levels=c("lo", "hi"),
                         labels=c("Lower boundary", "Upper boundary")), 
         model=factor(str_sub(bound, 6, 8), levels=c("old", "new")),
         era=factor(era, levels=c("H", "C"), 
                    labels=c("Historical", "Contemporary"))) %>%
  group_by(boundary, mtn, model, era) %>%
  summarise(mnDiff=mean(diff), seDiff=sd(diff)/sqrt(max(sim))) %>%
  ggplot(aes(x=mtn, y=mnDiff, ymin=mnDiff-2*seDiff, ymax=mnDiff+2*seDiff,
             colour=model)) +
  geom_hline(yintercept=0, linetype=2) + 
  geom_point(position=position_dodge(width=0.25)) + 
  geom_linerange(position=position_dodge(width=0.25)) + 
  facet_grid(era~boundary) + 
  scale_y_continuous(labels=scales::percent) +
  scale_colour_brewer("Model\nversion", type="qual", palette=2, 
                      breaks=c("old", "new"), 
                      labels=c(bquote(delta~fixed), 
                               bquote(delta~free))) +
  labs(x="Mountain Range", 
       y="Change in RMSE relative to empirical (mean ± 2 SE)")
ggsave("figs/RMSE_newVold_pctChg.pdf", width=6, height=5)





RMSE %>% filter(!is.na(obs.lo)) %>%
  group_by(mtn, era, sim) %>%
  summarise(obs.lo.rmse=sqrt(mean((obs.lo-true.lo)^2, na.rm=T)),
            mod.old.lo.rmse=sqrt(mean((mod.old.lo-true.lo)^2, na.rm=T)),
            mod.new.lo.rmse=sqrt(mean((mod.new.lo-true.lo)^2, na.rm=T)),
            obs.hi.rmse=sqrt(mean((obs.hi-true.hi)^2, na.rm=T)),
            mod.old.hi.rmse=sqrt(mean((mod.old.hi-true.hi)^2, na.rm=T)),
            mod.new.hi.rmse=sqrt(mean((mod.new.hi-true.hi)^2, na.rm=T))) %>%
  write_csv("out/newVold_RMSE.csv")


