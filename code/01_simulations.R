# Simulated data for evaluating model performance and sensitivity
# Tim Szewczyk


library(tidyverse); library(rjags); library(ggmcmc); theme_set(theme_bw())
source("code/00_fn.R")
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

sets <- expand.grid(b=c(25, 50, 100, 200),
                    mtn=c("FR", "SJ"), 
                    effort=c(0.5, 1, 2)) # proportional to empirical

n.sim <- 50
tmax <- 200  # number of years to simulate
yrs.obs <- (-29:0)+tmax  # years to draw samples from

beta <- c(0.2, -0.3, -0.3)  # mean slopes (avg b0, b1, b2 among species)
beta.sd <- c(0.2, 0.5, 0.1)  # sd in slopes (sd of b0, b1, b2 among species)


for(set in 1:nrow(sets)) {
  
  # parameters for 'set' details
  b <- sets$b[set]  
  mtn <- sets$mtn[set]
  type <- "all"
  effort <- sets$effort[set]
  
  Y.real.sample <- read_csv("data/sample_els_H.csv") %>%
    mutate(elBin=el %/% b * b,
           Sample=sp.i$Sample[match(sp, sp.i$Abbrev)]) 
  if(mtn=="FR") {
    Y.real.sample <- Y.real.sample %>% 
      filter(county %in% c("Boulder", "Larimer"))
  } else if(mtn=="SJ") {
    Y.real.sample <- Y.real.sample %>%
      filter(county %in% c("Dolores", "La Plata", "Montezuma", "San Juan"))
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
  gg.NZ <- gg.delta <- gg.beta <- rmse.ls <- vector("list", n.sim)
  
  for(s in 1:n.sim) {
    
    ##-- simulate communities
    comm.true <- simulate_communities(J, els, pDet.set, beta, beta.sd, tmax)
    comm.true$rng.large <- cbind(apply(comm.true$rng[,1,yrs.obs], 1, 
                                       min, na.rm=T),
                                 apply(comm.true$rng[,2,yrs.obs], 1, 
                                       max, na.rm=T))
    comm.true$rng.med <- cbind(apply(comm.true$rng[,1,yrs.obs], 1,
                                     median, na.rm=T),
                               apply(comm.true$rng[,2,yrs.obs], 1, 
                                     median, na.rm=T))
    
    
    ##-- sample from communities
    obs <- sample_community(J, els, yrs.obs, comm.true, Ytot)
    
    
    ##-- fit model
    jags_d <- list(J=J, 
                   n.el=n.els, 
                   y=obs$Y, 
                   Y=Ytot, 
                   delta_shp=as.matrix(comm.true$pDet.par), 
                   LAMBDA=obs$spAbund,
                   interpPatchy=c(scale(obs$interpPatchy)),
                   distAway=matrix(scale(c(obs$binsAway*b)), ncol=J))
    pars <- c("lambda", "Z", "beta")
    mod <- jags.model(file="code/00_multinom_b_global.txt", data=jags_d,
                      n.chains=4, n.adapt=2000, quiet=T,
                      inits=list(Z=matrix(1, n.els, J)))
    out <- coda.samples(mod, variable.names=pars, n.iter=1000, thin=5)
    
    
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
                          delta=c(comm.true$pDet.el))
    gg.NZ[[s]] <- full_join(ggs(out, "Z") %>%
                              mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
                                       str_remove("Z\\[") %>% as.numeric,
                                     spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
                                       str_remove("\\]") %>% as.numeric) %>%
                              rename(Z.mod=value), 
                            ggs(out, "lambda") %>%
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
      mutate(sim=s, b=b, mtn=mtn, effort=effort)
    
    # gg.delta[[s]] <- ggs(out, "delta") %>%
    #   mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
    #            str_remove("delta\\[") %>% as.numeric,
    #          spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
    #            str_remove("\\]") %>% as.numeric) %>% 
    #   full_join(., true.df, by=c("bin", "spp")) %>%
    #   mutate(sim=s, b=b, mtn=mtn, effort=effort)
    
    gg.beta[[s]] <- ggs(out, "beta") %>%
      mutate(sim=s, b=b, mtn=mtn, effort=effort)
    
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
      summarise(obs.lo=sqrt(mean((obs.lo-true.lo)^2, na.rm=T)),
                obs.hi=sqrt(mean((obs.hi-true.hi)^2, na.rm=T)),
                mod.lo=sqrt(mean((mod.lo-true.lo)^2, na.rm=T)),
                mod.hi=sqrt(mean((mod.hi-true.hi)^2, na.rm=T))) %>%
      mutate(sim=s, b=b, mtn=mtn, effort=effort)
    cat("\n--------------\n", 
        "Finished", s, "of", n.sim, "in set", set,
        "\n--------------\n")
  }
  
  set_pars <- paste(paste0(names(sets), c(as.character(b), 
                                          as.character(mtn), 
                                          as.character(type), 
                                          as.character(effort*100))), 
                    collapse="_")
  
  write_csv(do.call('rbind', rmse.ls), 
            paste0("out/RMSE_", set_pars, ".csv"))
  write_csv(do.call('rbind', gg.NZ),
            paste0("out/ggNZ_", set_pars, ".csv"))
  write_csv(do.call('rbind', gg.delta),
            paste0("out/ggDelta_", set_pars, ".csv"))
  write_csv(do.call('rbind', gg.beta),
            paste0("out/ggBeta_", set_pars, ".csv"))
}


RMSE <- map_dfr(dir("out", "RMSE", full.names=T), read_csv)


RMSE %>% replace(x=., list=.==0, NA) %>%
  group_by(b, type, effort) %>%
  summarise(mnDiff.lo=mean(mod.lo-obs.lo, na.rm=T), 
            mnDiff.hi=mean(mod.hi-obs.hi, na.rm=T), 
            pctDiff.lo=mean((mod.lo-obs.lo)/obs.lo, na.rm=T)*100,
            pctDiff.hi=mean((mod.hi-obs.hi)/obs.hi, na.rm=T)*100)

ggplot(RMSE, aes(x=(mod.lo-obs.lo)/obs.lo, colour=factor(b))) + 
  geom_vline(xintercept=0, linetype=2) + 
  geom_density() + facet_grid(effort~type, scales="free")
ggplot(RMSE, aes(x=(mod.hi-obs.hi)/obs.hi, colour=factor(b))) + 
  geom_vline(xintercept=0, linetype=2) + 
  geom_density() + facet_grid(effort~type, scales="free")

RMSE %>% mutate(diff.lo=mod.lo-obs.lo, diff.hi=mod.hi-obs.hi) %>%
  select(5:11) %>%
  pivot_longer(6:7, names_to="boundary", values_to="diff") %>%
  mutate(boundary=factor(str_sub(boundary, -2L, -1L), levels=c("lo", "hi")), 
         b=factor(b, levels=c("25", "50", "100", "200")),
         effort=factor(effort, levels=c(0.5, 1, 2), 
                       labels=c("50%", "100%", "200%"))) %>%
  group_by(b, boundary, mtn, type, effort) %>%
  summarise(mnDiff=mean(diff), seDiff=sd(diff)/sqrt(max(sim))) %>%
  ggplot(aes(x=effort, y=mnDiff, ymin=mnDiff-2*seDiff, ymax=mnDiff+2*seDiff,
             colour=type)) +
  geom_hline(yintercept=0, linetype=2) + 
  geom_point(position=position_dodge(width=0.25)) + 
  geom_linerange(position=position_dodge(width=0.25)) + 
  facet_grid(mtn~boundary) + 
  scale_colour_brewer("Trap Type", type="qual", palette=2) +
  labs(x="Elevational bin size", y="Change in RMSE (m)")

RMSE %>% mutate(pct.lo=(mod.lo-obs.lo)/obs.lo*100, 
                pct.hi=(mod.hi-obs.hi)/obs.hi*100,
                pct.lo=na_if(pct.lo, Inf),
                pct.hi=na_if(pct.hi, Inf)) %>%
  select(5:11) %>%
  pivot_longer(6:7, names_to="boundary", values_to="diff") %>%
  mutate(boundary=factor(str_sub(boundary, -2L, -1L), levels=c("lo", "hi")), 
         b=factor(b, levels=c("25", "50", "100", "200")),
         effort=factor(effort, levels=c(0.5, 1, 2), 
                       labels=c("50%", "100%", "200%"))) %>%
  group_by(b, boundary, mtn, type, effort) %>%
  summarise(mnDiff=mean(diff, na.rm=T), 
            seDiff=sd(diff, na.rm=T)/sqrt(max(sim))) %>%
  ggplot(aes(x=effort, y=mnDiff, ymin=mnDiff-2*seDiff, ymax=mnDiff+2*seDiff,
             colour=type)) +
  geom_hline(yintercept=0, linetype=2) + 
  geom_point(position=position_dodge(width=0.25)) + 
  geom_linerange(position=position_dodge(width=0.25)) + 
  facet_grid(mtn~boundary) + 
  scale_colour_brewer("Trap Type", type="qual", palette=2) +
  labs(x="Elevational bin size", y="Change in RMSE (%)")

RMSE %>% pivot_longer(1:4, names_to="est", values_to="RMSE") %>%
  mutate(boundary=str_sub(est, -2L, -1L), est=str_sub(est, 1, 3)) %>%
  ggplot(aes(x=boundary, y=RMSE, fill=paste(mtn, est))) + geom_boxplot() + 
  facet_grid(type~b) + 
  scale_fill_manual(values=c("#a6611a", "#018571", "#dfc27d", "#80cdc1"))

# 
# do.call('rbind', gg.NZ) %>%
#   ggplot(aes(x=Elevation)) + facet_grid(sim~spp) +
#   geom_ribbon(aes(xmin=Elevation, xmax=Elevation,
#                   ymin=0, ymax=as.numeric(Z)), alpha=0.25) +
#   geom_line(aes(y=as.numeric(interp.rng)), colour="blue") +
#   geom_line(aes(y=as.numeric(mn_prPres>0.05)), colour="red")
do.call('rbind', gg.NZ) %>% 
  ggplot(aes(mn_N.mod, log(N))) + geom_point(alpha=0.05)



ggBeta <- map_dfr(dir("out", "ggBeta", full.names=T), read_csv)

ggplot(filter(ggBeta, Parameter=="beta[1]"), 
       aes(x=value, colour=factor(b))) + 
