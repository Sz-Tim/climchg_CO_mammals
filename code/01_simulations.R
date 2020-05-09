# Simulated data for evaluating model performance and sensitivity
# Tim Szewczyk


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

sets <- expand.grid(b=c(200, 100, 50),
                    mtn=c("FR", "SJ")[2], 
                    type="all",#c("Large", "Sherman+Shrew"),
                    effort=c(0.5, 1, 2),
                    era=c("H", "C")[2]) # proportional to empirical)
                    

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
  if(type=="Large") {
    Y.real.sample <- filter(Y.real.sample, Sample=="Large")
  } else if(type=="Sherman+Shrew") {
    Y.real.sample <- filter(Y.real.sample, Sample!="Large")
  }
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
  gg.NZ <- gg.delta <- gg.beta <- rmse.ls <- vector("list", n.sim)
  
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
                   delta_shp=as.matrix(comm.true$pDet.par), 
                   LAMBDA=obs$spAbund,
                   interpPatchy=c(scale(obs$interpPatchy)),
                   distAway=matrix(scale(c(obs$binsAway*b)), ncol=J))
    pars <- c("lambda", "Z")
    cl <- makeCluster(4)
    out <- jags.parfit(cl=cl, data=jags_d, params=pars,
                       model="code/00_multinom_b_global.txt", n.thin=5,
                       inits=list(Z=matrix(1, n.els, J)),
                       n.chains=4, n.adapt=1000, n.update=1000)
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
      mutate(sim=s, b=b, mtn=mtn, type=type, effort=effort, era=era)
    
    # gg.delta[[s]] <- ggs(out, "delta") %>%
    #   mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
    #            str_remove("delta\\[") %>% as.numeric,
    #          spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
    #            str_remove("\\]") %>% as.numeric) %>% 
    #   full_join(., true.df, by=c("bin", "spp")) %>%
    #   mutate(sim=s, b=b, mtn=mtn, effort=effort)
    
    # gg.beta[[s]] <- ggs(out, "beta") %>%
      # mutate(sim=s, b=b, mtn=mtn, type=type, effort=effort)
    
    rmse.ls[[s]] <- data.frame(spp=1:J, 
                               true.lo=comm.true$rng.med[,1],
                               true.hi=comm.true$rng.med[,2],
                               obs.lo=na_if(obs$rng[1,], 0),
                               obs.hi=na_if(obs$rng[2,], 0)) %>%
      full_join(gg.NZ[[s]] %>% filter(mn_prPres >= 0.05) %>%
                  arrange(spp, bin) %>% group_by(spp) %>%
                  summarise(mod.lo=first(bin), 
                            mod.hi=last(bin)),
                by="spp") %>%
      mutate_at(2:7, ~.*b+min(els)) %>% 
      # summarise(
      #   # pOut.obs.lo=sum(obs.lo < true.lo, na.rm=T)/n(),
      # #           pOut.obs.hi=sum(obs.hi > true.hi, na.rm=T)/n(),
      # #           pOut.mod.lo=sum(mod.lo < true.lo, na.rm=T)/n(),
      # #           pOut.mod.hi=sum(mod.hi > true.hi, na.rm=T)/n(),
      # #           pExact.obs.lo=sum(obs.lo == true.lo, na.rm=T)/n(),
      # #           pExact.obs.hi=sum(obs.hi == true.hi, na.rm=T)/n(),
      # #           pExact.mod.lo=sum(mod.lo == true.lo, na.rm=T)/n(),
      # #           pExact.mod.hi=sum(mod.hi == true.hi, na.rm=T)/n(),
      #           obs.lo=sqrt(mean((obs.lo-true.lo)^2, na.rm=T)),
      #           obs.hi=sqrt(mean((obs.hi-true.hi)^2, na.rm=T)),
      #           mod.lo=sqrt(mean((mod.lo-true.lo)^2, na.rm=T)),
      #           mod.hi=sqrt(mean((mod.hi-true.hi)^2, na.rm=T))) %>%
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
            paste0("out/RMSE_", set_pars, ".csv"))
  write_csv(do.call('rbind', gg.NZ),
            paste0("out/ggNZ_", set_pars, ".csv"))
  # write_csv(do.call('rbind', gg.delta),
  #           paste0("out/ggDelta_", set_pars, ".csv"))
  # write_csv(do.call('rbind', gg.beta),
            # paste0("out/ggBeta_", set_pars, ".csv"))
}


RMSE <- map_dfr(dir("out", "RMSE", full.names=T), read_csv)


RMSE %>% filter(!is.na(obs.lo)) %>%
  group_by(b, mtn, type, effort, era, sim) %>%
  summarise(obs.lo.rmse=sqrt(mean((obs.lo-true.lo)^2, na.rm=T)),
            mod.lo.rmse=sqrt(mean((mod.lo-true.lo)^2, na.rm=T)),
            obs.hi.rmse=sqrt(mean((obs.hi-true.hi)^2, na.rm=T)),
            mod.hi.rmse=sqrt(mean((mod.hi-true.hi)^2, na.rm=T))) %>%
  mutate(diff.lo=mod.lo.rmse-obs.lo.rmse,
         diff.hi=mod.hi.rmse-obs.hi.rmse) %>%
  select(1:6,11:12) %>% ungroup %>%
  pivot_longer(7:8, names_to="boundary", values_to="diff") %>%
  mutate(boundary=factor(str_sub(boundary, -2L, -1L), levels=c("lo", "hi")), 
         b=factor(b, levels=c("25", "50", "100", "200")),
         effort=factor(effort, levels=c(0.5, 1, 2), 
                       labels=c("50%", "100%", "200%"))) %>%
  group_by(b, boundary, mtn, effort, era) %>%
  summarise(mnDiff=mean(diff), seDiff=sd(diff)/sqrt(max(sim))) %>%
  ggplot(aes(x=effort, y=mnDiff, ymin=mnDiff-2*seDiff, ymax=mnDiff+2*seDiff,
             colour=b, shape=mtn)) +
  geom_hline(yintercept=0, linetype=2) + 
  geom_point(position=position_dodge(width=0.25)) + 
  geom_linerange(position=position_dodge(width=0.25)) + 
  facet_grid(era~boundary) + 
  scale_colour_brewer("Bin size", type="qual", palette=2) +
  labs(x="Effort", y="Change in RMSE (m)")
  

RMSE %>% filter(!is.na(obs.lo)) %>%
  group_by(b, mtn, type, effort, era, sim) %>%
  summarise(obs.lo.rmse=sqrt(mean((obs.lo-true.lo)^2, na.rm=T)),
            mod.lo.rmse=sqrt(mean((mod.lo-true.lo)^2, na.rm=T)),
            obs.hi.rmse=sqrt(mean((obs.hi-true.hi)^2, na.rm=T)),
            mod.hi.rmse=sqrt(mean((mod.hi-true.hi)^2, na.rm=T))) %>%
  mutate(diff.lo=(mod.lo.rmse-obs.lo.rmse)/obs.lo.rmse,
         diff.hi=(mod.hi.rmse-obs.hi.rmse)/obs.hi.rmse) %>%
  select(1:6,11:12) %>% ungroup %>%
  pivot_longer(7:8, names_to="boundary", values_to="diff") %>%
  mutate(boundary=factor(str_sub(boundary, -2L, -1L), levels=c("lo", "hi"),
                         labels=c("Lower boundary", "Upper boundary")), 
         b=factor(b, levels=c("25", "50", "100", "200")),
         effort=factor(effort, levels=c(0.5, 1, 2), 
                       labels=c("50%", "100%", "200%")),
         era=factor(era, levels=c("H", "C"), 
                    labels=c("Historical", "Contemporary"))) %>%
  group_by(b, boundary, mtn, effort, era) %>%
  summarise(mnDiff=mean(diff, na.rm=T), 
            seDiff=sd(diff, na.rm=T)/sqrt(max(sim))) %>%
  ggplot(aes(x=effort, y=mnDiff, ymin=mnDiff-2*seDiff, ymax=mnDiff+2*seDiff,
             colour=b, shape=mtn)) +
  geom_hline(yintercept=0, linetype=2) + 
  geom_point(position=position_dodge(width=0.25), size=1.5) + 
  geom_linerange(position=position_dodge(width=0.25)) + 
  facet_grid(era~boundary) + 
  scale_shape_manual("Mountain\nRange", values=c(1,5)) +
  scale_y_continuous(labels=scales::percent) +
  scale_colour_brewer("Bin size", type="qual", palette=2) +
  labs(x="Sampling effort relative to empirical", 
       y="Change in RMSE relative to empirical (mean ± 2 SE)")
ggsave("figs/RMSE_sims_pctChg.pdf", width=6, height=5)


write_csv(RMSE, "out/simulations_out.csv")

RMSE %>% filter(!is.na(obs.lo)) %>%
  group_by(b, mtn, type, effort, era, sim) %>%
  summarise(obs.lo.rmse=sqrt(mean((obs.lo-true.lo)^2, na.rm=T)),
            mod.lo.rmse=sqrt(mean((mod.lo-true.lo)^2, na.rm=T)),
            obs.hi.rmse=sqrt(mean((obs.hi-true.hi)^2, na.rm=T)),
            mod.hi.rmse=sqrt(mean((mod.hi-true.hi)^2, na.rm=T))) %>%
  write_csv("out/simulations_RMSE.csv")
  
  mutate(diff.lo=(mod.lo.rmse-obs.lo.rmse)/obs.lo.rmse,
         diff.hi=(mod.hi.rmse-obs.hi.rmse)/obs.hi.rmse) %>%
  select(1:6,11:12) %>% ungroup %>%
  pivot_longer(7:8, names_to="boundary", values_to="diff") %>%
  mutate(boundary=factor(str_sub(boundary, -2L, -1L), levels=c("lo", "hi"),
                         labels=c("Lower boundary", "Upper boundary")), 
         b=factor(b, levels=c("25", "50", "100", "200")),
         effort=factor(effort, levels=c(0.5, 1, 2), 
                       labels=c("50%", "100%", "200%")),
         era=factor(era, levels=c("H", "C"), 
                    labels=c("Historical", "Contemporary"))) %>%
  group_by(b, boundary, mtn, effort, era) %>%
  summarise(mnDiff=mean(diff, na.rm=T), 
            seDiff=sd(diff, na.rm=T)/sqrt(max(sim)),
            minDiff=min(diff, na.rm=T), 
            maxDiff=max(diff, na.rm=T))

# RMSE %>% pivot_longer(1:4, names_to="est", values_to="RMSE") %>%
#   mutate(boundary=str_sub(est, -2L, -1L), est=str_sub(est, 1, 3)) %>%
#   ggplot(aes(x=boundary, y=RMSE, fill=paste(mtn, est))) + geom_boxplot() + 
#   facet_grid(effort~b) + 
#   scale_fill_manual(values=c("#a6611a", "#018571", "#dfc27d", "#80cdc1"))

# 
# do.call('rbind', gg.NZ) %>%
#   ggplot(aes(x=Elevation)) + facet_grid(sim~spp) +
#   geom_ribbon(aes(xmin=Elevation, xmax=Elevation,
#                   ymin=0, ymax=as.numeric(Z)), alpha=0.25) +
#   geom_line(aes(y=as.numeric(interp.rng)), colour="blue") +
#   geom_line(aes(y=as.numeric(mn_prPres>0.05)), colour="red")
# do.call('rbind', gg.NZ) %>% 
#   ggplot(aes(mn_N.mod, log(N))) + geom_point(alpha=0.05)



ggBeta <- map_dfr(dir("out", "ggBeta", full.names=T), read_csv)

ggplot(filter(ggBeta, Parameter=="beta[1]"), 
       aes(x=value, group=paste(sim, Chain), colour=sim)) + 
  geom_density() + facet_grid(mtn~effort, scales="free")
ggplot(filter(ggBeta, Parameter=="beta[2]"), 
       aes(x=value, group=paste(sim, Chain), colour=sim)) + 
  geom_density() + facet_grid(mtn~effort, scales="free")






pDet.full

sapply(1:20, FUN=function(x) mean(rbeta(20, 14.9, 19.5)))







RMSE.sum %>% group_by(effort) %>% 
  summarise(mn.lo=mean((mod.lo.rmse-obs.lo.rmse)/obs.lo.rmse)*100,
            mn.hi=mean((mod.hi.rmse-obs.hi.rmse)/obs.hi.rmse)*100) 

RMSE.sum %>% filter(b==50 & effort==1) %>% 
  group_by(mtn, era) %>% 
  summarise(obs.lo=mean((obs.lo.rmse)),
            obs.hi=mean((obs.hi.rmse)),
            mod.lo=mean((mod.lo.rmse)),
            mod.hi=mean((mod.hi.rmse))) %>%
  group_by(mtn) %>% 
  summarise(obs.lo=first(obs.lo)-last(obs.lo),
            mod.lo=first(mod.lo)-last(mod.lo),
            obs.hi=first(obs.hi)-last(obs.hi),
            mod.hi=first(mod.hi)-last(mod.hi))


for(mtn.i in c("FR", "SJ")) {
  for(era.i in c("H", "C")) {
    cat("-------------\n", mtn.i, era.i, "\n-------------\n")
    RMSE.sum %>% ungroup %>% mutate(b=as.character(b),
                                    effort=as.character(effort)) %>%
      filter(mtn==mtn.i & era==era.i) %>%
      aov((mod.lo.rmse-obs.lo.rmse)/obs.lo.rmse ~ effort + b, data=.) %>%
      TukeyHSD() %>% print
  }
}



for(mtn.i in c("FR", "SJ")) {
  for(era.i in c("H", "C")) {
    cat("-------------\n", mtn.i, era.i, "\n-------------\n")
    RMSE.sum %>% ungroup %>% mutate(b=as.character(b),
                                    effort=as.character(effort)) %>%
      filter(mtn==mtn.i & era==era.i) %>%
      aov((mod.hi.rmse-obs.hi.rmse)/obs.hi.rmse ~ effort + b, data=.) %>%
      TukeyHSD() %>% print
  }
}

