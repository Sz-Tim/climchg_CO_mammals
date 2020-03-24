



########
## Setup
########

# load libraries
library(rjags); library(ggmcmc); library(tidyverse); theme_set(theme_bw())

# read in data
vole.df <- read_csv("data/by_taxon/voles_counts.csv") %>%
  gather(set, Detections, 3:11) %>%
  mutate(Detected=Detections>0,
         Time=str_split_fixed(set, "_", 2)[,1],
         Region=str_split_fixed(set, "_", 2)[,2])
vole.tot <- filter(vole.df, Species=="TOTAL")
vole.df <- filter(vole.df, Species!="TOTAL") %>% 
  group_by(Elev, Time, Region) %>%
  mutate(total=sum(Detections))
vole.rng <- vole.df %>% group_by(Species, Time, Region) %>%
  filter(Detected) %>%
  summarise(minDet=min(Elev),
            maxDet=max(Elev))
vole.df <- left_join(vole.df, vole.rng, by=c("Species", "Time", "Region")) %>%
  group_by(Species, Time, Region) %>%
  mutate(interpRange=(Elev<=maxDet & Elev >=minDet)) %>%
  ungroup

# set up parameters, storage objects
n.el <- n_distinct(vole.df$Elev)
n.spp <- n_distinct(vole.df$Species)
tr <- expand.grid(Time=c("A", "D"), Region=c("FR", "SJ", "CO"), 
                  stringsAsFactors=F)
Y <- binsAway <- map(1:nrow(tr), ~matrix(0, nrow=n.el, ncol=n.spp))
Ytot <- map(1:nrow(tr), ~rep(0, n.el))
interpPatchy <- spAbund <- map(1:nrow(tr), ~rep(0, n.spp))
delta <- runif(n.spp, 0, 0.5)

# calculate observations, binsAway, patchiness
for(i in 1:nrow(tr)) {
  voles.i <- vole.df %>% filter(Time==tr$Time[i] & Region==tr$Region[i])
  for(s in 1:n.spp) {
    sp.i <- voles.i %>% filter(Species==unique(Species)[s])
    Y[[i]][,s] <- sp.i$Detections
    interpPatchy[[i]][s] <- sum(sp.i$Detections==0 & 
                                  sp.i$interpRange)/sum(sp.i$interpRange)
    if(is.na(interpPatchy[[i]][s])) interpPatchy[[i]][s] <- 0
    rng.sp <- c(which(sp.i$Elev==sp.i$minDet[1]), which(sp.i$Elev==sp.i$maxDet[1]))
    for(j in which(!sp.i$interpRange)) {
      binsAway[[i]][j,s] <- ifelse(sp.i$Detected[j], 0, min(abs(j - rng.sp)))
    }
  }
  Ytot[[i]] <- rowSums(Y[[i]])
  spAbund[[i]] <- colSums(Y[[i]])
}



########
## Fit JAGS model
########
NZ.ls <- setNames(vector("list", nrow(tr)), paste(tr$Time, tr$Region, sep="_"))
for(i in 1:nrow(tr)) {
  # observations
  obs.df <- data.frame(el=rep(1:n.el, times=n.spp),
                       spp=rep(1:n.spp, each=n.el),
                       Y=c(Y[[i]]),
                       binsAway=c(binsAway[[i]]))
  
  # data to feed model
  jags_d <- list(n.spp=n.spp, 
                 n.el=n.el, 
                 Y=Y[[i]], 
                 Ytot=Ytot[[i]], 
                 delta=delta, 
                 spAbund=spAbund[[i]],
                 binsAway=binsAway[[i]], 
                 interpPatchy=interpPatchy[[i]])
  
  
  # parameters to store
  pars <- c("Z", "N", "b")
  
  # fit model
  mod <- jags.model(file="code/multinom_jags.txt", data=jags_d,
                    n.chains=3, n.adapt=5000, 
                    inits=list(Z=matrix(1, n.el, n.spp)))
  out <- coda.samples(mod, variable.names=pars, n.iter=10000, thin=20)
  
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
           spp=as.numeric(spp))
  gg.NZ.sum <- gg.NZ %>% group_by(el, spp) %>%
    summarise(mn_prPres=mean(NZ>0),
              prPres_gt05=as.numeric(mn_prPres>0.05)) %>%
    full_join(obs.df, c("el", "spp")) %>%
    mutate(Time=tr$Time[i],
           Region=tr$Region[i])
  NZ.ls[[i]] <- gg.NZ.sum
}

all.df <- do.call("rbind", NZ.ls) %>%
  mutate(spp=factor(spp, levels=1:n.spp, labels=unique(vole.df$Species)))
time.diff <- all.df %>% select(el, spp, Region, Time, mn_prPres) %>% 
  spread(Time, mn_prPres)

ggplot(all.df, aes(x=el, y=mn_prPres, colour=Time, fill=Time)) + 
  geom_line() + 
  geom_area(data=filter(all.df, Time=="A"), aes(y=as.numeric(Y>0)), 
            alpha=0.25, colour=NA) + 
  geom_area(data=filter(all.df, Time=="D"), aes(y=as.numeric(Y>0)), 
            alpha=0.25, colour=NA) + 
  facet_grid(Region~spp) + ylim(0,1) +
  labs(x="Elevational bin", y="Predicted probability of presence")

ggplot(all.df, aes(x=el, y=as.numeric(mn_prPres>0.05), colour=Time, fill=Time)) + 
  geom_line() + 
  geom_area(data=filter(all.df, Time=="A"), aes(y=as.numeric(Y>0)), 
            alpha=0.25, colour=NA) + 
  geom_area(data=filter(all.df, Time=="D"), aes(y=as.numeric(Y>0)), 
            alpha=0.25, colour=NA) + 
  facet_grid(Region~spp) + ylim(0,1) +
  labs(x="Elevational bin", y="Predicted probability of presence")

ggplot(time.diff, aes(x=el, y=D-A)) + geom_line() +
  facet_grid(Region~spp) + ylim(-1,1) +
  labs(x="Elevational bin", y="Current - historic")

ggplot(time.diff, aes(x=A, y=D, colour=el)) + geom_point() +
  facet_grid(Region~spp) + ylim(-1,1) +
  labs(x="Current", y="Historic")
