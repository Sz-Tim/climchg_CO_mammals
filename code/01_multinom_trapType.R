






# This script fits a JAGS model where the number of individuals of each species
# detected at each elevation is drawn from a multinomial distribution, where the
# probability for each species is weighted by: 
#   1) the overall number of detections
#   2) the individual-level detection probability
#   3) the number of elevational bins away from its nearest detection
#   4) the patchiness of its detections within its interpolated range.
#
# More specifically, at each elevation i, Y[i] individuals are drawn from the
# species pool, where each species j is drawn with probability:
# lambda[i,j]*delta[j] where lambda[i,j] is the latent relative abundance of the
# species, and delta[j] is the probability that a given individual of the
# species is detected when present.
#
# Thus, at each elevation, we observe draws from the species pool, where the
# probability that any given individual belongs to species j is weighted by the
# species' relative abundance and the individual-level detection probability.
# That is,
#
# y[i,1:J] ~ multinomial(lambda[i,1:J]*delta[1:J], Y[i])
#
# where lambda[i,j] is the true abundance of each species j at elevation i,
# delta[j] is the individual-level detection probability for each species j, 
# Y[i] is the total number of detections at elevation i, and J is the total
# number of species potentially occurring at elevation i (i.e., the number of 
# species detected in the mountain range that could be detected using the).
# 
# Z[i,j] is drawn from a bernoulli distribution, where pr(Z[i,j]==1) = psi[i,j].
# psi[i,j] = antilogit(a[s] + 
#                      b[1]*distAway[i,j] + 
#                      b[2]*interpPatchiness[j] 
# 
# The model for the sampling-based probability of presence therefore includes:
#   - distance in meters from interpolated range (outside range only)
#   - patchiness within interpolated range (bins present / total bins)
#
# The variables are: 
# n.el: number of elevational bins
# J: number of species
# delta: individual-level detection probability, length=J
# a: intercept (species-specific)
# b: slopes (global)
# LAMBDA: total detections along gradient, length=J
# Z: latent presence/absence based on sampling only, nrows=n.el, ncol=J
# N: latent relative abundances, nrows=n.el, ncol=J
# y: detections, nrows=n.el, ncol=J
# Y: total detections at each elevation, length=n.el
# binsAway: number of bins away from detected range, nrows=n.el, ncol=J
# distAway: number of meters away from detected range, nrows=n.el, ncol=J
# interpPatchy: proportion non-detections in interpolated range, length=J
# rng: true range boundaries, nrow=2, ncol=n.spp





########
## Setup
########
#--- load libraries
library(rjags); library(dclone)
library(ggmcmc); library(tidyverse); theme_set(theme_bw())
all_species <- T

#--- time period, elevational bin size, mtn, and trap type combinations
if(all_species) {
  sets <- expand.grid(period=c("H", "C"),
                      bin=50, 
                      mtn=c("FR", "SJ"))
} else {
  sets <- expand.grid(period=c("H", "C"),
                      bin=50, 
                      mtn=c("FR", "SJ"), 
                      trap=c("Shrew", "Sherman", "Large")) 
}
sets.out <- a.out <- b.out <- vector("list", nrow(sets))

for(set in 1:nrow(sets)) {
  p <- sets$period[set]
  b <- sets$bin[set]  # elevational bin size
  m <- sets$mtn[set]  # mountain range
  if(all_species) {
    tr <- c("Shrew", "Sherman", "Large")
  } else {
    tr <- sets$trap[set]  # trap type
  }
  #--- read in data
  # species information: filter by trapping method
  sp.i <- read_csv("data/by_trapType/Mamm_Summary_Data.csv") %>%
    filter(Sample %in% tr)
  pDet <- read_csv("data/prDet_processed.csv") %>%
    mutate(Abbrev=sp.i$Abbrev[match(Species, sp.i$Species)])
  # historic data: rename columns, select species
  Y.df <- read_csv(paste0("data/by_trapType/Mamm_", m, "_", b, "m.csv")) %>%
    select("Elevation", contains(paste0("_", p))) %>%
    setNames(., str_remove(names(.), paste0("_", p))) %>%
    select("Elevation", one_of(sp.i$Abbrev[sp.i$Sample %in% tr]))
  # interpolated ranges
  Y.rng <- map_dfc(2:ncol(Y.df), ~range(which(Y.df[,.]>0))) %>% as.matrix
  
  #--- reshape data for JAGS
  Y <- as.matrix(Y.df[,-1])  # remove Elevation column
  n.el <- nrow(Y)
  n.spp <- ncol(Y)
  Ytot <- rowSums(Y)
  spAbund <- colSums(Y)
  delta_shp <- as.matrix(pDet[match(names(spAbund), pDet$Abbrev),9:10])
  interpPatchy <- rep(0, n.spp)
  binsAway <- interpRng <- matrix(0, nrow=n.el, ncol=n.spp)
  # calculate interpPatchy & binsAway for each species
  for(s in 1:n.spp) {
    if(is.infinite(Y.rng[1,s])) {
      interpPatchy[s] <- 1
      binsAway[,s] <- n.el
    } else {
      interpRng.s <- Y.rng[1,s]:Y.rng[2,s]
      interpRng[,s] <- 1:n.el %in% interpRng.s
      interpPatchy[s] <- sum(Y[interpRng.s,s]==0)/length(interpRng.s)
      for(j in (1:n.el)[-interpRng.s]) {
        binsAway[j,s] <- min(abs(j - Y.rng[,s]))
      }
    }
  }
  
  # store data for JAGS in list
  jags_d <- list(J=n.spp, 
                 n.el=n.el, 
                 y=Y, 
                 Y=Ytot, 
                 delta_shp=delta_shp, 
                 LAMBDA=spAbund,
                 interpPatchy=c(scale(interpPatchy)),
                 distAway=matrix(scale(c(binsAway*b)), ncol=n.spp))
  # observed data
  obs.df <- data.frame(bin=rep(1:n.el, times=n.spp),
                       Elevation=rep(Y.df$Elevation, times=n.spp),
                       spp=rep(1:n.spp, each=n.el),
                       Abbrev=rep(names(spAbund), each=n.el),
                       Y=c(Y),
                       Ytot=rep(Ytot, times=n.spp),
                       binsAway=c(binsAway),
                       distAway=c(binsAway)*b,
                       distAway_sc=c(scale(c(binsAway)*b)),
                       interpPatchy=rep(interpPatchy, each=n.el),
                       spAbund=rep(spAbund, each=n.el),
                       delta=rep(delta_shp[,1]/rowSums(delta_shp), each=n.el))
  
  
  
  ########
  ## Run JAGS model
  ########
  
  mod <- c("simple", "full")[2]
  #--- parameters to store
  if(mod=="simple") pars <- c("Z", "lambda", "beta", "alpha", "sigma") 
  if(mod=="full") pars <- c("lambda", "Z", "beta", "a", "sigma")
  
  #--- fit model
  cl <- makeCluster(4)
  out <- jags.parfit(cl=cl, data=jags_d, params=pars, 
                     model=ifelse(mod=="simple", 
                                  "code/dev/multinom_b_global_simple.txt",
                                  "code/00_multinom_b_global.txt"),
                     inits=list(Z=matrix(1, n.el, n.spp)),
                     n.chains=4, n.adapt=5000, n.update=1000, 
                     n.iter=1000, n.thin=5)
  stopCluster(cl)
  
  
  
  
  ########
  ## Summarize output
  ########
  gg.a <- ggs(out, "a")
  gg.b <- ggs(out, "beta")
  gg.Z <- ggs(out, "Z") %>%
    mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
             str_remove("Z\\[") %>% as.numeric,
           spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
             str_remove("\\]") %>% as.numeric) %>%
    rename(Z=value)
  if(mod=="simple") {
    gg.NZ <- gg.Z %>%
      full_join(obs.df, c("bin", "spp")) %>%
      mutate(period=p,
             mtn=m,
             binSize=b,
             set=paste0(m, "_", b, "m", "_", p),
             spp=as.numeric(spp))
    gg.NZ.sum <- gg.NZ %>% group_by(bin, spp) %>%
      summarise(mn_prPres=mean(Z>0)) %>%
      full_join(obs.df, c("bin", "spp")) %>%
      mutate(period=p,
             mtn=m,
             binSize=b,
             set=paste0(m, "_", b, "m", "_", p))
  } else {
    gg.N <- ggs(out, "lambda") %>%
      mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
               str_remove("lambda\\[") %>% as.numeric,
             spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
               str_remove("\\]") %>% as.numeric) %>%
      rename(N=value)
    gg.N.sum <- gg.N %>% group_by(bin, spp) %>%
      summarise(mn_N=mean(N),
                mn_prPres=mean(N>0)) %>%
      full_join(obs.df, c("bin", "spp")) %>%
      mutate(period=p,
             mtn=m,
             binSize=b,
             set=paste0(m, "_", b, "m", "_", p))
    # gg.NZ <- gg.N %>% mutate(NZ=N*Z>0, 
    # spp=as.numeric(spp))
    gg.NZ <- full_join(gg.Z, gg.N, c("bin", "spp", "Chain", "Iteration")) %>%
      mutate(NZ=N*Z,
             spp=as.numeric(spp))
    gg.NZ.sum <- gg.NZ %>% group_by(bin, spp) %>%
      summarise(mn_N=mean(N),
                mn_prPres=mean(NZ>0)) %>%
      full_join(obs.df, c("bin", "spp")) %>%
      mutate(period=p,
             mtn=m,
             binSize=b,
             set=paste0(m, "_", b, "m", "_", p))
  }
  
  
  sets.out[[set]] <- gg.NZ.sum
  b.out[[set]] <- gg.b
  a.out[[set]] <- gg.a
}

out.all <- do.call("rbind", sets.out) %>% 
  arrange(Abbrev, Elevation, period, mtn)
write_csv(out.all, paste0("out/all_", sets$bin[1], "m_bins_ALL_SPECIES.csv"))


########
## Plot output
########
#--- observations
FR_obs_H <- read_csv("data/by_trapType/Mamm_FR_RawData.csv") %>%
  select(contains("_H")) %>% setNames(., str_remove(names(.), "_H")) %>%
  gather(Abbrev, Elevation) %>% filter(!is.na(Elevation)) %>% unique
FR_obs_C <- read_csv("data/by_trapType/Mamm_FR_RawData.csv") %>%
  select(contains("_C")) %>% setNames(., str_remove(names(.), "_C")) %>%
  gather(Abbrev, Elevation) %>% filter(!is.na(Elevation)) %>% unique
SJ_obs_H <- read_csv("data/by_trapType/Mamm_SJ_RawData.csv") %>%
  select(contains("_H")) %>% setNames(., str_remove(names(.), "_H")) %>%
  gather(Abbrev, Elevation) %>% filter(!is.na(Elevation)) %>% unique
SJ_obs_C <- read_csv("data/by_trapType/Mamm_SJ_RawData.csv") %>%
  select(contains("_C")) %>% setNames(., str_remove(names(.), "_C")) %>%
  gather(Abbrev, Elevation) %>% filter(!is.na(Elevation)) %>% unique

#--- model output
out.all <- read_csv(paste0("out/all_", sets$bin[1], "m_bins_ALL_SPECIES.csv"))

#--- probability of presence
ggplot(filter(out.all, mtn=="FR"), 
       aes(x=Elevation)) + ylim(0,1) + facet_wrap(~Abbrev) +
  geom_point(data=FR_obs_H, aes(y=1), alpha=0.4, shape=1, colour="#7b3294") +
  geom_point(data=FR_obs_C, aes(y=0.9), alpha=0.4, shape=1, colour="#008837") +
  geom_line(aes(y=mn_prPres, group=period, colour=period)) +
  scale_colour_manual("Time\nPeriod", values=c("H"="#7b3294", "C"="#008837")) +
  labs(x="Elevation", y="Predicted probability of presence")
ggsave(paste0("figs/FR_prP_", sets$bin[1], "m_ALL_SPECIES.pdf"), 
       width=12, height=10, dpi=300)
ggplot(filter(out.all, mtn=="SJ"), 
       aes(x=Elevation)) + ylim(0,1) + facet_wrap(~Abbrev) +
  geom_point(data=SJ_obs_H, aes(y=1), alpha=0.4, shape=1, colour="#7b3294") +
  geom_point(data=SJ_obs_C, aes(y=0.9), alpha=0.4, shape=1, colour="#008837") +
  geom_line(aes(y=mn_prPres, group=period, colour=period)) +
  scale_colour_manual("Time\nPeriod", values=c("H"="#7b3294", "C"="#008837")) +
  labs(x="Elevation", y="Predicted probability of presence")
ggsave(paste0("figs/SJ_prP_", sets$bin[1], "m_ALL_SPECIES.pdf"), 
       width=12, height=10, dpi=300)

#--- probability of presence ≥ 0.05
ggplot(filter(out.all, mtn=="FR"), 
       aes(x=Elevation)) + ylim(0,1) + facet_wrap(~Abbrev) +
  geom_point(data=FR_obs_H, aes(y=1), alpha=0.4, shape=1, colour="#7b3294") +
  geom_point(data=FR_obs_C, aes(y=0.9), alpha=0.4, shape=1, colour="#008837") +
  geom_line(aes(y=as.numeric(mn_prPres>0.05), group=period, colour=period)) +
  scale_colour_manual("Time\nPeriod", values=c("H"="#7b3294", "C"="#008837")) +
  labs(x="Elevation", y=">5% Predicted probability of presence")
ggsave(paste0("figs/FR_pr_gr05_", sets$bin[1], "m_ALL_SPECIES.pdf"), 
       width=12, height=10, dpi=300)
ggplot(filter(out.all, mtn=="SJ"), 
       aes(x=Elevation)) + ylim(0,1) + facet_wrap(~Abbrev) +
  geom_point(data=SJ_obs_H, aes(y=1), alpha=0.4, shape=1, colour="#7b3294") +
  geom_point(data=SJ_obs_C, aes(y=0.9), alpha=0.4, shape=1, colour="#008837") +
  geom_line(aes(y=as.numeric(mn_prPres>0.05), group=period, colour=period)) +
  scale_colour_manual("Time\nPeriod", values=c("H"="#7b3294", "C"="#008837")) +
  labs(x="Elevation", y=">5% Predicted probability of presence")
ggsave(paste0("figs/SJ_pr_gr05_", sets$bin[1], "m_ALL_SPECIES.pdf"), 
       width=12, height=10, dpi=300)

#--- effect of distance & detectability
ggplot(filter(out.all, binsAway>0), aes(distAway, mn_prPres, colour=set)) + 
  geom_point(alpha=0.5, shape=1) + 
  scale_x_continuous("Distance from detection", breaks=c(0, 1000, 2000)) + 
  scale_y_continuous("Probability of presence", breaks=c(0,0.5,1), limits=c(0,1)) +
  facet_wrap(~delta*Abbrev) + stat_smooth(se=F, method="loess")
ggsave(paste0("figs/distanceByDelta", sets$bin[1], "m_ALL_SPECIES.pdf"), 
       width=12, height=12, dpi=300)


dist.ctr <- attr(scale(c(binsAway*b)), "scaled:center")
dist.scl <- attr(scale(c(binsAway*b)), "scaled:scale")
intP.ctr <- attr(scale(c(interpPatchy)), "scaled:center")
intP.scl <- attr(scale(c(interpPatchy)), "scaled:scale")

logit <- function(p) {log(p/(1-p))}
antilogit <- function(x) {exp(x)/(1+exp(x))}

dist.seq <- seq(min(jags_d$distAway), (1000-dist.ctr)/dist.scl, length.out=100)
intP.seq <- seq(min(jags_d$interpPatchy), max(jags_d$interpPatchy), length.out=100)

a.df <- filter(gg.a, Parameter=="a[37]")
b1.df <- filter(gg.b, Parameter=="beta[1]")
b2.df <- filter(gg.b, Parameter=="beta[2]")
b3.df <- filter(gg.b, Parameter=="beta[3]")
plot(NA, NA, xlim=c(0, max(dist.seq)*dist.scl+dist.ctr), ylim=c(0,1), 
     xlab="Distance from interpolated range", 
     ylab="Relative probability of presence",
     main="Distance effect | interpPatchy")
for(i in 1:nrow(b1.df)) {
  lines(dist.seq*dist.scl + dist.ctr, antilogit(dist.seq*b1.df$value[i]),
        col=rgb(0,0,0,0.01))
}


plot(NA, NA, xlim=c(0, max(obs.df$interpPatchy)), ylim=c(0,1), 
     xlab="Patchiness of interpolated range", 
     ylab="Relative probability of presence",
     main="Patchiness effect | distance")
for(i in 1:nrow(b2.df)) {
  lines(intP.seq*intP.scl + intP.ctr, antilogit(intP.seq*b2.df$value[i]),
        col=rgb(0,0,0,0.01))
}


plot(NA, NA, xlim=c(0, max(dist.seq)*dist.scl+dist.ctr), ylim=c(0,1), 
     xlab="Distance from interpolated range", 
     ylab="Relative probability of presence",
     main="Patchiness:Distance interaction | delta")
for(i in 1:length(intP.int)) {
  for(j in 1:nrow(b3.df)) {
    lines(dist.seq*dist.scl + dist.ctr, 
          antilogit(dist.seq*intP.int[i]*b3.df$value[j]), 
          col=rgb(i==1,i==2,i==3,0.05))
  }
}

s <- 8
plot(NA, NA, xlim=c(0, max(dist.seq)*dist.scl+dist.ctr), ylim=c(0,1), 
     xlab="Distance from interpolated range", 
     ylab="Relative probability of presence",
     main="Detection:Distance interaction | patchiness")
for(i in 1:length(delt.int)) {
  for(j in 1:nrow(b2.df)) {
    lines(dist.seq*dist.scl + dist.ctr, 
          antilogit(filter(gg.a, Parameter=="a[37]")$value[j] + 
                      dist.seq*b1.df$value[j] + 
                      interpPatchy[s]*b2.df$value[j] +
                      logit(delt.int[i])*b3.df$value[j] +
                      dist.seq*logit(delt.int[i])*b4.df$value[j]), 
          col=rgb(i==1,i==2,i==3,0.05))
  }
}
legend("topright", title="Detection", col=c("red", "green", "blue"), lwd=1,
       legend=round(delt.int, 3))



plot(NA, NA, xlim=c(0, max(dist.seq)*dist.scl+dist.ctr), ylim=c(0,1), 
     xlab="Distance from interpolated range", 
     ylab="Relative probability of presence",
     main="Patchiness & distance effect | delta")
for(i in 1:length(intP.int)) {
  for(j in 1:nrow(b3.df)) {
    lines(dist.seq*dist.scl + dist.ctr, 
          antilogit(a.df$value[j] + dist.seq*b1.df$value[j] +
                      intP.int[i]*b2.df$value[j] +
                      dist.seq*intP.int[i]*b3.df$value[j]), 
          col=rgb(i==1,i==2,i==3,0.05))
  }
}

