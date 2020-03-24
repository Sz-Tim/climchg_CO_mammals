library(tidyverse); theme_set(theme_bw())
vole.df <- read_csv("data/voles_counts.csv") %>%
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
  mutate(interpRange=(Elev<=maxDet & Elev >=minDet))

ggplot(vole.tot, aes(x=Elev, y=Detections, colour=Time)) + geom_line() + 
  facet_wrap(~Region)



# Detection = nDetect / nSamplesInInterpolatedRange
vole.det.all <- vole.df %>% ungroup %>% 
  filter(interpRange) %>% group_by(Species, Region, Time, Elev) %>%
  summarise(propSamples=sum(Detections)/sum(total))
vole.det.reg <- vole.df %>% ungroup %>% 
  filter(interpRange) %>% group_by(Species, Region, Elev) %>%
  summarise(propSamples=sum(Detections)/sum(total))
vole.det.reg.avg <- vole.df %>% ungroup %>%
  filter(interpRange) %>% group_by(Species, Time, Region) %>%
  summarise(propSamples=sum(Detections)/sum(total),
            propPresence=mean(Detected))
vole.det.el <- vole.df %>% ungroup %>% 
  filter(interpRange) %>% group_by(Species, Elev) %>%
  summarise(propSamples=sum(Detections)/sum(total))
vole.det <- vole.df %>% ungroup %>% 
  filter(interpRange) %>% group_by(Species) %>%
  summarise(propSamples=sum(Detections)/sum(total))

vole.det.reg.avg %>% filter(Species=="Mgap" & Time!="BC" & Region=="FR")

ggplot(filter(vole.det.reg.avg, Time!="BC"),
       aes(Time, propPresence, colour=Region, group=Region)) + geom_point() +
  facet_wrap(~Species) + ylim(0,1) + geom_line()

ggplot(vole.det.all, aes(x=Elev, y=propSamples, colour=Time)) + ylim(0,1) +
  stat_smooth(method="loess", span=2, se=F) + geom_point() + 
  facet_grid(Region~Species) +
  labs(title="nDetect/nSamples within interpolated ranges", 
       x="Elevation", y="Proportion of Samples")
ggplot(vole.det.reg, aes(x=Elev, y=propSamples)) + 
  stat_smooth(method="loess", span=2, se=F) + geom_point() + ylim(0,1) +
  facet_grid(Region~Species) +
  labs(title="nDetect/nSamples within interpolated ranges", 
       x="Elevation", y="Proportion of Samples")
ggplot(vole.det.el, aes(x=Elev, y=propSamples)) + 
  stat_smooth(method="loess", span=2, se=F) + geom_point() + ylim(0,1) +
  facet_wrap(~Species) +
  labs(title="nDetect/nSamples within interpolated ranges", 
       x="Elevation", y="Proportion of Samples")


ggplot(vole.df, aes(Elev, Detections, colour=Time)) + geom_line() +
  facet_grid(Region~Species)
ggplot(vole.df, aes(Elev, Detections/total, fill=Species)) + geom_area() +
  facet_grid(Region~Time)

ggplot(vole.df, aes(Elev, Detected, colour=Time)) + 
  geom_jitter(width=0, height=0.1, shape=1) +
  facet_grid(Region~Species)

ggplot(filter(vole.rng, Time != "BC")) + 
  geom_point(aes(x=minDet, y=Time)) + geom_point(aes(x=maxDet, y=Time)) + 
  geom_segment(aes(x=minDet, xend=maxDet, y=Time, yend=Time)) + 
  facet_grid(Region~Species)

ggplot(filter(vole.df, Time!="BC" & interpRange)) + 
  geom_point(aes(x=Elev, y=Time, colour=Detected)) + facet_grid(Region~Species)


vole.df %>% filter(Time!="BC") %>% 
  summarise(nMissing=sum(interpRange & !Detected))

vole.df %>% filter(Time!="BC") %>% 
  summarise(nMissing=sum(interpRange & !Detected),
            propMissing=nMissing/sum(interpRange)) %>% 
  ggplot(aes(x=Region, y=propMissing, colour=Time)) + 
  geom_point() + facet_wrap(~Species)
