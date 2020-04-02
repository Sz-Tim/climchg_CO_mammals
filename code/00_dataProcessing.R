library(fitdistrplus); library(tidyverse)



# aggregating occurrence data

f.path <- dir("data/orig/", "*ata.xlsx", full.names=T)

obs <- vector("list", length(f.path))
for(f in f.path) {
  f.sheets <- readxl::excel_sheets(f) %>% 
    grep(pattern="Sampling", value=T, invert=T) %>% 
    grep(pattern="UNK", value=T, invert=T)
  obs[[f]] <- map_dfr(f.sheets, 
                 ~readxl::read_xlsx(f, .x, range=readxl::cell_cols("A:B")) %>%
                   setNames(c("county", "el")) %>% mutate(sp=.x))
}
obs <- do.call('rbind', obs)

write_csv(obs, "data/sample_els_H.csv")




# probability of detection
# raw data: rows = species, columns = sites
prDet.raw <- readxl::read_xlsx("data/orig/Prob_Detection.xlsx", 1) %>%
  mutate(Genus=str_split_fixed(Species, " ", 2)[,1]) 
prDet.df <- prDet.raw %>%
  select(Order, Family, Genus, Species, 5:36) %>% 
  pivot_longer(5:36, names_to="Site", values_to="prDet", values_drop_na=T) 
# fit beta distribution for species with > 2 values
prDet.sum <- prDet.df %>%  group_by(Genus, Species) %>%
  summarise(tryGenus=n()<2, useFam=F, mn=mean(prDet), sd=sd(prDet), 
            shp1=ifelse(tryGenus, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[1]),
            shp2=ifelse(tryGenus, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[2]))
# fit beta distribution for genera with > 2 values
prDet.genus <- prDet.df %>%  group_by(Family, Genus) %>%
  summarise(useFam=n()<2, mn=mean(prDet), sd=sd(prDet),
            shp1=ifelse(useFam, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[1]),
            shp2=ifelse(useFam, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[2]))
# fit beta distribution for families with > 2 values (all)
prDet.fam <- prDet.df %>%  group_by(Family) %>%
  summarise(useOrder=n()<2, mn=mean(prDet), sd=sd(prDet),
            shp1=ifelse(useOrder, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[1]),
            shp2=ifelse(useOrder, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[2]))

# add all taxa, info
prDet.sp <- prDet.raw %>% select(Order, Family, Genus, Species) %>%
  left_join(., filter(prDet.sum, !is.na(shp1)), by=c("Species", "Genus")) %>%
  mutate(tryGenus=replace(tryGenus, is.na(tryGenus), TRUE))
prDet.genus <- prDet.raw %>% group_by(Genus, Family, Order) %>% summarise() %>%
  left_join(., prDet.genus %>% filter(!is.na(shp1)), 
            by=c("Family", "Genus")) %>%
  mutate(useFam=replace(useFam, is.na(useFam), TRUE))

# merge dataframes
prDet.full <- rbind(prDet.sp %>% filter(!tryGenus), # species values
                    left_join(prDet.sp %>% 
                                filter(tryGenus) %>% 
                                select(-useFam, -mn, -sd, -shp1, -shp2), 
                              rbind(prDet.genus %>%  # genus values 
                                      filter(!useFam) %>% 
                                      ungroup,
                                    left_join(prDet.genus %>% 
                                                filter(useFam) %>% 
                                                ungroup %>%
                                                select(-mn, -sd, -shp1, -shp2),
                                              prDet.fam %>% # family values
                                                select(-useOrder), 
                                              by="Family")
                              ),
                              by=c("Order", "Family", "Genus")
                    )
)

write_csv(prDet.full, "data/prDet_processed.csv")


cb.palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
for(f in 1:n_distinct(prDet.full$Family)) {
  fam <- unique(prDet.full$Family)[f]
  gen <- unique(filter(prDet.full, Family==fam & tryGenus)$Genus)
  gen <- setNames(seq_along(gen), gen)
  gray.palette <- paste0("gray", round(seq(30, 80, length.out=length(gen))))
  pDet.f <- filter(prDet.full, Family==fam) %>%
    mutate(lty=rep(c(1,5), each=length(cb.palette))[1:n()], 
           lty=ifelse(tryGenus, 2, lty),
           lty=ifelse(useFam, 3, lty),
           col=rep(cb.palette, 10)[1:n()], 
           col=ifelse(tryGenus, gray.palette[gen[Genus]], col),
           col=ifelse(useFam, "black", col),
           lwd=ifelse(tryGenus, 2, 1.5),
           lwd=ifelse(useFam, 1.5, lwd)) %>%
    arrange(Species)
  pdf(paste0("out/pDet/pDet_priors_", fam, ".pdf"), width=7, height=7)
  plot(NA, NA, xlim=c(0,1), ylim=c(0,20), xlab="pDet", ylab="Density", main=fam)
  walk(1:nrow(pDet.f), ~curve(dbeta(x, pDet.f$shp1[.], pDet.f$shp2[.]), 
                              from=0, to=1, add=T, col=pDet.f$col[.], 
                              lty=pDet.f$lty[.], lwd=pDet.f$lwd[.]))
  legend("topright", pDet.f$Species, col=pDet.f$col, lty=pDet.f$lty, bty="n", lwd=pDet.f$lwd)
  dev.off()
}
