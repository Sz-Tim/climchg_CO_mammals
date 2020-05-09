Files:

####----------------------------------------------------------------------------
1. species_predictions_50m_all_spp.csv

Summary of elevational distribution for each species

Abbrev: species abbreviation
mtn: mountain range
period: time period / era
Source: obs (observed) modPres (posterior mean > 0.05, meaning <95% prAbsence)
elMin: minimum elevation observed or predicted
elMax: maximum elevation observed or predicted
elMid: median of elMax and elMin
minChg: change in elMin across periods, by source and mtn (negative: C < H)
maxChg: change in elMax across periods, by source and mtn (negative: C < H)
midChg: change in elMid across periods, by source and mtn (negative: C < H)
minDir: direction of minChg, with threshold of ±100m (exactly 100 = change)
minDir: direction of minChg, with threshold of ±100m (exactly 100 = change)
minDir: direction of minChg, with threshold of ±100m (exactly 100 = change)




####----------------------------------------------------------------------------
2. species_predictions_50m_spp_10obs.csv

Summary of elevational distribution for each species with at least 10 observations in the historic dataset for either SJ or FR

Abbrev: species abbreviation
mtn: mountain range
period: time period / era
Source: obs (observed) modPres (posterior mean > 0.05, meaning <95% prAbsence)
elMin: minimum elevation observed or predicted
elMax: maximum elevation observed or predicted
elMid: median of elMax and elMin
minChg: change in elMin across periods, by source and mtn (negative: C < H)
maxChg: change in elMax across periods, by source and mtn (negative: C < H)
midChg: change in elMid across periods, by source and mtn (negative: C < H)
minDir: direction of minChg, with threshold of ±100m (exactly 100 = change)
minDir: direction of minChg, with threshold of ±100m (exactly 100 = change)
minDir: direction of minChg, with threshold of ±100m (exactly 100 = change)




####----------------------------------------------------------------------------
3. all_50m_bins_ALL_SPECIES_NEW.csv

Full predictions for each species x elevational bin x mtn x period

bin: elevational bin index
spp: species index
mn_N: posterior mean for lambda (relative abundance)
mn_prPres: posterior mean for presence = pr(lambda * Z) > 0
Elevation: elevation
Abbrev: species abbreviation
Y: observations of each species in each bin (by mtn x period)
Ytot: total observations in bin (by mtn x period)
binsAway: number of bins outside interpolated range; 0 = inside observed range
distAway: binsAway * binsize (i.e., times 50)
distAway_sc: scale(distAway) -- standardized for use in model
interpPatchy: nBinsObs / nBinsInterpolated
spAbund: total number of observations for species (by mtn x period)
delta: mean individual detection probability from species' empirical distribution
period: time period / era
mtn: mountain range
binSize: elevational bin size (50m)
set: mtn_binSize_period




####----------------------------------------------------------------------------
4. prDet_processed.csv

Individual detection probabilities for each species. If there were at least 2 observations for a species, a Beta distribution was fit using the R function fitdistrplus::fitdist(). For species with < 2 values, the whole genus was pooled. For genera with < 2 values, the whole family was pooled. The parameters for the Beta distribution serve as the prior distribution for each species in the full model. 

Order
Family
Genus
Species
tryGenus: TRUE if nObs per species < 2; use values for whole genus
useFam: TRUE if nObs per genus < 2; use values for whole family
mn: mean value
sd: standard deviation 
shp1: alpha = first shape parameter for Beta distribution
shp2: beta = second shape parameter for Beta distribution




####----------------------------------------------------------------------------
5. sample_els_H.csv

Historic elevations across all species.

County
el
sp




####----------------------------------------------------------------------------
6. sample_els_C.csv

Contemporary elevations across all species.

County
el
sp




####----------------------------------------------------------------------------
7. simulations_out.csv

Output from sensitivity analysis runs. Results for 50 simulated samples and model runs for communities on each mountain range and era, with sampling effort at 0.5, 1, or 2 x the true effort (preserving the relative sampling distribution across elevations) and elevational bin size at 50m, 100m, or 200m. Communities were simulated, so RMSE is compared to the true range boundaries. 

spp: simulated species number
true.lo: true lower boundary
true.hi: true upper boundary
obs.lo: sampled lower boundary
obs.hi: sampled upper boundary
mod.lo: modeled lower boundary
mod.hi: modeled upper boundary
sim: simulation number (1-50)
b: elevational bin size (50, 100, 200)
mtn: mountain range (FR, SJ)
type: "all", meaning all species were included
effort: proportional effort compared to empirical sampling distribution (.5, 1, 2)
era: sampling era (H: history, C: contemporary)




####----------------------------------------------------------------------------
8. simulations_RMSE.csv

Root mean squared error from sensitivity analysis runs, calculated from 'simulations_out.csv'. 

b: elevational bin size (50, 100, 200)
mtn: mountain range (FR, SJ)
type: "all", meaning all species were included
effort: proportional effort compared to empirical sampling distribution (.5, 1, 2)
era: sampling era (H: history, C: contemporary)
sim: simulation number (1-50)
obs.lo.rmse: sqrt( mean( (obs.lo-true.lo)^2 ) )
mod.lo.rmse: sqrt( mean( (mod.lo-true.lo)^2 ) )
obs.hi.rmse: sqrt( mean( (obs.hi-true.hi)^2 ) )
mod.hi.rmse: sqrt( mean( (mod.hi-true.hi)^2 ) )


