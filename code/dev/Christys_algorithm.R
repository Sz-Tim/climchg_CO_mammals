

# Algorithm for range based on sampling:
#  1. Species rarity metric: #/FR total, #/SJ total, #/rest total
#   - idea of rarity is good, but not sure about this implementation...
#   - definitely need some way to separate relative abundance from detectability
#  2. How often it is undetected in an elevation bin in its interpolated range
#   - using proportion of bins ignores differential sampling effort across bins
#   - this is a function of relative abundance and individual detectability
#  3. What is the detectability (# specimens of spp/# total specimens) across the range? [something better?]
#   - blurs relative abundance with detectability
#   - really need detectability on an individualized basis
#   - probably not possible to do with historic data 
#  4. What is the sampling effort (# total specimens) in the adjacent bins to the current range limits? 
#   - seems good
#  5. Lower probability the farther away from a known detection (e.g., 2 bins away lower prob than 1 bin etc)
#   - yes, but current implementation (weight by 2x, etc) is no good
#   - could use total probability including intermediate elevations
#   - could calculate sequentially, multiplying by calculated total
#     - effectively the same thing, just computationally cheaper
#   - could weight for each species, based on patchiness within interpolated range
#  6. Lower probability of detection outside absolute CO range limits
#   - is there any way to do this that isn't just arbitrary?
    
# Probability of occurring outside the detected elevational range 
# psi = pr(missed) = probability of occurring outside of the detected range


# Pr(missing) = pr(presentCOrng)  #??
#             + nUndetectedInterpEl/nInterpEl
#             + (0.5 - nSampleOutsideRng/total)  #??
#             - rarityRank/nSpp
#             - detectabilityMetric
#             - distanceFromDetectedRng

#    Pr(missing) = [Pr(presentCOrng)+Prop of bins missing in interpolated range+ (0.5 - #samples outside RL/total)] – [raritysp + Detectability metric (across range or near range edges)+ (#bins from detection)]



# Detectability vs. relative abundance


