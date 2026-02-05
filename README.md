Inhibitory Close Pairs (ICP) Geostatistical Sampling Implementation

Overview
This R script implements an Inhibitory Close Pairs (ICP) geostatistical sampling design based on the methodology described in BMJ Global Health (2022;7:e009500). The ICP approach allows for epidemiological surveys with similar precision to classical survey designs but with significantly smaller sample sizes, reduced travel distance, and shorter completion times.

 Prerequisites
 Required R Packages
```r
install.packages(c("geosphere", "ggplot2", "sp", "sf", "dplyr", "parallel", "data.table"))
```
 Required Data Files, Available Here
1. `kenya_villages.shp` - Village locations with longitude/latitude
2. `wards.shp` - Ward boundaries to attach villages to wards
3. `County.shp` - County boundaries for visualization

ICP Sampling Functions
 `select_primary_points()`
Selects primary sampling points with minimum distance constraints:
- Algorithm:
  1. Randomly selects candidate points
  2. Checks distance to existing points using Haversine formula
  3. Uses fallback strategy if no suitable point found
  4. Ensures geographical spread across study area

 `add_close_pairs()`
Adds close pairs to each primary point:

- Algorithm:
  1. For each primary point, finds villages within radius ζ
  2. Randomly selects one as close pair
  3. If none within radius, selects closest available village
  4. Prevents duplicate selection

Sampling Parameters
```r
 Default parameters (Adjust these based on study requirements)
n_primary <- 150       #Number of primary points
k <- 150               #Number of close pairs
delta_km <- 40         #δ (delta): Minimum distance between primary points. Larger values increase geographical spread but may make sampling more difficult.
zeta_km <- 15          #ζ (zeta): Radius for close pair selection. Smaller values capture local variation; larger values provide more flexibility in pair selection.
```

Some outputs
1. sampling_map.pdf: Visual representation of sampling distribution
2. ICP_sampled_villages.csv: Complete list of sampled villages with attributes

Advantages of ICP Design (Dileepa et al., 2022)
1. Reduced Sample Size: Achieves similar precision with 25% fewer samples than classical designs
2. Geographical Balance: Prevents oversampling in dense areas and undersampling in sparse areas
3. Efficient Travel: Reduces travel distance by >40% and survey time by >75% 
4. High Sensitivity/Specificity: >80% accuracy in detecting high-risk areas with sampling fractions >20%

Troubleshooting
1. No villages within radius: Script automatically selects closest village and provides warning

References

Ediriweera DS, de Silva T, Kasturiratne A, et al. Geographically regulated designs of incidence surveys can match the precision of classical survey designs whilst requiring smaller sample sizes: the case of snakebite envenoming in Sri Lanka. BMJ Global Health 2022;7:e009500.

```txt
pdf attached
