## landscape_spatial-temporal_ANPP
Analysis for spatial-temporal fluctuations in ANPP at CPER

## directory structure

# analyses - 
  
  R scripts to generated final analyses and figures
  
  01_site_map_figure.R: Code to generate site figure 1
  
  02_f2_precipitation_anpp_variation.R: Code to run analysis of spatial vs temporal variation in precipitation and ANPP and generate figure 2 
  
  03_f3_global_sensitivity.R: Code to run analysis of sensitivity of ANPP to growing season precipitation across functional groups and sites
  
  04_legacy_results.R: Code to variable selection analysis to identify imporant precipitation legacy variables
  
  05_f6_legacy_fig.R: Code to generate legacy variable figure 6
  
  06_supplemental_figures.R

# data - 
The main datasets used here are archived at Ag Data Commons (Merchant et al., 2026). To run these scripts download stored data and save into the data directory.

Citation: Merchant, T., Hajek, O. L., Augustine, D. J., Porensky, L. M., Kearney, S. P., & Hoover, D. L. (2026). Data from: Spatial and temporal drivers of landscape-level variation of aboveground net primary production in a semi-arid grassland [Dataset]. Ag Data Commons. https://doi.org/10.15482/USDA.ADC/32065638

Additional data used to generate the site figure
  
pastures.shp = shapefile of pastures at CPER
  id = pasture id
 

carm.shp = shapefile of CARM pastures at CPER
  pasture = pasture id

manual_spatialcoords.csv = historical locations of manual rain gauges at CPER
  pasture = pasture id
  UTM.E = easting coordinate
  UTM.N = northing coordinate

autogauge_spatialcoords.csv = locations of tipping bucket automated rain gauges at CPER
  pasture = pasture id
  UTM.E = easting coordinate
  UTM.N = northing coordinate
