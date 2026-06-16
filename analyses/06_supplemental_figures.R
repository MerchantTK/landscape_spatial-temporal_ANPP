#library
library(tidyr)
library(dplyr)
library(ggplot2)
library(sp)
library(stars)
library(raster)
library(patchwork)
#palettes
library(rcartocolor)

library(ggmap)
library(maps)
library(mapdata)
library(sfheaders)
library(lme4)
require(sjPlot)

precip.palette <- c(  '#F78154','#389DE5', '#5DDCAC' ,'#9368B7')
#set up shape files 

past.df <- read.csv('cper_broom_shp_v2.csv')

past.sf <- past.df %>%
  st_as_sf(coords = c("long", "lat")) %>%
  group_by(id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

#Supplemental site fig 1
ppt.1218 <- read.csv('ppt_full_12-18.csv')
can.coords <- read.csv('Spatialcoords.csv')%>%
  dplyr::select(Pasture, UTM.E, UTM.N)%>%
  #rename_with(tolower)%>%
  st_as_sf(coords = c('UTM.E', 'UTM.N'))

inst.coords <- read.csv('CPER_Instrumentation_Locations.csv')%>%
  st_as_sf(coords = c('Easting.UTM', 'Northing.UTM'))


blockpairs <- read.csv('plotpairs.csv')%>%
  dplyr::select(pasture, block)%>%
  mutate(pasture = ifelse(pasture == 'NH', '10S',pasture))%>%
  unique()
biomass <- read.csv('AGM_Biomass_Widecln_attr_2024-03-26.csv')%>%
  rename_all(tolower)%>%
  dplyr::select('pasture', 'x', 'y')%>%
  st_as_sf(coords = c('x', 'y'))

unique(past.sf$id)
past.blocks.sf <- past.sf%>%
  rename( pasture = id)%>%
  left_join(blockpairs, by = 'pasture')

ggplot()+
  geom_sf(data = past.blocks.sf, fill = 'grey40')+
  geom_sf(data = past.blocks.sf, aes(fill = as.factor(block)))+
  #geom_sf_label(data = past.blocks.sf, aes(label = pasture))
  geom_sf(data = biomass, color = 'blue3', size = 2)+
  geom_sf(data = can.coords, color = 'red3', size = 2)+
  geom_sf(data = inst.coords, color = 'green', size = 2)+
  #scale_fill_viridis_d(option = 'turbo')+
  labs(fill = 'Block')+
  theme_void(base_size = 12)+
  ggspatial::annotation_scale(location = 'bl')+
  ggspatial::annotation_north_arrow(pad_y = unit(1, "cm"))
?ggspatial::annotation_north_arrow


## Breakdown of biomass percentage by ecosite 
#Supplemental figure 3
#full dataframe
veg.precip.df <- read.csv('data/legacy_full_df.csv')

veg.precip.df%>%
  dplyr::select(pairblock, ecosite, water.yr, wspg, c3pg, forb)%>%
  pivot_longer(-c(pairblock, ecosite, water.yr), values_to = 'mass', names_to = 'fun.group')%>%
  ggplot()+
  stat_summary(fun = mean , geom = 'col', aes(x = water.yr, y = mass, fill = fun.group), position = position_stack())+
  facet_wrap(~ecosite, ncol = 1)+
  scale_fill_manual(values = precip.palette, labels = c('Cool', 'Forb', 'Warm'))+
  theme_bw()+
  labs(fill = 'Functional Group', y = bquote('ANPP kg ha'^-1*'mm'^-1), x = 'Harvest year')

