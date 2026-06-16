####generate the site figure 1
#TKM


`%notin%` <- Negate(`%in%`) 
#library
library(tidyverse)
library(sf)
#library(ggmap)
library(maps)
library(mapdata)
library(spData)
library(sfheaders)
library(egg)
require(ggpmisc)
library(mapme.biodiversity)
install.packages('speciesgeocodeR')

precip.palette <- c(  '#F78154','#389DE5', '#5DDCAC' ,'#9368B7')
#set up shape files 

past.sf <- st_read('data/pastures.shp')



##Site figure
carm.sf <- st_read('data/CARM_data/CARM_data.shp')

st_crs(past.sf) <- st_crs(carm.sf)


#State geometries 
data(us_states)
us_states

central_plains_states <- us_states%>%
  filter(NAME %in%  c('Colorado', 'Wyoming', 'Nebraska', 'Kansas', 'Oklahoma', 'New Mexico', 'Texas'))%>%
  st_transform(crs = st_crs(carm.sf))
# st_as_sf( coords = c('long', 'lat'))%>%
#group_by(region)%>%
# summarise(geometry = st_combine(geometry)) %>%
#st_cast("POLYGON")



##to add ecoregions download from https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip
# eco_regions <- st_read('data/ecoregions_wwf/wwf_terr_ecos.shp')%>%
#   filter(REALM == 'NA')
# short_grass <- eco_regions%>%
#   filter(ECO_NAME %in% c('Western short grasslands'))



carm.point.sf <- st_transform(carm.point.sf, crs = st_crs(central_plains_states))



mapreference <- ggplot()+
  geom_sf(data = central_plains_states, fill = 'white', alpha = 0.05, linewidth = 1)+
  # geom_sf(data = short_grass, fill = 'grey50', alpha = 0.6)+
  geom_sf(data = carm.point.sf, fill = 'red',color = 'red', size = 0.5, shape = 22)+
  annotate(
    'text', x = -Inf, y = Inf, hjust = -2.5, vjust = 2.3, label = 'A', size = 4
  )+
  theme_void()

mapreference
### FIGURE 1 plot### 

#plots 

biomass.coords <- read.csv('data/biomass_coords.csv')%>%
  st_as_sf(coords = c('x', 'y'), crs = st_crs(past.sf))

can.coords <- read.csv('data/manual_spatialcoords.csv')%>%
  st_as_sf(coords = c('UTM.E', 'UTM.N'), crs = st_crs(past.sf))

inst.coords <- read.csv('data/autogage_spatialcoords.csv')%>%
  st_as_sf(coords = c('UTM.E', 'UTM.N'), crs = st_crs(past.sf))


sitemap <- ggplot()+
  geom_sf(data = past.sf, fill = 'grey50')+
  geom_sf(data = carm.sf, aes(fill = forcats::fct_relevel(ecosite, c('Loamy Plains', 'Salt Flat', 'Sandy Plains', 'Other'))))+
  geom_sf(data = biomass.coords, color = 'black', size = 2, shape = 21)+
  
  geom_sf(data = can.coords, color = 'black', size = 2)+
  geom_sf(data = inst.coords, color = 'black', size = 2)+
  #geom_sf(data = biomass, shape = 3, size = 1.5)+
  labs(fill = 'Ecosite')+
  scale_fill_manual(values = precip.palette)+
  theme_void(base_size = 12)+
  annotate(
    'text', x = -Inf, y = Inf, hjust = -15, vjust = 2.3, label = 'B', size = 4
  )+
  ggspatial::annotation_scale(location = 'bl')+
  ggspatial::annotation_north_arrow(pad_y = unit(1, "cm"))+
  theme_void()+
  theme(legend.position = 'left')

fig.1 <- sitemap + inset_element(mapreference, 0, 0.55, 0.44, 1, align_to = 'full' ) 

#ggsave('figures/sitemap.pdf', fig.1, width = 13, height =13.5, units = 'cm')
