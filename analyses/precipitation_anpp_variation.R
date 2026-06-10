#create figure 2 precipitation and anpp variation 
`%notin%` <- Negate(`%in%`) 
#library
library(tidyverse)
library(sf)
library(sp)
library(stars)
library(raster)
library(patchwork)


#library(ggmap)
library(maps)
library(mapdata)
library(sfheaders)
library(lme4)
require(sjPlot)
library(egg)
require(ggpmisc)
library(ggspatial)

precip.palette <- c(  '#F78154','#389DE5', '#5DDCAC' ,'#9368B7')


#precip at pasture level 
monthly.ppt <- read.csv( 'data/monthly_ppt_past.csv')%>%
  as.data.frame()

#precip at the plotpair level
precip.vars <- read.csv('data/CPER_precip_variables.csv')

#precip
past.df <- read.csv('data/cper_broom_shp_v2.csv')

past.sf <- past.df %>%
  st_as_sf(coords = c("long", "lat")) %>%
  group_by(id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

######blockpairs
ppt.vars.sf <- monthly.ppt %>%
  as.data.frame()%>%
  mutate(water.yr = ifelse(month < 8, year, year + 1))%>%
  filter(water.yr < 2024 & water.yr > 2010)%>%
  group_by(id, water.yr)%>%
  summarise(gs.ppt = sum(month.ppt[month %in% 4:7]),
            cs.ppt = sum(month.ppt[month %in% c(11,12,1,2,3)]),
            f.ppt =  sum(month.ppt[month %in% 8:10]))%>%
  merge(past.sf, by = 'id')%>%
  st_as_sf()

###############

#figure 2 spatial vairablity in precip
mapped.precip.plt <- ggplot()+
  geom_sf(data = ppt.vars.sf, aes(fill = gs.ppt))+
  scale_fill_viridis_c(option = 'turbo', direction = -1)+
  facet_wrap(~water.yr, ncol = 5)+
  labs(title = 'Growing Season Precipitation', fill = 'GS PPT (mm)')+
  #theme()
  theme_void(base_size = 13)+
  theme(text=element_text(family="sans"), plot.title = element_text(hjust = 0.5),
        legend.position = c(.8, .15),
        legend.direction = 'horizontal',
        legend.title.position = 'top')
ggspatial::annotation_scale(location = 'br')


ggsave('figures/figure_2.pdf', mapped.precip.plt, width = 18, height = 15, units = 'cm')
