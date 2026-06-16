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

past.sf <- st_read('data/pastures.shp')
#-----------------------------------------------#
#Supplemental site fig 1
#-----------------------------------------------#

biomass <- read.csv('data/biomass_coords.csv')%>%
  st_as_sf(coords = c('x', 'y'), crs = st_crs(past.sf))

can.coords <- read.csv('data/manual_spatialcoords.csv')%>%
  st_as_sf(coords = c('UTM.E', 'UTM.N'), crs = st_crs(past.sf))

inst.coords <- read.csv('data/autogage_spatialcoords.csv')%>%
  st_as_sf(coords = c('UTM.E', 'UTM.N'), crs = st_crs(past.sf))

blockpairs <- read.csv('plotpairs.csv')%>%
  dplyr::select(pasture, block)%>%
  mutate(pasture = ifelse(pasture == 'NH', '10S',pasture))%>%
  unique()


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

#-----------------------------------------------#
## Breakdown of biomass percentage by ecosite 
#Supplemental figure 3
#-----------------------------------------------#
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


#-----------------------------------------------#
######sensitivity to water year precip########
#supplemental fig 4
#-----------------------------------------------#


###
m.global.t.mass.wyr <- lm(data = veg.precip.df, total.mass.dm ~ wyr.ppt*ecosite )
m.global.wspg.wyr <- lm(data = veg.precip.df, wspg.dm ~ wyr.ppt*ecosite)
m.global.c3pg.wyr <- lm(data = veg.precip.df, c3pg.dm ~ wyr.ppt*ecosite)
m.global.forb.wyr <- lm(data = veg.precip.df, forb.dm ~ wyr.ppt*ecosite)



#get grouping letters 
trends.tm.wyr  <- emtrends(m.global.t.mass.wyr , var = 'wyr.ppt', ~ ecosite)
letters.t.mass.wyr  <- multcomp::cld(trends.tm.wyr , Letters = letters)%>%
  mutate(response = 't.mass')
trends.wspg.wyr  <- emtrends(m.global.wspg.wyr , var = 'wyr.ppt', ~ ecosite)
letters.wspg.wyr <- multcomp::cld(trends.wspg.wyr , Letters = letters)%>%
  mutate(response = 'wspg')
trends.c3pg.wyr  <- emtrends(m.global.c3pg.wyr , var = 'wyr.ppt', ~ ecosite)
letters.c3pg.wyr  <- multcomp::cld(trends.c3pg.wyr , Letters = letters)%>%
  mutate(response = 'c3pg')
trends.forb.wyr  <- emtrends(m.global.forb.wyr , var = 'wyr.ppt', ~ ecosite)
letters.forb.wyr  <- multcomp::cld(trends.forb.wyr , Letters = letters)%>%
  mutate(response = 'forb')

letters.df.wyr  <- letters.t.mass.wyr %>%
  bind_rows(letters.wspg.wyr )%>%
  bind_rows(letters.c3pg.wyr )%>%
  bind_rows(letters.forb.wyr )

#get contrasts
contrast.t.mass.wyr  <- as.data.frame(contrast(emmeans::emtrends(m.global.t.mass.wyr , var = c('wyr.ppt'), ~ecosite), method = "pairwise")) %>%
  mutate(ecosite1 = str_trim(str_extract(contrast, "^[^-]+")),
         ecosite2 = str_trim(str_extract(contrast, "[^-]+$")),
         response = 't.mass') %>%
  filter(p.value < 0.05) # Filter only significant contrasts

contrast.wspg.wyr  <- as.data.frame(contrast(emmeans::emtrends(m.global.wspg.wyr , var = c('wyr.ppt'), ~ecosite), method = "pairwise"))%>%
  mutate(ecosite1 = str_trim(str_extract(contrast, "^[^-]+")),
         ecosite2 = str_trim(str_extract(contrast, "[^-]+$")),
         response = 'wspg') %>%
  filter(p.value < 0.05) # Filter only significant contrasts

contrast.c3pg.wyr  <- as.data.frame(contrast(emmeans::emtrends(m.global.c3pg.wyr , var = c('wyr.ppt'), ~ecosite), method = "pairwise"))%>%
  mutate(ecosite1 = str_trim(str_extract(contrast, "^[^-]+")),
         ecosite2 = str_trim(str_extract(contrast, "[^-]+$")),
         response = 'c3pg') %>%
  filter(p.value < 0.05) # Filter only significant contrasts

contrast.forb.wyr  <- as.data.frame(contrast(emmeans::emtrends(m.global.forb.wyr , var = c('wyr.ppt'), ~ecosite), method = "pairwise"))%>%
  mutate(ecosite1 = str_trim(str_extract(contrast, "^[^-]+")),
         ecosite2 = str_trim(str_extract(contrast, "[^-]+$")),
         response = 'forb') %>%
  filter(p.value < 0.05) # Filter only significant contrasts

contrast_df.wyr <- bind_rows(contrast.t.mass.wyr , contrast.wspg.wyr )%>%
  bind_rows(contrast.c3pg.wyr )%>%
  bind_rows(contrast.forb.wyr )

slps.global.wyr  <- as.data.frame(emmeans::emtrends(m.global.t.mass.wyr , var = c('wyr.ppt'), ~ecosite))%>%
  mutate(response = 't.mass')%>%
  bind_rows(as.data.frame(emmeans::emtrends(m.global.wspg.wyr , var = c('wyr.ppt'), ~ecosite))%>%
              mutate(response = 'wspg'))%>%
  bind_rows(as.data.frame(emmeans::emtrends(m.global.c3pg.wyr , var = c('wyr.ppt'), ~ecosite))%>%
              mutate(response = 'c3pg'))%>%
  bind_rows(as.data.frame(emmeans::emtrends(m.global.forb.wyr , var = c('wyr.ppt'), ~ecosite))%>%
              mutate(response = 'forb'))
slps.global.wyr %>%
  #filter(ecosite == 'Loamy')%>%
  ggplot( aes(x = fct_relevel(response, c('t.mass', 'wspg', 'c3pg', 'forb')), fill = ecosite))+
  geom_segment(aes(y = lower.CL, yend = upper.CL), position = position_dodge(0.5), linewidth = 1)+
  geom_point(aes(y = wyr.ppt.trend), position = position_dodge(0.5), shape = 21, size = 4)+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw()+
  labs(x = 'Functional group', y = 'ANPP Sensitivity')+
  scale_x_discrete(labels = c('Total', "Warm season", "Cool season", 'Forb'))+
  theme(text = element_text(size = 15))


require('ggpubr')
library('ggplot2')
ggplot(data = slps.global.wyr )+
  geom_segment(aes(x = ecosite, y = lower.CL, yend = upper.CL), position = position_dodge(0.5), linewidth = 1)+
  geom_point(aes(x = ecosite, y = wyr.ppt.trend, fill = ecosite), position = position_dodge(0.5), shape = 21, size = 4)+
  geom_text(data = letters.df.wyr , aes(x = ecosite, y = upper.CL + 1,  label = .group))+
  # geom_bracket(aes(xmin = ecosite1, xmax = ecosite2, y.position = 6.2), position = position_dodge(width = 0.5))+
  #geom_linerange(data = contrast_df, aes(ymin = ecosite1, ymax = ecosite2, x = 6.2, xmax = 6.2), position = position_dodge(width = 0.5), linewidth = 1)+
  #geom_labelsegment(data = contrast_df, aes(y = ecosite1, yend = ecosite2, x = 8, xend = 8, label= round(p.value*1000)/1000), position = position_dodge(width = 0.5))+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  #coord_flip()+
  theme_bw()+
  labs(x = 'Functional group', y = bquote(ANPP~Sensitivity~kg~ha^-1~mm^-1), fill = "Ecosite")+
  #scale_x_discrete(labels = c('Total', "Warm season", "Cool season", 'Forb'))+
  facet_grid(cols = vars(fct_relevel(response, c('t.mass', 'wspg', 'c3pg', 'forb'))), labeller = as_labeller(fun.grp.labels), switch = 'x')+
  theme(text = element_text(size = 12),
        panel.spacing = unit(0,'lines'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.spacing = unit(0.1, 'cm'))


