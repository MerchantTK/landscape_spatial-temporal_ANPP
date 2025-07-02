####generate manuscript figs
#TKM
#5.6.24

`%notin%` <- Negate(`%in%`) 
#library
library(tidyverse)
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
#set up shape files 

past.df <- read.csv('data/cper_broom_shp_v2.csv')

past.sf <- past.df %>%
  st_as_sf(coords = c("long", "lat")) %>%
  group_by(id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

#precip at pasture level 
monthly.ppt <- read.csv( 'data/monthly_ppt_past.csv')%>%
  as.data.frame()
#precip at the plotpair level
precip.vars <- read.csv('data/CPER_precip_variables.csv')

#plots 
biomass.coords <- read.csv('data/AGM_Biomass_Widecln_attr_2024-03-26.csv')%>%
  rename_all(tolower)%>%
  filter(is.na(x) == F)%>%
  mutate(year = yearsampled)%>%
  ungroup()%>%
  unite(plotpasture, plot, pasture, remove = F)%>%
  dplyr::select(pasture, plot, plotpasture, ecosite, x, y )%>%
  unique()
plot.coords.sf <- st_as_sf(biomass.coords, coords =  c('x','y'), crs = st_crs(past.sf))

#full dataframe
veg.precip.df <- read.csv('data/legacy_full_df.csv')

##Site figure




ecosites.sf <- st_read('data/AGM_ecosites_by_pastures2/AGM_ecosites_by_pastures2.shp')%>%
  mutate(ecosite = ifelse(ECOSITE1 %notin% c('Loamy Plains', 'Sandy Plains', 'Salt Flat'), 'Other', ECOSITE1))

st_crs(ecosites.sf) <- st_crs(past.sf)

carm.sf <- st_read('data/CARM_data/CARM_data.shp')%>%
  mutate(ecosite = ifelse(ECOSITE1 %notin% c('Loamy Plains', 'Sandy Plains', 'Salt Flat'), 'Other', ECOSITE1))
st_crs(carm.sf) <- st_crs(past.sf)




#ggsave('sitemap.png', sitemap, width = 6, height = 5)


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




spatial.cv <- ppt.vars.sf%>%
  st_drop_geometry()%>%
  group_by(id)%>%
  mutate(mean.gs.ppt = mean(gs.ppt, na.rm = T),
         gs.ppt.dm = gs.ppt - mean.gs.ppt)%>%
  ungroup()%>%
  mutate(no.spatial.gs.ppt= gs.ppt.dm + mean(gs.ppt, na.rm = T))%>%
  #filter(id %in% veg.precip.df$id)%>%
  group_by(water.yr)%>%
  summarise(cv= sd(no.spatial.gs.ppt, na.rm = T)/abs(mean(no.spatial.gs.ppt, na.rm = T)),
            sd = sd(no.spatial.gs.ppt, na.rm = T),
            variation = 'Spatio-temporal')



temp.cv <- ppt.vars.sf%>%
 # st_drop_geometry()%>%
  group_by(id)%>%
  summarise(cv = sd(gs.ppt, na.rm = T)/mean(gs.ppt, na.rm = T),
            mgsp = mean(gs.ppt, na.rm = T),
            variation = 'Temporal')

cv.cross <- temp.cv%>%
  st_drop_geometry()%>%
  dplyr::select(cv, variation)%>%
  bind_rows(spatial.cv)



ggplot(temp.cv, aes(fill = cv))+
  geom_sf()+
  scale_fill_viridis_c(option = 'turbo')

ggplot(temp.cv, aes(fill = mgsp))+
  geom_sf()+
  scale_fill_viridis_c(option = 'turbo', direction = -1)



monthly.ppt %>%
  as.data.frame()%>%
  mutate(water.yr = ifelse(month < 8, year, year + 1))%>%
  filter(year == 2023)%>%
  filter(month > 4 & month < 9)%>%
  merge(past.sf, by = 'id')%>%
  st_as_sf()%>%
  ggplot()+
  scale_fill_viridis_c(option = 'turbo', direction = -1)+
  geom_sf(aes(fill = month.ppt))+
  facet_wrap(~month)


mean(ppt.vars.sf$f.ppt, na.rm = T)
mean(ppt.vars.sf$cs.ppt, na.rm = T)
mean(ppt.vars.sf$gs.ppt, na.rm = T)


#cv of mass 
spatial.mass.cv <- veg.precip.df%>%
  mutate(no.spatial.total.mass = total.mass.dm + mean(total.mass, na.rm = T))%>%
  group_by(water.yr)%>%
  summarise(cv= sd(no.spatial.total.mass, na.rm = T)/abs(mean(no.spatial.total.mass, na.rm = T)),
            sd = sd(no.spatial.total.mass, na.rm = T),
            variation = 'Spatio-temporal')



temp.mass.cv <- veg.precip.df%>%
  mutate(no.spatial.total.mass = total.mass.dm + mean(total.mass, na.rm = T))%>%
  group_by(pairblock)%>%
  summarise(cv = sd(no.spatial.total.mass, na.rm = T)/abs(mean(no.spatial.total.mass, na.rm = T)),
            mgsp = abs(mean(no.spatial.total.mass, na.rm = T)),
            variation = 'Temporal')


cv.mass.cross <- temp.mass.cv%>%
  st_drop_geometry()%>%
  dplyr::select(cv, variation)%>%
  bind_rows(spatial.mass.cv)

mass.cv.means <- aggregate(cv.mass.cross, cv ~ variation, FUN = 'mean' )

precip.cv.means <- aggregate(cv.cross, cv ~ variation, FUN = 'mean' )

mass.cv.plt <- ggplot(cv.mass.cross, aes(x = cv))+
  geom_density( aes(x = cv, fill = variation), alpha = 0.75)+
  xlim(0,1)+
  geom_text(label = 'Mass', aes(x = 0.75, y = 4.75), size = 4)+
  labs(y = 'Density', x = 'CV', fill = 'Variation', color = 'Variation')+
  geom_vline(data = mass.cv.means, aes(xintercept = cv, color = variation), linewidth = 1)+
  ylim(0,25)+
  geom_text( aes(x = -Inf, y = Inf, label = 'G'),
             hjust = -1, vjust = 1)+
  scale_fill_manual(values = c('#008080', 'tomato2'))+
  scale_color_manual(values = c('#008080', 'tomato2'))+
  theme_bw(base_size = 11)+
  theme(text=element_text(family="sans"))

precip.cv.plt <- ggplot(cv.cross, aes(x = cv))+
  geom_density( aes(x = cv, fill = variation), alpha = 0.75)+
  geom_text(label = 'Precipitation', aes(x = 0.75, y = 4.75), size = 4)+
  geom_vline(data = precip.cv.means, aes(xintercept = cv, color = variation), linewidth = 1)+
  scale_fill_manual(values = c('#008080', 'tomato2'))+
  scale_color_manual(values = c('#008080', 'tomato2'))+
  xlim(0,1)+
  geom_text( aes(x = -Inf, y = Inf, label = 'F'),
            hjust = -1, vjust = 1)+
  labs(y = 'Density', x = 'CV', fill = 'Variation', color = 'Variation')+
  theme_bw(base_size = 11)+
  theme(text=element_text(family="sans"))

cvplt <- precip.cv.plt + mass.cv.plt + plot_layout(nrow = 2, axes = 'collect', guides = 'collect')

#ggsave('cvplot.png', cvplt, width = 4.5, height = 3.5)



ggplot(temp.mass.cv, aes(x = cv.temp))+
  geom_density(fill = 'tomato2')+
  geom_density( aes(x = cv.sp),fill = '#008080')+
  xlim(0,1)


spatial.cv%>%
  merge(spatial.mass.cv, by = c('water.yr', 'variation'))%>%
  rename(cv.precip = cv.x,
         cv.mass = cv.y,
         sd.precip = sd.x,
         sd.mass = sd.y)%>%
  ggplot(aes(x = cv.precip, y = cv.mass))+
  geom_point()+
  ggpmisc::stat_poly_eq(ggpmisc::use_label('p.value'), method = 'lm')+
  geom_smooth(method = 'lm')
  
#figure 2 spatial vairablity in precip
mapped.precip.plt <- ggplot()+
  geom_sf(data = ppt.vars.sf, aes(fill = gs.ppt))+
  scale_fill_viridis_c(option = 'turbo', direction = -1)+
  facet_wrap(~water.yr, ncol = 5)+
  labs(title = 'Growing Season Precipitation', fill = 'GS PPT')+
  #theme()
  theme_void(base_size = 13)+
  theme(text=element_text(family="sans"), plot.title = element_text(hjust = 0.5),
        legend.position = c(.8, .15),
        legend.direction = 'horizontal',
        legend.title.position = 'top')
  ggspatial::annotation_scale(location = 'br')
  
mapped.precip.plt

ppt.vars.sf%>%
  group_by(id, geometry)%>%
  summarize(mean.gs = mean(gs.ppt, na.rm = T))%>%
  ggplot()+
  geom_sf( aes(fill = mean.gs))+
  scale_fill_viridis_c(option = 'turbo', direction = -1)+
  #facet_wrap(~water.yr, ncol = 5)+
  labs(title = 'Mean Growing Season Precipitation', fill = 'GS PPT')+
  #theme()
  theme_void(base_size = 13)+
  theme(#text=element_text(family="sans"), plot.title = element_text(hjust = 0.5),
        #legend.position = c(.8, .15),
        legend.direction = 'horizontal',
        legend.title.position = 'top')
ggspatial::annotation_scale(location = 'br')

#ggsave('mapfig.png', mapped.precip.plt, width = 6, height = 5)

# timeline of gs precip
timeline.precip.plt <- ppt.vars.sf%>%
  ggplot(aes(x = water.yr, y = gs.ppt, group = id))+
  geom_line()+
  labs(x = 'Water year', y = 'GS PPT mm')+
  geom_text(aes(x = -Inf, y = Inf, label = 'A'),
            hjust = -1, vjust = 1.5)+
  scale_x_continuous(breaks = seq(2012, 2022, 2))+
  theme_bw(base_size = 11)+
  theme(text=element_text(family="sans"), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank())



#timeline of veg vars 
timeline.veg.labels <- veg.precip.df%>%
  dplyr::select(total.mass, wspg, c3pg, forb, water.yr, pairblock, ecosite)%>%
  pivot_longer(-c(water.yr, pairblock, ecosite), names_to = 'fun.grp', values_to = 'mass')%>%
  distinct(fun.grp)%>%
  mutate(label = LETTERS[2:5],
         x = -Inf,       # X position inside each panel
         y = Inf)       # Y position inside each pan)

timeline.veg  <- veg.precip.df%>%
  dplyr::select(total.mass, wspg, c3pg, forb, water.yr, pairblock, ecosite)%>%
  pivot_longer(-c(water.yr, pairblock, ecosite), names_to = 'fun.grp', values_to = 'mass')%>%
  filter(fun.grp != 'total.mass')%>%
  mutate(fun.grp = forcats::fct_relevel(fun.grp, c( 'wspg', 'c3pg', 'forb')))%>%
  ggplot()+
  # geom_line(aes(x = as.numeric(as.character(water.yr)), y = mass, group = pairblock, color = ecosite), linewidth = 0.7)+
  stat_summary(geom = 'pointrange', fun.data = mean_cl_normal, aes(x = as.numeric(as.character(water.yr)), y = mass, group = ecosite, color = ecosite))+
  stat_summary(geom = 'line', fun.data = mean_cl_normal, aes(x = as.numeric(as.character(water.yr)), y = mass, group = ecosite, color = ecosite))+
  geom_text(data = timeline.veg.labels[2:4,], aes(x = x, y = y, label = label),
            hjust = -1, vjust = 1)+
  facet_grid(rows = vars(forcats::fct_relevel(fun.grp, c( 'wspg', 'c3pg', 'forb'))), labeller = as_labeller(fun.grp.labels))+
  labs(x = 'Water year', y =  bquote('ANPP kg ha'^-1), color = 'Ecosite')+
  scale_x_continuous(breaks = seq(2010, 2022, 2), limits = c(2011, 2023))+
  scale_color_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  theme(text=element_text(family="sans"), 
        #axis.title = element_text(vjust = 10),
        legend.position = 'bottom')

timeline.t.mass  <- veg.precip.df%>%
  dplyr::select(total.mass, wspg, c3pg, forb, water.yr, pairblock, ecosite)%>%
  pivot_longer(-c(water.yr, pairblock, ecosite), names_to = 'fun.grp', values_to = 'mass')%>%
  filter(fun.grp == 'total.mass')%>%
 # mutate(fun.grp = forcats::fct_relevel(fun.grp, c( 'wspg', 'c3pg', 'forb')))%>%
  ggplot()+
  # geom_line(aes(x = as.numeric(as.character(water.yr)), y = mass, group = pairblock, color = ecosite), linewidth = 0.7)+
  stat_summary(geom = 'pointrange', fun.data = mean_cl_normal, aes(x = as.numeric(as.character(water.yr)), y = mass, group = ecosite, color = ecosite))+
  stat_summary(geom = 'line', fun.data = mean_cl_normal, aes(x = as.numeric(as.character(water.yr)), y = mass, group = ecosite, color = ecosite))+
  geom_text(data = timeline.veg.labels[1,], aes(x = x, y = y, label = label),
            hjust = -1, vjust = 1)+
  facet_grid(rows = vars(fun.grp), labeller = as_labeller(fun.grp.labels), scale = 'free_y')+
  labs(x = 'Water year', y = bquote('ANPP kg ha'^-1), color = 'Ecosite')+
  scale_x_continuous(breaks = seq(2010, 2022, 2), limits = c(2011, 2023))+
  scale_color_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  theme(text=element_text(family="sans"), 
        #axis.title = element_text(vjust = 10),
        legend.position = 'none')

layout <- '
ACD
BCD'
timeline.precip.plt + timeline.t.mass + timeline.veg + cvplt + plot_layout(design = layout)



fun.grp.labels <- c( 
  `total.mass` = 'Total Mass',
  `wspg` = 'Warm',
  `c3pg` = 'Cool',
  `forb` = 'Forb')

timeline.t.mass.plt <- veg.precip.df%>%
  dplyr::select(total.mass, wspg, c3pg, forb, water.yr, pairblock, ecosite)%>%
  pivot_longer(-c(water.yr, pairblock, ecosite), names_to = 'fun.grp', values_to = 'mass')%>%
  filter(fun.grp == 'total.mass')%>%
 # mutate(fun.grp = forcats::fct_relevel(fun.grp, c( 'wspg', 'c3pg', 'forb')))%>%
  ggplot()+
  geom_line(aes(x = water.yr, y = mass, group = pairblock, color = ecosite), linewidth = 0.7)+
  facet_grid(rows = vars(fun.grp), scales = 'free_y', labeller = as_labeller(fun.grp.labels))+
  labs(x = 'Water year', y = 'Biomass production (kg/ha)', color = 'Ecosite')+
  scale_x_continuous(breaks = seq(2010, 2022, 2), limits = c(2011, 2023))+
  scale_color_manual(values = precip.palette)+
  guides(color = 'none')+
  theme_bw(base_size = 13)+
  theme(text=element_text(family="sans"))




timeline.veg.plt <- veg.precip.df%>%
  dplyr::select(total.mass, wspg, c3pg, forb, water.yr, pairblock, ecosite)%>%
  pivot_longer(-c(water.yr, pairblock, ecosite), names_to = 'fun.grp', values_to = 'mass')%>%
  filter(fun.grp != 'total.mass')%>%
  mutate(fun.grp = forcats::fct_relevel(fun.grp, c( 'wspg', 'c3pg', 'forb')))%>%
  ggplot()+
  geom_line(aes(x = water.yr, y = mass, group = pairblock, color = ecosite), linewidth = 0.7)+
  facet_grid(rows = vars(fun.grp), labeller = as_labeller(fun.grp.labels))+
  labs(x = 'Water year', y = 'Biomass production (kg/ha)', color = 'Ecosite')+
  scale_x_continuous(breaks = seq(2010, 2022, 2), limits = c(2011, 2023))+
  scale_color_manual(values = precip.palette)+
  theme_bw(base_size = 13)+
  theme(text=element_text(family="sans"), 
        #axis.title = element_text(vjust = 10),
        legend.position = 'bottom')



layout <- "
AAAAAAAAA
AAAAAAAAA
BBBBBBBBB
BBBBBBBBB
BBBBBBBBB
CCCCCCCCC
CCCCCCCCC
CCCCCCCCC
CCCCCCCCC
"

fig2_timeline <-  timeline.precip.plt + timeline.t.mass.plt + timeline.veg.plt+ plot_layout(design = layout, axes= 'collect', axis_titles = 'collect') +
plot_annotation(tag_levels = "A")


#ggsave('fig_timeline.png', fig2_timeline, height = 7.5, width = 7.5)






#model set table: total mass######
m.tm.gs <- lmer(total.mass ~ gs.ppt*ecosite +  (1|pairblock) + (1|water.yr),
                data = veg.precip.df)

m.tm.gs.cs <- lmer(total.mass ~ gs.ppt*ecosite + f.ppt*ecosite  +  (1|pairblock) + (1|water.yr),
                   data = veg.precip.df)

m.tm.gs.f <- lmer(total.mass ~ gs.ppt*ecosite + f.ppt*ecosite  +  (1|pairblock) + (1|water.yr),
                  data = veg.precip.df)

m.tm.gs.gs1 <- lmer(total.mass ~ gs.ppt*ecosite + lag.gs*ecosite  +  (1|pairblock) + (1|water.yr),
                    data = veg.precip.df)

m.tm.gs.f.cs <- lmer(total.mass ~ gs.ppt*ecosite + f.ppt*ecosite + lag.cs*ecosite +  (1|pairblock) + (1|water.yr),
                     data = veg.precip.df)

m.tm.gs.f.gs1 <- lmer(total.mass ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite +  (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)


m.tm.gs.cs.gs1 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite +  (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

m.tm.gs.cs.f.gs1 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + f.ppt*ecosite + (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)
AIC( m.tm.gs, m.tm.gs.cs, m.tm.gs.f, m.tm.gs.gs1, m.tm.gs.f.cs, m.tm.gs.f.gs1, m.tm.gs.cs.gs1, m.tm.gs.cs.f.gs1)
#add 2 season lags

#
m.tm.gs.cs.f.gs1.cs1 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite  + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                             data = veg.precip.df)
m.tm.gs.cs.f.gs1.f1 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite  + lag.f*ecosite + (1|pairblock) + (1|water.yr),
                            data = veg.precip.df)
m.tm.gs.cs.f.gs1.gs2 <- lmer(total.mass ~ ecosite*gs.ppt + ecosite*cs.ppt + ecosite*f.ppt + ecosite*lag.gs  + ecosite*lag2.gs + (1|pairblock) + (1|water.yr),
                             data = veg.precip.df)

#
m.tm.gs.cs.f.cs1 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)

m.tm.gs.cs.f.f1 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                        data = veg.precip.df)

m.tm.gs.cs.f.gs2 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)

#
m.tm.gs.cs.gs1.cs1 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                           data = veg.precip.df)

m.tm.gs.cs.gs1.f1 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)

m.tm.gs.cs.gs1.gs2 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                           data = veg.precip.df)

#


m.tm.gs.f.gs1.cs1 <- lmer(total.mass ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)

m.tm.gs.f.gs1.f1 <- lmer(total.mass ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)

m.tm.gs.f.gs1.gs2 <- lmer(total.mass ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)


#
m.tm.gs.cs.cs1 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

m.tm.gs.cs.f1 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite +  lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)

m.tm.gs.cs.gs2 <- lmer(total.mass ~ gs.ppt*ecosite + cs.ppt*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

#

m.tm.gs.f.cs1 <- lmer(total.mass ~ gs.ppt*ecosite + f.ppt*ecosite  + lag.cs*ecosite+ (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)

m.tm.gs.f.f1 <- lmer(total.mass ~ gs.ppt*ecosite + f.ppt*ecosite  + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                     data = veg.precip.df)

m.tm.gs.f.gs2 <- lmer(total.mass ~ gs.ppt*ecosite + f.ppt*ecosite  + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)



#

m.tm.gs.gs1.cs1 <- lmer(total.mass ~ gs.ppt*ecosite  + lag.gs*ecosite + lag.cs*ecosite+ (1|pairblock) + (1|water.yr),
                        data = veg.precip.df)

m.tm.gs.gs1.f1 <- lmer(total.mass ~ gs.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

m.tm.gs.gs1.gs2 <- lmer(total.mass ~ gs.ppt*ecosite +  lag.gs*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                        data = veg.precip.df)



aic.tm.df <-AIC( m.tm.gs, m.tm.gs.cs, m.tm.gs.f, m.tm.gs.gs1, m.tm.gs.f.cs, m.tm.gs.f.gs1, m.tm.gs.cs.gs1, m.tm.gs.cs.f.gs1,
                 m.tm.gs.cs.f.gs1.cs1, m.tm.gs.cs.f.gs1.f1, m.tm.gs.cs.f.gs1.gs2,
                 m.tm.gs.cs.f.cs1, m.tm.gs.cs.f.f1, m.tm.gs.cs.f.gs2,
                 m.tm.gs.cs.gs1.cs1, m.tm.gs.cs.gs1.f1, m.tm.gs.cs.gs1.gs2,
                 m.tm.gs.f.gs1.cs1, m.tm.gs.f.gs1.f1, m.tm.gs.f.gs1.gs2,
                 m.tm.gs.cs.cs1, m.tm.gs.cs.f1, m.tm.gs.cs.gs2,
                 m.tm.gs.f.cs1, m.tm.gs.f.f1, m.tm.gs.f.gs2,
                 m.tm.gs.gs1.cs1, m.tm.gs.gs1.f1, m.tm.gs.gs1.gs2)%>%
  mutate(delta.aic.tm =  AIC - min(AIC) )

#model checks on best model 
plot(m.tm.gs.cs.f.gs1.gs2)
qqmath(resid(m.tm.gs.cs.f.gs1.gs2))
plot(m.tm.gs.cs.f.gs1.gs2, rstudent(.) ~ hatvalues(.))


#write.csv(aic.tm.df, 'total_mass_AIC.csv')

#model set table: warm season ######
m.ws.gs <- lmer(wspg ~ gs.ppt*ecosite +  (1|pairblock) + (1|water.yr),
                data = veg.precip.df)

m.ws.gs.cs <- lmer(wspg ~ gs.ppt*ecosite + f.ppt*ecosite  +  (1|pairblock) + (1|water.yr),
                   data = veg.precip.df)

m.ws.gs.f <- lmer(wspg ~ gs.ppt*ecosite + f.ppt*ecosite  +  (1|pairblock) + (1|water.yr),
                  data = veg.precip.df)

m.ws.gs.gs1 <- lmer(wspg ~ gs.ppt*ecosite + lag.gs*ecosite  +  (1|pairblock) + (1|water.yr),
                    data = veg.precip.df)

m.ws.gs.f.cs <- lmer(wspg ~ gs.ppt*ecosite + f.ppt*ecosite + lag.cs*ecosite +  (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)

m.ws.gs.f.gs1 <- lmer(wspg ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite +  (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)


m.ws.gs.cs.gs1 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite +  (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

m.ws.gs.cs.f.gs1 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + f.ppt*ecosite + (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)
AIC( m.ws.gs, m.ws.gs.cs, m.ws.gs.f, m.ws.gs.gs1, m.ws.gs.f.cs, m.ws.gs.f.gs1, m.ws.gs.cs.gs1, m.ws.gs.cs.f.gs1)
#add 2 season lags

#
m.ws.gs.cs.f.gs1.cs1 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite  + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)
m.ws.gs.cs.f.gs1.f1 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite  + lag.f*ecosite + (1|pairblock) + (1|water.yr),
                             data = veg.precip.df)
m.ws.gs.cs.f.gs1.gs2 <- lmer(wspg ~ ecosite*gs.ppt + cs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite  + lag2.gs*ecosite + (1|pairblock) + (1|water.yr),
                             data = veg.precip.df)

#
m.ws.gs.cs.f.cs1 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)

m.ws.gs.cs.f.f1 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)

m.ws.gs.cs.f.gs2 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)

#
m.ws.gs.cs.gs1.cs1 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                           data = veg.precip.df)

m.ws.gs.cs.gs1.f1 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)

m.ws.gs.cs.gs1.gs2 <- lmer(wspg ~ ecosite*gs.ppt + ecosite*cs.ppt + ecosite*lag.gs + ecosite*lag2.gs+ (1|pairblock) + (1|water.yr),
                           data = veg.precip.df)

#


m.ws.gs.f.gs1.cs1 <- lmer(wspg ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)

m.ws.gs.f.gs1.f1 <- lmer(wspg ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)

m.ws.gs.f.gs1.gs2 <- lmer(wspg ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)


#
m.ws.gs.cs.cs1 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                           data = veg.precip.df)

m.ws.gs.cs.f1 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite +  lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)

m.ws.gs.cs.gs2 <- lmer(wspg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                           data = veg.precip.df)

#

m.ws.gs.f.cs1 <- lmer(wspg ~ gs.ppt*ecosite + f.ppt*ecosite  + lag.cs*ecosite+ (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)

m.ws.gs.f.f1 <- lmer(wspg ~ gs.ppt*ecosite + f.ppt*ecosite  + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                     data = veg.precip.df)

m.ws.gs.f.gs2 <- lmer(wspg ~ gs.ppt*ecosite + f.ppt*ecosite  + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)



#

m.ws.gs.gs1.cs1 <- lmer(wspg ~ gs.ppt*ecosite  + lag.gs*ecosite + lag.cs*ecosite+ (1|pairblock) + (1|water.yr),
                        data = veg.precip.df)

m.ws.gs.gs1.f1 <- lmer(wspg ~ gs.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

m.ws.gs.gs1.gs2 <- lmer(wspg ~ gs.ppt*ecosite +  lag.gs*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                        data = veg.precip.df)



aic.ws.df <-AIC( m.ws.gs, m.ws.gs.cs, m.ws.gs.f, m.ws.gs.gs1, m.ws.gs.f.cs, m.ws.gs.f.gs1, m.ws.gs.cs.gs1, m.ws.gs.cs.f.gs1,
     m.ws.gs.cs.f.gs1.cs1, m.ws.gs.cs.f.gs1.f1, m.ws.gs.cs.f.gs1.gs2,
     m.ws.gs.cs.f.cs1, m.ws.gs.cs.f.f1, m.ws.gs.cs.f.gs2,
     m.ws.gs.cs.gs1.cs1, m.ws.gs.cs.gs1.f1, m.ws.gs.cs.gs1.gs2,
     m.ws.gs.f.gs1.cs1, m.ws.gs.f.gs1.f1, m.ws.gs.f.gs1.gs2,
     m.ws.gs.cs.cs1, m.ws.gs.cs.f1, m.ws.gs.cs.gs2,
     m.ws.gs.f.cs1, m.ws.gs.f.f1, m.ws.gs.f.gs2,
     m.ws.gs.gs1.cs1, m.ws.gs.gs1.f1, m.ws.gs.gs1.gs2)%>%
  mutate(delta.aic.ws =  AIC - min(AIC) )

#write.csv(aic.ws.df, 'warm_season_AIC.csv')

#model checks on best model 
plot(m.ws.gs.cs.gs1.gs2 )
qqmath(resid(m.ws.gs.cs.gs1.gs2 ))
plot(m.ws.gs.cs.gs1.gs2 , rstudent(.) ~ hatvalues(.))



#model set table: cool season ######
m.c3.gs <- lmer(c3pg ~ gs.ppt*ecosite +  (1|pairblock) + (1|water.yr),
                data = veg.precip.df)

m.c3.gs.cs <- lmer(c3pg ~ gs.ppt*ecosite + f.ppt*ecosite  +  (1|pairblock) + (1|water.yr),
                   data = veg.precip.df)

m.c3.gs.f <- lmer(c3pg ~ gs.ppt*ecosite + f.ppt*ecosite  +  (1|pairblock) + (1|water.yr),
                  data = veg.precip.df)

m.c3.gs.gs1 <- lmer(c3pg ~ gs.ppt*ecosite + lag.gs*ecosite  +  (1|pairblock) + (1|water.yr),
                    data = veg.precip.df)

m.c3.gs.f.cs <- lmer(c3pg ~ gs.ppt*ecosite + f.ppt*ecosite + lag.cs*ecosite +  (1|pairblock) + (1|water.yr),
                     data = veg.precip.df)

m.c3.gs.f.gs1 <- lmer(c3pg ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite +  (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)


m.c3.gs.cs.gs1 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite +  (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

m.c3.gs.cs.f.gs1 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + f.ppt*ecosite + (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)
AIC( m.c3.gs, m.c3.gs.cs, m.c3.gs.f, m.c3.gs.gs1, m.c3.gs.f.cs, m.c3.gs.f.gs1, m.c3.gs.cs.gs1, m.c3.gs.cs.f.gs1)


#add 2 season lags

#
m.c3.gs.cs.f.gs1.cs1 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite  + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                             data = veg.precip.df)
m.c3.gs.cs.f.gs1.f1 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite  + lag.f*ecosite + (1|pairblock) + (1|water.yr),
                            data = veg.precip.df)
m.c3.gs.cs.f.gs1.gs2 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite  + lag2.gs*ecosite + (1|pairblock) + (1|water.yr),
                             data = veg.precip.df)

#
m.c3.gs.cs.f.cs1 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)

m.c3.gs.cs.f.f1 <- lmer(c3pg ~ ecosite*gs.ppt + ecosite*cs.ppt + ecosite*f.ppt + ecosite*lag.f+ (1|pairblock) + (1|water.yr),
                        data = veg.precip.df)

m.c3.gs.cs.f.gs2 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)

#
m.c3.gs.cs.gs1.cs1 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                           data = veg.precip.df)

m.c3.gs.cs.gs1.f1 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)

m.c3.gs.cs.gs1.gs2 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                           data = veg.precip.df)

#


m.c3.gs.f.gs1.cs1 <- lmer(c3pg ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)

m.c3.gs.f.gs1.f1 <- lmer(c3pg ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)

m.c3.gs.f.gs1.gs2 <- lmer(c3pg ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)


#
m.c3.gs.cs.cs1 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

m.c3.gs.cs.f1 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite +  lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)

m.c3.gs.cs.gs2 <- lmer(c3pg ~ gs.ppt*ecosite + cs.ppt*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

#

m.c3.gs.f.cs1 <- lmer(c3pg ~ gs.ppt*ecosite + f.ppt*ecosite  + lag.cs*ecosite+ (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)

m.c3.gs.f.f1 <- lmer(c3pg ~ gs.ppt*ecosite + f.ppt*ecosite  + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                     data = veg.precip.df)

m.c3.gs.f.gs2 <- lmer(c3pg ~ gs.ppt*ecosite + f.ppt*ecosite  + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)



#

m.c3.gs.gs1.cs1 <- lmer(c3pg ~ gs.ppt*ecosite  + lag.gs*ecosite + lag.cs*ecosite+ (1|pairblock) + (1|water.yr),
                        data = veg.precip.df)

m.c3.gs.gs1.f1 <- lmer(c3pg ~ gs.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

m.c3.gs.gs1.gs2 <- lmer(c3pg ~ gs.ppt*ecosite +  lag.gs*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                        data = veg.precip.df)



aic.c3.df <-AIC( m.c3.gs, m.c3.gs.cs, m.c3.gs.f, m.c3.gs.gs1, m.c3.gs.f.cs, m.c3.gs.f.gs1, m.c3.gs.cs.gs1, m.c3.gs.cs.f.gs1,
                 m.c3.gs.cs.f.gs1.cs1, m.c3.gs.cs.f.gs1.f1, m.c3.gs.cs.f.gs1.gs2,
                 m.c3.gs.cs.f.cs1, m.c3.gs.cs.f.f1, m.c3.gs.cs.f.gs2,
                 m.c3.gs.cs.gs1.cs1, m.c3.gs.cs.gs1.f1, m.c3.gs.cs.gs1.gs2,
                 m.c3.gs.f.gs1.cs1, m.c3.gs.f.gs1.f1, m.c3.gs.f.gs1.gs2,
                 m.c3.gs.cs.cs1, m.c3.gs.cs.f1, m.c3.gs.cs.gs2,
                 m.c3.gs.f.cs1, m.c3.gs.f.f1, m.c3.gs.f.gs2,
                 m.c3.gs.gs1.cs1, m.c3.gs.gs1.f1, m.c3.gs.gs1.gs2)%>%
  mutate(delta.aic.c3 =  AIC - min(AIC) )


#write.csv(aic.c3.df, 'cool_season_AIC.csv')


#model checks on best model 
plot(m.c3.gs.cs.f.f1   )
qqmath(resid(m.c3.gs.cs.f.f1  ))
plot(m.c3.gs.cs.f.f1   , rstudent(.) ~ hatvalues(.))


#model set table: forbs ######
m.fo.gs <- lmer(forb ~ gs.ppt*ecosite +  (1|pairblock) + (1|water.yr),
                data = veg.precip.df)

m.fo.gs.cs <- lmer(forb ~ gs.ppt*ecosite + f.ppt*ecosite  +  (1|pairblock) + (1|water.yr),
                   data = veg.precip.df)

m.fo.gs.f <- lmer(forb ~ gs.ppt*ecosite + f.ppt*ecosite  +  (1|pairblock) + (1|water.yr),
                  data = veg.precip.df)

m.fo.gs.gs1 <- lmer(forb ~ gs.ppt*ecosite + lag.gs*ecosite  +  (1|pairblock) + (1|water.yr),
                    data = veg.precip.df)

m.fo.gs.f.cs <- lmer(forb ~ gs.ppt*ecosite + f.ppt*ecosite + lag.cs*ecosite +  (1|pairblock) + (1|water.yr),
                     data = veg.precip.df)

m.fo.gs.f.gs1 <- lmer(forb ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite +  (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)


m.fo.gs.cs.gs1 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite +  (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

m.fo.gs.cs.f.gs1 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + f.ppt*ecosite + (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)
AIC( m.fo.gs, m.fo.gs.cs, m.fo.gs.f, m.fo.gs.gs1, m.fo.gs.f.cs, m.fo.gs.f.gs1, m.fo.gs.cs.gs1, m.fo.gs.cs.f.gs1)


#add 2 season lags

#
m.fo.gs.cs.f.gs1.cs1 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite  + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                             data = veg.precip.df)
m.fo.gs.cs.f.gs1.f1 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite  + lag.f*ecosite + (1|pairblock) + (1|water.yr),
                            data = veg.precip.df)

summary(m.fo.gs.cs.f.gs1.f1 )
m.fo.gs.cs.f.gs1.gs2 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite  + lag2.gs*ecosite + (1|pairblock) + (1|water.yr),
                             data = veg.precip.df)

#
m.fo.gs.cs.f.cs1 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)

m.fo.gs.cs.f.f1 <- lmer(forb ~ ecosite*gs.ppt + ecosite*cs.ppt + ecosite*f.ppt + ecosite*lag.f+ (1|pairblock) + (1|water.yr),
                        data = veg.precip.df)

m.fo.gs.cs.f.gs2 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + f.ppt*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)

#
m.fo.gs.cs.gs1.cs1 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                           data = veg.precip.df)

m.fo.gs.cs.gs1.f1 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)

m.fo.gs.cs.gs1.gs2 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.gs*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                           data = veg.precip.df)

#


m.fo.gs.f.gs1.cs1 <- lmer(forb ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)

m.fo.gs.f.gs1.f1 <- lmer(forb ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                         data = veg.precip.df)

m.fo.gs.f.gs1.gs2 <- lmer(forb ~ gs.ppt*ecosite + f.ppt*ecosite + lag.gs*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                          data = veg.precip.df)


#
m.fo.gs.cs.cs1 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + lag.cs*ecosite + (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

m.fo.gs.cs.f1 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite +  lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)

m.fo.gs.cs.gs2 <- lmer(forb ~ gs.ppt*ecosite + cs.ppt*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

#

m.fo.gs.f.cs1 <- lmer(forb ~ gs.ppt*ecosite + f.ppt*ecosite  + lag.cs*ecosite+ (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)

m.fo.gs.f.f1 <- lmer(forb ~ gs.ppt*ecosite + f.ppt*ecosite  + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                     data = veg.precip.df)

m.fo.gs.f.gs2 <- lmer(forb ~ gs.ppt*ecosite + f.ppt*ecosite  + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                      data = veg.precip.df)



#

m.fo.gs.gs1.cs1 <- lmer(forb ~ gs.ppt*ecosite  + lag.gs*ecosite + lag.cs*ecosite+ (1|pairblock) + (1|water.yr),
                        data = veg.precip.df)

m.fo.gs.gs1.f1 <- lmer(forb ~ gs.ppt*ecosite + lag.gs*ecosite + lag.f*ecosite+ (1|pairblock) + (1|water.yr),
                       data = veg.precip.df)

m.fo.gs.gs1.gs2 <- lmer(forb ~ gs.ppt*ecosite +  lag.gs*ecosite + lag2.gs*ecosite+ (1|pairblock) + (1|water.yr),
                        data = veg.precip.df)



aic.fo.df <-AIC( m.fo.gs, m.fo.gs.cs, m.fo.gs.f, m.fo.gs.gs1, m.fo.gs.f.cs, m.fo.gs.f.gs1, m.fo.gs.cs.gs1, m.fo.gs.cs.f.gs1,
                 m.fo.gs.cs.f.gs1.cs1, m.fo.gs.cs.f.gs1.f1, m.fo.gs.cs.f.gs1.gs2,
                 m.fo.gs.cs.f.cs1, m.fo.gs.cs.f.f1, m.fo.gs.cs.f.gs2,
                 m.fo.gs.cs.gs1.cs1, m.fo.gs.cs.gs1.f1, m.fo.gs.cs.gs1.gs2,
                 m.fo.gs.f.gs1.cs1, m.fo.gs.f.gs1.f1, m.fo.gs.f.gs1.gs2,
                 m.fo.gs.cs.cs1, m.fo.gs.cs.f1, m.fo.gs.cs.gs2,
                 m.fo.gs.f.cs1, m.fo.gs.f.f1, m.fo.gs.f.gs2,
                 m.fo.gs.gs1.cs1, m.fo.gs.gs1.f1, m.fo.gs.gs1.gs2)%>%
  mutate(delta.aic.fo =  AIC - min(AIC) )

#write.csv(aic.fo.df, 'forb_AIC.csv')

#model checks on best model 
plot(m.c3.gs.cs.f.f1   )
qqmath(resid(m.c3.gs.cs.f.f1  ))
plot(m.c3.gs.cs.f.f1   , rstudent(.) ~ hatvalues(.))



#### table of AIC####
models.tbl <- cbind(aic.tm.df, aic.ws.df, aic.c3.df, aic.fo.df)

#write.csv(models.tbl, 'all_models.csv')


####table of best models####

#table of best models


tab_model( m.tm.gs.cs.f.gs1.gs2, m.ws.gs.cs.gs1.gs2 , m.c3.gs.cs.f.f1, m.fo.gs.cs.f.f1,
          pred.labels = c('intercept (Loamy)', 'Salt Flats',  'Sandy', 'GrowSeason[t0]', 'ColdSeason[t0]', 'Fall[t0]', 'GrowSeason[t-1]', 'GrowSeason[t-2]',  'Salt Flats:GrowSeason[t0]', 'Sandy:GrowSeason[t0]', 'Salt Flats:ColdSeason[t0]', 'Sandy:ColdSeason[t0]', 'Salt Flats:Fall[t0]', 'Sandy:Fall[t0]','Salt Flats:GrowSeason[t0]', 'Sandy:GrowSeason[t-1]', 'Salt Flats:GrowSeason[t-2]', 'Sandy:GrowSeason[t-2]', 'Fall[t-1]', 'Salt Flats:Fall[t-1]', 'Sandy:Fall[t-1]' ),
          dv.labels = c('Total', 'Warm', 'Cool', 'Forb'))
###### Best model output

var_names <- c(`gs.ppt` = 'GrowSeason[t0]',
               `cs.ppt` = 'ColdSeason[t0]',
               `f.ppt` = 'Fall[t0]',
               `lag.gs` = 'GrowSeason[t-1]',
               `lag.f` = 'Fall[t-1]',
               `lag2.gs` = 'GrowSeason[t-2]')

#warm season 
aic.ws.df 

sig.ws.gs <- emmeans::emtrends( m.ws.gs.cs.gs1.gs2 , ~ecosite, var = 'gs.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

gs.ws1 <- effects::predictorEffect('gs.ppt', m.ws.gs.cs.gs1.gs2  , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = gs.ppt)%>%
  mutate(var = 'gs.ppt')%>%
  merge(sig.ws.gs, by = 'ecosite')

sig.ws.cs <- emmeans::emtrends( m.ws.gs.cs.gs1.gs2 , ~ecosite, var = 'cs.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

cs.ws1 <- effects::predictorEffect('cs.ppt',m.ws.gs.cs.gs1.gs2  , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = cs.ppt)%>%
  mutate(var = 'cs.ppt')%>%
  merge(sig.ws.cs, by = 'ecosite')


sig.ws.lgs <- emmeans::emtrends( m.ws.gs.cs.gs1.gs2, ~ecosite, var = 'lag.gs') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

lgs.ws1 <- effects::predictorEffect('lag.gs',m.ws.gs.cs.gs1.gs2  , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = lag.gs)%>%
  mutate(var = 'lag.gs')%>%
  merge(sig.ws.lgs, by = 'ecosite')


sig.ws.lgs2 <- emmeans::emtrends( m.ws.gs.gs1.gs2, ~ecosite, var = 'lag2.gs') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

lgs2.ws1 <- effects::predictorEffect('lag2.gs', m.ws.gs.gs1.gs2  , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = lag2.gs)%>%
  mutate(var = 'lag2.gs')%>%
  merge(sig.ws.lgs2, by = 'ecosite')

ws1.df <- bind_rows(gs.ws1, lgs.ws1)%>%
  bind_rows(lgs2.ws1)%>%
  bind_rows(cs.ws1)

ws.slp.plt <- ggplot(ws1.df, aes(x = var.value, y = fit, color = ecosite))+
  geom_line(size = 2)+
  geom_ribbon(aes(x = var.value, ymin = lower, ymax = upper, fill = ecosite), alpha = 0.3)+
  labs(x = 'Precipitation mm', y = bquote('Predicted ANPP kg ha'^1), title = 'Warm Season', fill = 'Ecosite', color = 'Ecosite')+
  guides(linetype = 'none', fill = 'none', color = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c( 'lag2.gs','lag.gs', 'cs.ppt', 'gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))

var_names


#forest plot for warm season
?get_model_data()

ws.m.data <- get_model_data(m.ws.gs.cs.gs1.gs2, type = 'est')%>%
  filter(term %notin% c('ecositeSalt Flats', 'ecositeSandy'))%>%
  separate(term, sep = ':', into = c('ecosite', 'precip.param'), fill = 'left')%>%
  separate(ecosite, sep = 7, into = c('x','ecosite'), fill = 'right')%>%
  mutate(ecosite = ifelse(is.na(ecosite), 'Loamy', ecosite),
         precip.param = forcats::fct_relevel(precip.param, c('lag2.gs', 'lag.gs', 'cs.ppt', 'gs.ppt')))%>%
  group_by(precip.param)%>%
  mutate(estimate = ifelse(ecosite == 'Loamy', estimate, estimate + estimate[ecosite == 'Loamy']),
         conf.low = ifelse(ecosite == 'Loamy', conf.low, conf.low + estimate[ecosite == 'Loamy']),
         conf.high = ifelse(ecosite == 'Loamy', conf.high, conf.high + estimate[ecosite == 'Loamy']))



ws.m.plt <- ggplot(data = ws.m.data)+
  geom_point(aes(y = estimate, x = precip.param, fill = ecosite), size = 2, shape = 21, position = position_dodge(0.5))+
  geom_segment(aes(y = conf.low, yend = conf.high, x = precip.param, color = ecosite), position = position_dodge(0.5))+
  geom_hline(aes(yintercept = 0))+
  geom_vline(xintercept = c(1.45, 2.5, 3.55))+
  scale_fill_manual(values = precip.palette)+
  scale_color_manual(values = precip.palette)+
  labs(x = 'Precipitation parameters', y = 'Estimates', fill = 'Ecosite', color = 'Ecosite', title = 'Warm season graminoid mass')+
  theme_bw(base_size = 11)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())


layout <- c(
  'aaaa
  bbbb'
)

ws.m.plt + ws.slp.plt + plot_layout(design = layout)


#plot best model effect plot
#m.ws.gs.cs.gs1.gs2 = best warm season model
require(sjPlot)
require(broom)


## cool season 
aic.c3.df 
#m.c3.gs.cs.f.gs1.f1

sig.c3.gs <- emmeans::emtrends( m.c3.gs.cs.f.f1 , ~ecosite, var = 'gs.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

gs.c31 <- effects::predictorEffect('gs.ppt', m.c3.gs.cs.f.f1   , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = gs.ppt)%>%
  mutate(var = 'gs.ppt')%>%
  merge(sig.c3.gs, by = 'ecosite')

sig.c3.cs <- emmeans::emtrends( m.c3.gs.cs.f.f1  , ~ecosite, var = 'cs.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

cs.c31 <- effects::predictorEffect('cs.ppt',m.c3.gs.cs.f.f1   , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = cs.ppt)%>%
  mutate(var = 'cs.ppt')%>%
  merge(sig.c3.cs, by = 'ecosite')

sig.c3.f <- emmeans::emtrends( m.c3.gs.cs.f.f1  , ~ecosite, var = 'f.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

f.c31 <- effects::predictorEffect('f.ppt',m.c3.gs.cs.f.f1   , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = f.ppt)%>%
  mutate(var = 'f.ppt')%>%
  merge(sig.c3.f, by = 'ecosite')




sig.c3.lf<- emmeans::emtrends(m.c3.gs.cs.f.f1 , ~ecosite, var = 'lag.f') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

lf.c31 <- effects::predictorEffect('lag.f',m.c3.gs.cs.f.f1  , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = lag.f)%>%
  mutate(var = 'lag.f')%>%
  merge(sig.c3.lf, by = 'ecosite')

c31.df <- bind_rows(gs.c31, lf.c31)%>%
  bind_rows(cs.c31)%>%
  bind_rows(f.c31)

cs.slp.plt <- ggplot(c31.df, aes(x = var.value, y = fit, color = ecosite))+
  geom_ribbon(aes(x = var.value, ymin = lower, ymax = upper, fill = ecosite), alpha = 0.3)+
  labs(x = 'Precipitation mm', y = bquote('Predicted ANPP kg ha'^1), title = 'Cool Season', fill = 'Ecosite', color = 'Ecosite')+
  geom_line(size = 2)+
  guides(linetype = 'none', fill = 'none', color = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c('lag.f',  'f.ppt', 'cs.ppt','gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))


#forest plot for warm season


cs.m.data <- get_model_data(m.c3.gs.cs.f.f1, type = 'est')%>%
  filter(term %notin% c('ecositeSalt Flats', 'ecositeSandy'))%>%
  separate(term, sep = ':', into = c('ecosite', 'precip.param'), fill = 'left')%>%
  separate(ecosite, sep = 7, into = c('x','ecosite'), fill = 'right')%>%
  mutate(ecosite = ifelse(is.na(ecosite), 'Loamy', ecosite),
         precip.param = forcats::fct_relevel(precip.param, c('gs.ppt', 'cs.ppt', 'f.ppt', 'lag.f')))%>%
  group_by(precip.param)%>%
  mutate(estimate = ifelse(ecosite == 'Loamy', estimate, estimate + estimate[ecosite == 'Loamy']),
         conf.low = ifelse(ecosite == 'Loamy', conf.low, conf.low + estimate[ecosite == 'Loamy']),
         conf.high = ifelse(ecosite == 'Loamy', conf.high, conf.high + estimate[ecosite == 'Loamy']))%>%
  arrange(precip.param)



cs.m.plt <- ggplot(data = cs.m.data)+
  geom_point(aes(y = estimate, x = precip.param, fill = ecosite), size = 2, shape = 21, position = position_dodge(0.5))+
  geom_segment(aes(y = conf.low, yend = conf.high, x = precip.param, color = ecosite), position = position_dodge(0.5))+
  geom_hline(aes(yintercept = 0))+
  geom_vline(xintercept = c(1.45, 2.5, 3.55))+
  scale_fill_manual(values = precip.palette)+
  scale_color_manual(values = precip.palette)+
  labs(x = 'Precipitation parameters', y = 'Estimates', fill = 'Ecosite', color = 'Ecosite', title = "Cool season graminoid mass")+
  theme_bw(base_size = 11)+
  theme( axis.title.x=element_blank(),
         axis.text.x=element_blank())


layout <- c(
  'aaaa
  bbbb'
)

cs.m.plt + cs.slp.plt + plot_layout(design = layout)



## forb 
aic.fo.df 


sig.fo.gs <- emmeans::emtrends( m.fo.gs.cs.f.f1  , ~ecosite, var = 'gs.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

gs.fo1 <- effects::predictorEffect('gs.ppt', m.fo.gs.cs.f.f1    , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = gs.ppt)%>%
  mutate(var = 'gs.ppt')%>%
  merge(sig.fo.gs, by = 'ecosite')

sig.fo.cs <- emmeans::emtrends(m.fo.gs.cs.f.f1  , ~ecosite, var = 'cs.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

cs.fo1 <- effects::predictorEffect('cs.ppt',m.fo.gs.cs.f.f1    , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = cs.ppt)%>%
  mutate(var = 'cs.ppt')%>%
  merge(sig.fo.cs, by = 'ecosite')

sig.fo.f <- emmeans::emtrends( m.fo.gs.cs.f.f1   , ~ecosite, var = 'f.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

f.fo1 <- effects::predictorEffect('f.ppt', m.fo.gs.cs.f.f1    , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = f.ppt)%>%
  mutate(var = 'f.ppt')%>%
  merge(sig.fo.f, by = 'ecosite')


sig.fo.lf<- emmeans::emtrends( m.fo.gs.cs.f.f1  , ~ecosite, var = 'lag.f') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

lf.fo1 <- effects::predictorEffect('lag.f', m.fo.gs.cs.f.f1    , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = lag.f)%>%
  mutate(var = 'lag.f')%>%
  merge(sig.fo.lf, by = 'ecosite')

fo1.df <- bind_rows(gs.fo1, lf.fo1)%>%
  bind_rows(cs.fo1)%>%
  bind_rows(f.fo1)

fo.slp.plt <- ggplot(fo1.df, aes(x = var.value, y = fit, color = ecosite))+
  geom_ribbon(aes(x = var.value, ymin = lower, ymax = upper, fill = ecosite), alpha = 0.3)+
  labs(x = 'Precipitation mm', y = bquote('Predicted ANPP kg ha'^-1), title = 'Warm Season', fill = 'Ecosite', color = 'Ecosite')+
  geom_line(  size = 2)+
  
  guides(linetype = 'none', fill = 'none', color = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c('lag.f',  'f.ppt', 'cs.ppt','gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))


#forest plot for warm season


fo.m.data <- get_model_data(m.fo.gs.cs.f.f1, type = 'est')%>%
  filter(term %notin% c('ecositeSalt Flats', 'ecositeSandy'))%>%
  separate(term, sep = ':', into = c('ecosite', 'precip.param'), fill = 'left')%>%
  separate(ecosite, sep = 7, into = c('x','ecosite'), fill = 'right')%>%
  mutate(ecosite = ifelse(is.na(ecosite), 'Loamy', ecosite),
         precip.param = forcats::fct_relevel(precip.param, c('gs.ppt', 'cs.ppt', 'f.ppt', 'lag.f')))%>%
  group_by(precip.param)%>%
  mutate(estimate = ifelse(ecosite == 'Loamy', estimate, estimate + estimate[ecosite == 'Loamy']),
         conf.low = ifelse(ecosite == 'Loamy', conf.low, conf.low + estimate[ecosite == 'Loamy']),
         conf.high = ifelse(ecosite == 'Loamy', conf.high, conf.high + estimate[ecosite == 'Loamy']))%>%
  arrange(precip.param)



fo.m.plt <- ggplot(data = fo.m.data)+
  geom_point(aes(y = estimate, x = precip.param, fill = ecosite), size = 2, shape = 21, position = position_dodge(0.5))+
  geom_segment(aes(y = conf.low, yend = conf.high, x = precip.param, color = ecosite), position = position_dodge(0.5))+
  geom_hline(aes(yintercept = 0))+
  geom_vline(xintercept = c(1.45, 2.5, 3.55))+
  scale_fill_manual(values = precip.palette)+
  scale_color_manual(values = precip.palette)+
  labs(x = 'Precipitation parameters', y = 'Estimates', fill = 'Ecosite', color = 'Ecosite', title = 'Forb mass')+
  theme_bw(base_size = 11)+
  theme( axis.title.x=element_blank(),
         axis.text.x=element_blank())


layout <- c(
  'aaaa
  bbbb'
)

fo.m.plt + fo.slp.plt + plot_layout(design = layout)



### total mass : Fall 

aic.tm.df
sig.tm.gs <- emmeans::emtrends( m.tm.gs.cs.f.gs1.gs2 , ~ecosite, var = 'gs.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

gs.tm1 <- effects::predictorEffect('gs.ppt', m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = gs.ppt)%>%
  mutate(var = 'gs.ppt')%>%
  merge(sig.tm.gs, by = 'ecosite')

sig.tm.cs <- emmeans::emtrends( m.tm.gs.cs.f.gs1.gs2 , ~ecosite, var = 'cs.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

cs.tm1 <- effects::predictorEffect('cs.ppt', m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = cs.ppt)%>%
  mutate(var = 'cs.ppt')%>%
  merge(sig.tm.cs, by = 'ecosite')

sig.tm.f <- emmeans::emtrends( m.tm.gs.cs.f.gs1.gs2 , ~ecosite, var = 'f.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

f.tm1 <- effects::predictorEffect('f.ppt', m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = f.ppt)%>%
  mutate(var = 'f.ppt')%>%
  merge(sig.tm.f, by = 'ecosite')

sig.tm.lgs <- emmeans::emtrends(m.tm.gs.cs.f.gs1.gs2, ~ecosite, var = 'lag.gs') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

lgs.tm1 <- effects::predictorEffect('lag.gs',m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = lag.gs)%>%
  mutate(var = 'lag.gs')%>%
  merge(sig.tm.lgs, by = 'ecosite')


sig.tm.lgs2<- emmeans::emtrends( m.tm.gs.cs.f.gs1.gs2, ~ecosite, var = 'lag2.gs') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

lgs2.tm1 <- effects::predictorEffect('lag2.gs', m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = lag2.gs)%>%
  mutate(var = 'lag2.gs')%>%
  merge(sig.tm.lgs2, by = 'ecosite')

tm1.df <- bind_rows(gs.tm1, lgs.tm1)%>%
  bind_rows(lgs2.tm1)%>%
  bind_rows(cs.tm1)%>%
  bind_rows(f.tm1)

ggplot(tm1.df, aes(x = var.value, y = fit, color = ecosite))+
  geom_ribbon(aes(x = var.value, ymin = lower, ymax = upper, fill = ecosite), alpha = 0.3)+
  labs(x = 'ppt', y = 'predicted cool season graminoid mass')+
  geom_line( aes(linetype = forcats::fct_relevel(sig, c('y','n'))), size = 2)+
  
  guides(linetype = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 12)+
  facet_wrap(~var, ncol = 2, scales = 'free_x')


### total mass: gs2
sig.tm.gs <- emmeans::emtrends( m.tm.gs.cs.f.gs1.gs2 , ~ecosite, var = 'gs.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

gs.tm1 <- effects::predictorEffect('gs.ppt', m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = gs.ppt)%>%
  mutate(var = 'gs.ppt')%>%
  merge(sig.tm.gs, by = 'ecosite')

sig.tm.cs <- emmeans::emtrends( m.tm.gs.cs.f.gs1.gs2 , ~ ecosite, var = 'cs.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

cs.tm1 <- effects::predictorEffect('cs.ppt',m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = cs.ppt)%>%
  mutate(var = 'cs.ppt')%>%
  merge(sig.tm.cs, by = 'ecosite')

sig.tm.f <- emmeans::emtrends( m.tm.gs.cs.f.gs1.gs2 , ~ecosite, var = 'f.ppt') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

f.tm1 <- effects::predictorEffect('f.ppt', m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = f.ppt)%>%
  mutate(var = 'f.ppt')%>%
  merge(sig.tm.f, by = 'ecosite')

sig.tm.lgs <- emmeans::emtrends( m.tm.gs.cs.f.gs1.gs2, ~ecosite, var = 'lag.gs') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

lgs.tm1 <- effects::predictorEffect('lag.gs',m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = lag.gs)%>%
  mutate(var = 'lag.gs')%>%
  merge(sig.tm.lgs, by = 'ecosite')


sig.tm.lgs2<- emmeans::emtrends( m.tm.gs.cs.f.gs1.gs2, ~ecosite, var = 'lag2.gs') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

lgs2.tm1 <- effects::predictorEffect('lag2.gs', m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = lag2.gs)%>%
  mutate(var = 'lag2.gs')%>%
  merge(sig.tm.lgs2, by = 'ecosite')

tm1.df <- bind_rows(gs.tm1, lgs.tm1)%>%
  bind_rows(f.tm1)%>%
  bind_rows(cs.tm1)%>%
  bind_rows(lgs2.tm1)

tm.slp.plt <- ggplot(tm1.df, aes(x = var.value, y = fit, color = ecosite))+
  geom_ribbon(aes(x = var.value, ymin = lower, ymax = upper, fill = ecosite), alpha = 0.3)+
  labs(x = 'Precipitation mm', y = bquote('Predicted ANPP kg ha'^-1), title = 'Total', fill = 'Ecosite', color = 'Ecosite')+
  geom_line(  size = 2)+
  guides(linetype = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c('lag2.gs',  'lag.gs', 'cs.ppt',  'f.ppt','gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))


#forest plot for total mass 


tm.m.data <- get_model_data(m.tm.gs.cs.f.gs1.gs2, type = 'est')%>%
  filter(term %notin% c('ecositeSalt Flats', 'ecositeSandy'))%>%
  separate(term, sep = ':', into = c('ecosite', 'precip.param'), fill = 'left')%>%
  separate(ecosite, sep = 7, into = c('x','ecosite'), fill = 'right')%>%
  mutate(ecosite = ifelse(is.na(ecosite), 'Loamy', ecosite),
         precip.param = forcats::fct_relevel(precip.param, c('gs.ppt', 'f.ppt', 'cs.ppt', 'lag.gs', 'lag2.gs')))%>%
  group_by(precip.param)%>%
  mutate(estimate = ifelse(ecosite == 'Loamy', estimate, estimate + estimate[ecosite == 'Loamy']),
         conf.low = ifelse(ecosite == 'Loamy', conf.low, conf.low + estimate[ecosite == 'Loamy']),
         conf.high = ifelse(ecosite == 'Loamy', conf.high, conf.high + estimate[ecosite == 'Loamy']))




tm.m.plt <- ggplot(data = tm.m.data)+
  geom_point(aes(y = estimate, x = precip.param, fill = ecosite), size = 2, shape = 21, position = position_dodge(0.5))+
  geom_segment(aes(y = conf.low, yend = conf.high, x = precip.param, color = ecosite), position = position_dodge(0.5))+
  geom_hline(aes(yintercept = 0))+
  geom_vline(xintercept = c(1.45, 2.5, 3.55,4.6))+
  scale_fill_manual(values = precip.palette)+
  scale_color_manual(values = precip.palette)+
  labs(x = 'Precipitation parameters', y = 'Estimates', fill = 'Ecosite', color = 'Ecosite', title = 'Total Mass')+
  theme_bw(base_size = 11)+
  theme( axis.title.x=element_blank(),
         axis.text.x=element_blank())

tm.m.plt + tm.slp.plt + plot_layout(design = layout)

layout <- c(
  'aaaa
  aaaa
  bbbb
  bbbb
  cccc
  cccc
  dddd
  dddd'
)



###combine all plots for figure 4
fig.4 <-  tm.slp.plt + ws.slp.plt  + cs.slp.plt + fo.slp.plt + 
  plot_layout( ncol = 2, guides = 'collect', axis_titles = 'collect') +
  plot_annotation(tag_levels = 'A')

fig.4

#ggsave( 'model_figure_4.png',fig.4, width = 9, height = 14.5)
### plot spatiotemporal relationships for each functional group and gs precip #####



veg.precip.df$water.yr.f <-  as.factor(veg.precip.df$water.yr)

#make global sensitivity model plot - forest plot
m.global.t.mass <- lm(data = veg.precip.df, total.mass.dm ~ gs.ppt*ecosite )
m.global.wspg <- lm(data = veg.precip.df, wspg.dm ~ gs.ppt*ecosite)
m.global.c3pg <- lm(data = veg.precip.df, c3pg.dm ~ gs.ppt*ecosite)
m.global.forb <- lm(data = veg.precip.df, forb.dm ~ gs.ppt*ecosite)


slps.global <- as.data.frame(emmeans::emtrends(m.global.t.mass, var = c('gs.ppt'), ~ecosite))%>%
  mutate(response = 't.mass')%>%
  bind_rows(as.data.frame(emmeans::emtrends(m.global.wspg, var = c('gs.ppt'), ~ecosite))%>%
              mutate(response = 'wspg'))%>%
  bind_rows(as.data.frame(emmeans::emtrends(m.global.c3pg, var = c('gs.ppt'), ~ecosite))%>%
                mutate(response = 'c3pg'))%>%
  bind_rows(as.data.frame(emmeans::emtrends(m.global.forb, var = c('gs.ppt'), ~ecosite))%>%
              mutate(response = 'forb'))
slps.global%>%
  filter(ecosite == 'Loamy')%>%
  ggplot( aes(x = fct_relevel(response, c('t.mass', 'wspg', 'c3pg', 'forb')), fill = ecosite))+
  geom_segment(aes(y = lower.CL, yend = upper.CL), position = position_dodge(0.5), linewidth = 1)+
  geom_point(aes(y = gs.ppt.trend), position = position_dodge(0.5), shape = 21, size = 4)+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw()+
  labs(x = 'Functional group', y = 'NPP Sensitivity')+
  scale_x_discrete(labels = c('Total', "Warm season", "Cool season", 'Forb'))+
  theme(text = element_text(size = 15))

ggplot(data = slps.global, aes(x = fct_relevel(response, c('t.mass', 'wspg', 'c3pg', 'forb')), fill = ecosite))+
  geom_segment(aes(y = lower.CL, yend = upper.CL), position = position_dodge(0.5), linewidth = 1)+
  geom_point(aes(y = gs.ppt.trend), position = position_dodge(0.5), shape = 21, size = 4)+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw()+
  labs(x = 'Functional group', y = bquote(ANPP~Sensitivity~kg~ha^-1~mm^-1), fill = "Ecosite")+
  scale_x_discrete(labels = c('Total', "Warm season", "Cool season", 'Forb'))+
  theme(text = element_text(size = 12))

#ggsave('sens_gsppt.png', width = 5.4, height = 3)



## Breakdown of biomass percentage by ecosite Supplemental figure
veg.precip.df%>%
  dplyr::select(pairblock, ecosite, water.yr, wspg, c3pg, forb)%>%
  pivot_longer(-c(pairblock, ecosite, water.yr), values_to = 'mass', names_to = 'fun.group')%>%
  ggplot()+
  stat_summary(fun = mean , geom = 'col', aes(x = water.yr, y = mass, fill = fun.group), position = position_stack())+
  facet_wrap(~ecosite, ncol = 1)+
  scale_fill_manual(values = precip.palette)+
  theme_bw()+
  labs(fill = 'Functional Group', y = 'Mass', x = 'Harvest year')
    
    