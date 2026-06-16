# precipitation and anpp variation 
#figures 2 and 3 

`%notin%` <- Negate(`%in%`) 
#library
library(tidyverse)
library(sf)
library(patchwork)
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


#CPER pasture shapefile
past.sf <- st_read('data/pastures.shp')

#precip

#precip at the plotpair level
precip.vars <- read.csv('data/CPER_precip_variables.csv')

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

###########  CV analysis + figure 3 ##########

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
  #st_drop_geometry()%>%
  group_by(id)%>%
  summarise(cv = sd(gs.ppt, na.rm = T)/mean(gs.ppt, na.rm = T),
            mgsp = mean(gs.ppt, na.rm = T),
            variation = 'Temporal')

cv.cross <- temp.cv%>%
  st_drop_geometry()%>%
  dplyr::select(cv, variation)%>%
  bind_rows(spatial.cv)


mean(ppt.vars.sf$f.ppt, na.rm = T)
mean(ppt.vars.sf$cs.ppt, na.rm = T)
mean(ppt.vars.sf$gs.ppt, na.rm = T)


#cv of ANPP
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
?annotate()
mass.cv.plt <- ggplot(cv.mass.cross, aes(x = cv))+
  geom_density( aes(x = cv, fill = variation), alpha = 0.75)+
  xlim(0,1)+
  labs(y = 'Density', x = 'CV', fill = 'Variation', color = 'Variation')+
  geom_vline(data = mass.cv.means, aes(xintercept = cv, color = variation), linewidth = 0.75)+
  ylim(0,25)+
  annotate('text', label = 'ANPP', x = 0.75, y = 10, size = 3)+
  annotate('text', x = -Inf, y = Inf, label = 'G',
           hjust = -1, vjust = 1.5)+
  scale_fill_manual(values = c('#008080', 'tomato2'))+
  scale_color_manual(values = c('#008080', 'tomato2'))+
  theme_bw(base_size = 11)+
  theme(text=element_text(family="sans"),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.margin =  margin(0, 0, 8, 0, "pt"))

precip.cv.plt <- ggplot(cv.cross, aes(x = cv))+
  geom_density( aes(x = cv, fill = variation), alpha = 0.75)+
  annotate('text', label = 'Precipitation', x = 0.7, y = 10, size = 3)+
  geom_vline(data = precip.cv.means, aes(xintercept = cv, color = variation), linewidth = 0.75)+
  scale_fill_manual(values = c('#008080', 'tomato2'))+
  scale_color_manual(values = c('#008080', 'tomato2'))+
  xlim(0,1)+
  annotate('text', x = -Inf, y = Inf, label = 'F',
           hjust = -1, vjust = 1.5)+
  labs(y = 'Density', x = 'CV', fill = 'Variation', color = 'Variation')+
  theme_bw(base_size = 11)+
  theme(text=element_text(family="sans"),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.margin =  margin(0, 0, 8, 0, "pt"))

# timeline of gs precip
timeline.precip.plt <- ppt.vars.sf%>%
  ggplot(aes(x = water.yr, y = gs.ppt, group = id,  color = id))+
  geom_line()+
  labs(x = 'Water year', y = 'GS PPT (mm)')+
  annotate('text', x = -Inf, y = Inf, label = 'A',
           hjust = -1, vjust = 1.5)+
  scale_x_continuous(breaks = seq(2012, 2022, 2))+
  scale_color_viridis_d()+
  guides(color = 'none')+
  theme_bw(base_size = 11)+
  theme(text=element_text(family="sans"), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.margin =  margin(0, 0, 8, 0, "pt"))



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
  stat_summary(geom = 'pointrange', fun.data = mean_cl_normal, aes(x = as.numeric(as.character(water.yr)), y = mass, group = ecosite, color = ecosite), size = 0.4)+
  stat_summary(geom = 'line', fun.data = mean_cl_normal, aes(x = as.numeric(as.character(water.yr)), y = mass, group = ecosite, color = ecosite))+
  geom_text(data = timeline.veg.labels[2:4,], aes(x = x, y = y, label = label),
            hjust = -1, vjust = 1.5)+
  facet_grid(rows = vars(forcats::fct_relevel(fun.grp, c( 'wspg', 'c3pg', 'forb'))), labeller = as_labeller(fun.grp.labels))+
  labs(x = 'Water year', y =   bquote('ANPP (kg ha'^-1*')'), color = 'Ecosite')+
  scale_x_continuous(breaks = seq(2010, 2022, 2), limits = c(2011, 2023))+
  scale_color_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  theme(text=element_text(family="sans"),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        #axis.title = element_text(vjust = 10),
        legend.position = 'bottom',
        plot.margin =  margin(0, 0, 8, 0, "pt"))

timeline.t.mass  <- veg.precip.df%>%
  dplyr::select(total.mass, wspg, c3pg, forb, water.yr, pairblock, ecosite)%>%
  pivot_longer(-c(water.yr, pairblock, ecosite), names_to = 'fun.grp', values_to = 'mass')%>%
  filter(fun.grp == 'total.mass')%>%
  # mutate(fun.grp = forcats::fct_relevel(fun.grp, c( 'wspg', 'c3pg', 'forb')))%>%
  ggplot()+
  # geom_line(aes(x = as.numeric(as.character(water.yr)), y = mass, group = pairblock, color = ecosite), linewidth = 0.7)+
  stat_summary(geom = 'pointrange', fun.data = mean_cl_normal, aes(x = as.numeric(as.character(water.yr)), y = mass, group = ecosite, color = ecosite), size = 0.4)+
  stat_summary(geom = 'line', fun.data = mean_cl_normal, aes(x = as.numeric(as.character(water.yr)), y = mass, group = ecosite, color = ecosite))+
  geom_text(data = timeline.veg.labels[1,], aes(x = x, y = y, label = label),
            hjust = -1, vjust = 1.5)+
  facet_grid(rows = vars(fun.grp), labeller = as_labeller(fun.grp.labels), scale = 'free_y')+
  labs(x = 'Water year', y = bquote('ANPP (kg ha'^-1*')'), color = 'Ecosite')+
  scale_x_continuous(breaks = seq(2010, 2022, 2), limits = c(2011, 2023))+
  scale_color_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  theme(text=element_text(family="sans"),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        #axis.title = element_text(vjust = 10),
        legend.position = 'none',
        plot.margin =  margin(0, 0, 8, 0, "pt"))

cvplt <- (precip.cv.plt / mass.cv.plt / guide_area() /  plot_spacer()) + plot_layout(nrow = 4, ncol = 1,  axes = 'collect', guides = 'collect', heights = c(2,2,1,1))

layout1 <- '
AC
BC
DD'


fig.3.1 <- timeline.precip.plt + timeline.t.mass + timeline.veg + guide_area() + plot_layout(design = layout1, heights = c(3,3,1), guides = 'collect')


fig.3 <- wrap_elements(fig.3.1) + wrap_elements(cvplt ) +  plot_layout(ncol = 2, widths = c(0.75, 0.25)) 


ggsave('figures/figure_3.svg', fig.3, height = 12, width = 20, units = 'cm')
