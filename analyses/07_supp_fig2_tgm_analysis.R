####generate manuscript figs
#TKM
#5.6.24

`%notin%` <- Negate(`%in%`) 
#library

library(tidyverse)
library(patchwork)
library(lme4)
require(sjPlot)

#palettes
precip.palette <- c(  '#F78154','#389DE5', '#5DDCAC' ,'#9368B7')


#-------------------
#####set up dataframe 

# #precip at the plot level
# precip.vars <- read.csv('data/CPER_precip_variables.csv')
# 
# biomass <- read.csv('data/AGM_Biomass_Widecln_attr_2024-03-26.csv')%>%
#   rename_all(tolower)%>%
#   filter(is.na(x) == F)%>%
#   mutate(year = yearsampled)%>%
#   dplyr::select(year, pasture, plot, transect, ag, bobu, c3pg, forb, wspg, sd, ss, x,y, ecosite)%>%
#   pivot_longer(-c(year, pasture, plot, transect,x,y, ecosite), names_to = 'funct_grp', values_to = 'mass.g')%>%
#   group_by(pasture, plot, year, funct_grp,x, y, ecosite)%>%
#   summarize(m.mass = mean(mass.g))%>%
#   ungroup()%>%
#   group_by(year, plot, pasture, x, y, ecosite)%>%
#   summarise(wspg = sum(m.mass[funct_grp %in% c('wspg', 'bobu')]),
#             c3pg = sum(m.mass[funct_grp %in% 'c3pg']),
#             ag = sum(m.mass[funct_grp %in% 'ag']),
#             sd = sum(m.mass[funct_grp %in% c('sd')]),
#             forb = sum(m.mass[funct_grp %in% c('forb')]),
#             total.mass = sum(m.mass[funct_grp %notin% c('sd')]))%>%
#   rename(water.yr = year)
# 
# veg.precip.separated <- left_join(precip.vars, biomass, by = c('water.yr', 'plot', 'pasture'))
# write.csv(veg.precip.separated,'data/legacy_separated_df.csv', row.names = FALSE)

veg.precip.separated <- read.csv('data/legacy_separated_df.csv')
past_treats <- read.csv('data/CPER_treatments')

past_trm <- past_treats%>%
  filter(treatment %in% 'TGM')

#filter to TRM treatments 
veg.precip.df.trm <- veg.precip.separated%>%
  filter(pasture %in% past_trm$pasture)

#Run models 
#model total mass######

m.tm.gs.cs.f.gs1.gs2 <- lmer(total.mass ~ ecosite*gs.ppt + ecosite*cs.ppt + ecosite*f.ppt + ecosite*lag.gs  + ecosite*lag2.gs + (1|plotpasture) + (1|water.yr),
                             data = veg.precip.df.trm)

#model checks on best model 
plot(m.tm.gs.cs.f.gs1.gs2)
#qqmath(resid(m.tm.gs.cs.f.gs1.gs2))
plot(m.tm.gs.cs.f.gs1.gs2, rstudent(.) ~ hatvalues(.))


#write.csv(aic.tm.df, 'total_mass_AIC.csv')

#model warm season ######

m.ws.gs.cs.gs1.gs2 <- lmer(wspg ~ ecosite*gs.ppt + ecosite*cs.ppt + ecosite*lag.gs + ecosite*lag2.gs+ (1|plotpasture) + (1|water.yr),
                           data = veg.precip.df.trm)

#



#model checks on best model 
plot(m.ws.gs.cs.gs1.gs2 )
#qqmath(resid(m.ws.gs.cs.gs1.gs2 ))
plot(m.ws.gs.cs.gs1.gs2 , rstudent(.) ~ hatvalues(.))



#model cool season ######


m.c3.gs.cs.f.f1 <- lmer(c3pg ~ ecosite*gs.ppt + ecosite*cs.ppt + ecosite*f.ppt + ecosite*lag.f+ (1|plotpasture) + (1|water.yr),
                        data = veg.precip.df.trm)


#model checks on best model 
plot(m.c3.gs.cs.f.f1   )
qqmath(resid(m.c3.gs.cs.f.f1  ))
plot(m.c3.gs.cs.f.f1   , rstudent(.) ~ hatvalues(.))

#model forbs ######


m.fo.gs.cs.f.f1 <- lmer(forb ~ ecosite*gs.ppt + ecosite*cs.ppt + ecosite*f.ppt + ecosite*lag.f+ (1|plotpasture) + (1|water.yr),
                        data = veg.precip.df.trm)

#model checks on best model 
plot(m.fo.gs.cs.f.f1   )

plot(m.fo.gs.cs.f.f1   , rstudent(.) ~ hatvalues(.))



#### table of AIC####
models.tbl <- cbind(aic.tm.df, aic.ws.df, aic.c3.df, aic.fo.df)

write.csv(models.tbl, 'all_models_trm.csv')


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

#####generate plot 
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
  labs(x = 'Precipitation (mm)', y = bquote('Predicted ANPP (kg ha'^-1*')'), title = 'Warm season', fill = 'Ecosite', color = 'Ecosite')+
  guides(linetype = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c( 'lag2.gs','lag.gs', 'cs.ppt', 'gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))

var_names




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
  labs(x = 'Precipitation (mm)', y = bquote('Predicted ANPP (kg ha'^-1*')'), title = 'Cool season', fill = 'Ecosite', color = 'Ecosite')+
  geom_line(size = 2)+
  
  guides(fill = 'none', color = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c('lag.f',  'f.ppt', 'cs.ppt','gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))


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
  labs(x = 'Precipitation (mm)', y = bquote('Predicted ANPP (kg ha'^-1*')'), title = 'Forb', fill = 'Ecosite', color = 'Ecosite')+
  geom_line(  size = 2)+
  
  guides(fill = 'none', color = 'none')+
  
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c('lag.f',  'f.ppt', 'cs.ppt','gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))
 



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
  labs(x = 'Precipitation (mm)', y = bquote('Predicted ANPP (kg ha'^-1*')'), title = 'Total', fill = 'Ecosite', color = 'Ecosite')+
  geom_line(  size = 2)+
  
  guides(fill = 'none', color = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c('lag2.gs',  'lag.gs', 'cs.ppt',  'f.ppt','gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))



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
fig.s2 <-  tm.slp.plt + ws.slp.plt  + cs.slp.plt + fo.slp.plt + 
  plot_layout( design = layout, guides = 'collect', axes = 'collect') +
  plot_annotation(tag_levels = 'A')

fig.s2

#ggsave( 'model_figure_4_trm.png',fig.4, width = 7, height = 8)
### plot spatiotemporal relationships for each functional group and gs precip #####


require(ggpmisc)
global.wspg.plt <- ggplot(data = veg.precip.df, aes(x = gs.ppt, y = wspg.dm))+
  geom_point(aes(fill = as.factor(water.yr)), shape = 21, size = 1.3)+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = ..eq.label..), size = 3)+
  stat_poly_eq(aes(label =  ..rr.label..), label.y = 0.8, size = 3)+
  scale_fill_viridis_d(option = 'turbo')+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = 'Warm ssn gram', fill = 'Harvest year')+
  scale_x_continuous(limits = c(90, 380))+
  scale_y_continuous(limits = c(-800, 2000))+
  theme_bw(base_size = 11)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing.x= unit(0, 'lines'))

global.c3pg.plt <- ggplot(data = veg.precip.df, aes(x = gs.ppt, y = c3pg.dm))+
  geom_point(aes(fill = as.factor(water.yr)), shape = 21, size = 1.3)+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = ..eq.label..), size = 3)+
  stat_poly_eq(aes(label =  ..rr.label..), label.y = 0.8, size = 3)+
  scale_fill_viridis_d(option = 'turbo')+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = 'Cool ssn gram' , fill = 'Harvest year')+
  scale_x_continuous(limits = c(90, 380))+
  scale_y_continuous(limits = c(-800, 2000))+
  theme_bw(base_size = 11)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing.x= unit(0, 'lines'))

global.forb.plt <- ggplot(data = veg.precip.df, aes(x = gs.ppt, y = forb.dm))+
  geom_point(aes(fill = as.factor(water.yr)), shape = 21, size = 1.3)+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = paste(..eq.label..)), size = 3)+
  stat_poly_eq(aes(label =  ..rr.label..), label.y = 0.8, size = 3)+
  scale_fill_viridis_d(option = 'turbo')+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = 'Forb', fill = 'Harvest year')+
  scale_x_continuous(limits = c(90, 380))+
  scale_y_continuous(limits = c(-800, 2000))+
  theme_bw(base_size = 11)+
  theme(panel.spacing.x= unit(0, 'lines'))




global.tm.plt <- ggplot(data = veg.precip.df, aes(x = gs.ppt, y = total.mass.dm))+
  geom_point(aes(fill = as.factor(water.yr)), shape = 21, size = 1.3)+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = paste(..eq.label..)), size = 3)+
  stat_poly_eq(aes(label =  ..rr.label..), label.y = 0.8, size = 3)+
  scale_fill_viridis_d(option = 'turbo')+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = 'Total mass (kg*ha-1)', fill = 'Harvest year')+
  scale_x_continuous(limits = c(90, 380))+
  #scale_y_continuous(limits = c(-500, 1200))+
  theme_bw(base_size = 11)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing.x= unit(0, 'lines'))

sens_fig <- global.tm.plt + plot_spacer() + global.wspg.plt +plot_spacer() + global.c3pg.plt + plot_spacer() +global.forb.plt +
  plot_layout(guides = "collect")+
  plot_layout(axes = "collect_x", ncol = 1, heights = c(3,-.55,3,-.55,3,-.55,3)) 

ggsave( 'sens_fig.png', sens_fig, width = 6, height = 7)
x <- summary(lm(total.mass.dm ~ gs.ppt, data = veg.precip.df.trm))
x$coefficients[2,4]


##compare to non spatial model 
veg.precip.df <- veg.precip.df%>%
  group_by(water.yr)%>%
  mutate(gs.ppt.m = gs.ppt[pairblock == '1_1'])

non.spat.m <-summary(lmer(data = veg.precip.df, total.mass ~ gs.ppt.m*ecosite + (1|plotpasture)))

non.spat.rmse <- sqrt(mean(non.spat.m$residuals^2))

spat.m <-summary(lmer(data = veg.precip.df, total.mass ~ gs.ppt*ecosite + (1|pairblock)))

spat.rmse <- sqrt(mean(spat.m$residuals^2))

ggplot(data = veg.precip.df, aes(x = x, y = y))+
  geom_text(aes(label = pairblock))

global.tm.plt.nsp <- ggplot(data = veg.precip.df, aes(x = gs.ppt.m, y = total.mass.dm))+
  geom_point(aes(fill = as.factor(water.yr)), shape = 21, size = 3)+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "/")))+
  scale_fill_viridis_d(option = 'turbo')+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = 'Total mass (kg/ha)', fill = 'Harvest year')+
  scale_x_continuous(limits = c(90, 380))+
  theme_bw(base_size = 15)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#so spatial precipitation doesn't improve fit overall
#but does it in specific years?
ggplot(data = veg.precip.df, aes(x = gs.ppt, y = total.mass.dm))+
  geom_point(aes(fill = as.factor(water.yr)), shape = 21, size = 3)+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")))+
  scale_fill_viridis_d(option = 'turbo')+
  facet_wrap(~water.yr)+
  labs(x = 'Growing season PPT (mm)', y = 'Total mass (kg/ha)', fill = 'Harvest year')+
  scale_x_continuous(limits = c(90, 380))+
  theme_bw(base_size = 15)

summary(lm(data = veg.precip.df, total.mass.dm ~gs.ppt))$r.squared
#plot CV by R2 of spatial only model


annual.cv.df <- veg.precip.df%>%
  group_by(water.yr)%>%
  summarise(cv = sd(gs.ppt)/mean(gs.ppt),
            sd = sd(gs.ppt),
            sd2 = sd(gs.ppt + lag.gs),
            mean = mean(gs.ppt),
            lag.mean = mean(lag.gs),
            mean.rt = mean/lag.mean,
            rsq = as.numeric(summary(lm(gs.ppt ~ total.mass.dm))$r.squared[1]))

summary(lm(data = annual.cv.df, rsq ~ mean.rt))
require(ggpmisc)

sdspplt <- ggplot(annual.cv.df, aes(x = sd, y = rsq))+
  geom_point(shape = 21,size = 3, aes(fill = as.factor(water.yr)))+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")))+
  scale_fill_viridis_d(option = 'turbo')+
  labs(x = 'PPT standard deviation', y  = 'R squared', fill = 'Year')+
  theme_bw()

sd2plt <- ggplot(annual.cv.df, aes(x = sd2, y = rsq))+
  geom_point(shape = 21,size = 3, aes(fill = sd))+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")))+
  scale_fill_viridis_c(option = 'turbo')+
  labs(x = 'PPT standard deviation 2-yr', y  = 'R squared')+
  theme_bw()

cvspplt <- ggplot(annual.cv.df, aes(x = cv, y = rsq))+
  geom_point(shape = 21,size = 3, aes(fill = sd))+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")))+
  scale_fill_viridis_c(option = 'turbo')+
  labs(x = 'PPT coef of variaiton', y  = 'R squared')+
  theme_bw()

 cvspplt + sdspplt + plot_layout(guides = 'collect', axes = 'collect')

anppt <- ggplot(annual.cv.df, aes(x = mean, y = rsq))+
  geom_point(shape = 21,size = 3, aes(fill = sd))+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")))+
  scale_fill_viridis_c(option = 'turbo')+
  labs(x = 'Annual ppt', y  = 'R squared')+
  theme_bw()

pptrt <- ggplot(annual.cv.df, aes(x = mean.rt, y = rsq))+
  geom_point(shape = 21,size = 3, aes(fill = sd))+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")))+
  scale_fill_viridis_c(option = 'turbo')+
  labs(x = 'Current vs previous ppt', y  = 'R squared')+
  theme_bw()

prevppt <- ggplot(annual.cv.df, aes(x = lag.mean, y = rsq))+
  geom_point(shape = 21,size = 3, aes(fill = sd))+
  geom_smooth(method = 'lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")))+
  scale_fill_viridis_c(option = 'turbo')+
  labs(x = 'Previous ppt', y  = 'R squared')+
  theme_bw()



anppt + prevppt + pptrt + plot_layout(guides = 'collect', axes = 'collect')
  

  #2021 case study
require(ggpmisc)
### 2021 case study 
total.mass.21 <- veg.precip.df%>%
  filter(ecosite != 'Salt Flats')%>%
  filter(water.yr == 2021)%>%
  group_by(ecosite)%>%
  mutate(r.sq = summary(lm(total.mass.dm ~ gs.ppt))$r.squared,
         sig = ifelse(r.sq < 0.1, 'n', 'y'))%>%
  ggplot( aes(x = gs.ppt, y = total.mass.dm))+
  geom_point(shape = 21, size = 2, fill = 'grey50')+
  geom_smooth(method = 'lm', aes(linetype = sig), color = 'black')+
  scale_linetype_manual(values = c( 'solid'))+
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")))+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = bquote('Total mass kg ha'^-1))+
  guides(linetype = 'none')+
  theme_bw(base_size = 12)+
  theme(axis.title.x = element_blank())

t <- summary(lm(data = veg.precip.df, wspg ~ gs.ppt))


wspg.21 <- veg.precip.df%>%
  filter(ecosite != 'Salt Flats')%>%
  filter(water.yr == 2021)%>%
  group_by(ecosite)%>%
  mutate(r.sq = summary(lm(wspg.dm ~ gs.ppt))$r.squared,
         sig = ifelse(r.sq < 0.1, 'n', 'y'))%>%
  ggplot( aes(x = gs.ppt, y =wspg.dm))+
  geom_point(shape = 21, size = 2, fill = precip.palette[4])+
  geom_smooth(method = 'lm', aes(linetype = sig), color = 'black')+
  scale_linetype_manual(values = c( 'solid'))+
  stat_poly_eq(aes(label = paste( ..rr.label.., sep = "~~~")))+
  facet_wrap(~ecosite)+
  labs( y = bquote('WS mass kg ha'^-1))+
  guides(linetype = 'none')+
  theme_bw(base_size = 12)+
  theme(axis.title.x = element_blank())


c3pg.21 <- veg.precip.df%>%
  filter(ecosite != 'Salt Flats')%>%
  filter(water.yr == 2021)%>%
  group_by(ecosite)%>%
  mutate(p.val = summary(lm(c3pg.dm ~ gs.ppt))$coefficients[2,4],
         sig = ifelse(p.val > 0.1, 'n', 'y'))%>%
  ggplot( aes(x = gs.ppt, y = c3pg.dm))+
  geom_point(shape = 21, size = 3, fill =  precip.palette[4])+
  geom_smooth(method = 'lm', aes(linetype = sig))+
  scale_linetype_manual(values = c('dashed', 'solid'))+
  stat_poly_eq(aes(label = paste( ..rr.label.., ..p.value.label.., sep = "~~~")))+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = 'CS mass (kg/ha)')+
  theme_bw(base_size = 12)

forb.21 <- veg.precip.df%>%
  
  filter(water.yr == 2021)%>%
  filter(ecosite != 'Salt Flats')%>%
  group_by(ecosite)%>%
  mutate(p.val = summary(lm(forb.dm ~ gs.ppt))$coefficients[2,4],
         sig = ifelse(p.val > 0.1, 'n', 'y'))%>%
  ggplot( aes(x = gs.ppt, y = forb.dm))+
  geom_point(shape = 21, size = 3, fill=  precip.palette[3] )+
  geom_smooth(method = 'lm', aes(linetype = sig))+
  scale_linetype_manual(values = c('dashed', 'solid'))+
  stat_poly_eq(aes(label = paste( ..rr.label.., ..p.value.label.., sep = "~~~")))+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = 'Forb mass (kg/ha)')+
  theme_bw(base_size = 12)

total.mass.21 + wspg.21 + c3pg.21 + forb.21 + plot_layout(axes = 'collect', ncol = 1)


total.mass.22 <- veg.precip.df%>%
  filter(water.yr == 2022)%>%
  filter(ecosite != 'Salt Flats')%>%
  group_by(ecosite)%>%
  mutate(r.sq = summary(lm(t.mass.resid ~ lag2.gs))$r.squared,
         sig = ifelse(r.sq < 0.1, 'n', 'y'))%>%
  ggplot( aes(x = lag.gs, y = t.mass.resid))+
  geom_point(shape = 22, size = 2, fill = 'grey50')+
  geom_smooth(method = 'lm', aes(linetype = sig), color = 'black')+
  scale_linetype_manual(values = c('dashed', 'solid'))+
  stat_poly_eq(aes(label = paste( ..rr.label.., sep = "~~~")))+
  facet_wrap(~ecosite)+
  labs(x = '', y = bquote('Total mass residual kg ha'^-1))+
  guides(linetype = 'none')+
  theme_bw(base_size = 12)+
  theme(axis.title.x = element_blank())

wspg.22 <- veg.precip.df%>%
  filter(water.yr == 2022)%>%
  filter(ecosite != 'Salt Flats')%>%
  group_by(ecosite)%>%
  mutate(r.sq = summary(lm(wspg.resid~ lag2.gs))$r.squared,
         sig = ifelse(r.sq < 0.1, 'n', 'y'))%>%
  ggplot( aes(x = lag.gs, y =wspg.resid))+
  geom_point(shape = 22, size = 2, fill = precip.palette[4])+
  geom_smooth(method = 'lm', aes(linetype = sig), color = 'black')+
  scale_linetype_manual(values = c( 'solid'))+
  stat_poly_eq(aes(label = paste( ..rr.label.., sep = "~~~")))+
  facet_wrap(~ecosite)+
  labs(x = '', y = bquote('WS residual kg ha'^-1))+
  guides(linetype = 'none')+
  theme_bw(base_size = 12)+
  theme(axis.title.x = element_blank())

c3pg.22 <- veg.precip.df%>%
  filter(water.yr == 2022)%>%
  filter(ecosite != 'Salt Flats')%>%
  group_by(ecosite)%>%
  mutate(p.val = summary(lm(c3pg.resid~ lag.gs))$coefficients[2,4],
         sig = ifelse(p.val > 0.1, 'n', 'y'))%>%
  ggplot( aes(x = lag.gs, y = c3pg.resid))+
  geom_point(shape = 22, size = 3, fill = precip.palette[4])+
  geom_smooth(method = 'lm', aes(linetype = sig))+
  scale_linetype_manual(values = c('dashed', 'solid'))+
  stat_poly_eq(aes(label = paste( ..rr.label.., sep = "~~~")))+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = 'CS residual (kg/ha)')+
  theme_bw(base_size = 12)

forb.22 <- veg.precip.df%>%
  filter(water.yr == 2022)%>%
  filter(ecosite != 'Salt Flats')%>%
  group_by(ecosite)%>%
  mutate(p.val = summary(lm(forb.resid~ lag.gs))$coefficients[2,4],
         sig = ifelse(p.val > 0.1, 'n', 'y'))%>%
  ggplot( aes(x = lag.gs, y = forb.resid))+
  geom_point(shape = 22, size = 3, fill= precip.palette[3] )+
  geom_smooth(method = 'lm', aes(linetype = sig))+
  scale_linetype_manual(values = c('dashed', 'solid'))+
  stat_poly_eq(aes(label = paste( ..rr.label.., sep = "~~~")))+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = ' forb residual(kg/ha)')+
  theme_bw(base_size = 12)

total.mass.22 + wspg.22 + c3pg.22 + plot_layout(ncol = 1)




total.mass.23 <- veg.precip.df%>%
  filter(water.yr == 2023)%>%
  filter(ecosite != 'Salt Flats')%>%
  group_by(ecosite)%>%
  mutate(r.sq = summary(lm(t.mass.resid ~ lag2.gs))$r.squared,
         sig = ifelse(r.sq < 0.1, 'n', 'y'))%>%
  ggplot( aes(x = lag2.gs, y = t.mass.resid))+
  geom_point(shape = 23, size = 2, fill = 'grey50')+
  geom_smooth(method = 'lm', aes(linetype = sig), color = 'black')+
  scale_linetype_manual(values = c('dashed', 'solid'))+
  stat_poly_eq(aes(label = paste( ..rr.label.., sep = "~~~")))+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = bquote('Total mass residual kg ha'^-1))+
  guides(linetype = 'none')+
  theme_bw(base_size = 12)

wspg.23 <- veg.precip.df%>%
  filter(water.yr == 2023)%>%
  filter(ecosite != 'Salt Flats')%>%
  group_by(ecosite)%>%
  mutate(r.sq = summary(lm(wspg.resid~ lag2.gs))$r.squared,
         sig = ifelse(r.sq < 0.1, 'n', 'y'))%>%
  ggplot( aes(x = lag2.gs, y =wspg.resid))+
  geom_point(shape = 23, size = 2, fill = precip.palette[4])+
  geom_smooth(method = 'lm', aes(linetype = sig), color = 'black')+
  scale_linetype_manual(values = c('dashed', 'solid'),
                        labels = c('N', 'Y'))+
  stat_poly_eq(aes(label = paste( ..rr.label.., sep = "~~~")))+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = bquote('WS residual kg ha'^-1))+
  guides(linetype = 'none')+
  theme_bw(base_size = 12)

c3pg.23 <- veg.precip.df%>%
  filter(water.yr == 2023)%>%
  filter(ecosite != 'Salt Flats')%>%
  group_by(ecosite)%>%
  mutate(p.val = summary(lm(c3pg.resid~ lag2.gs))$coefficients[2,4],
         sig = ifelse(p.val > 0.1, 'n', 'y'))%>%
  ggplot( aes(x = lag2.gs, y = c3pg.resid))+
  geom_point(shape = 23, size = 3, fill = precip.palette[4])+
  geom_smooth(method = 'lm', aes(linetype = sig))+
  scale_linetype_manual(values = c('dashed', 'solid'))+
  stat_poly_eq(aes(label = paste( ..rr.label.., sep = "~~~")))+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = 'CS residual (kg/ha)')+
  theme_bw(base_size = 12)

forb.23 <- veg.precip.df%>%
  filter(water.yr == 2023)%>%
  filter(ecosite != 'Salt Flats')%>%
  group_by(ecosite)%>%
  mutate(p.val = summary(lm(forb.resid ~ lag2.gs))$coefficients[2,4],
         sig = ifelse(p.val > 0.1, 'n', 'y'))%>%
  ggplot( aes(x = lag2.gs, y = forb.resid))+
  geom_point(shape = 23, size = 3, fill= precip.palette[3] )+
  geom_smooth(method = 'lm', aes(linetype = sig))+
  scale_linetype_manual( values = c('dashed', 'solid'))+
  stat_poly_eq(aes(label = paste( ..rr.label.., sep = "~~~")), position = position_dodge())+
  facet_wrap(~ecosite)+
  labs(x = 'Growing season PPT (mm)', y = ' forb residual (kg/ha)')+
  theme_bw(base_size = 12)

total.mass.21 +total.mass.22 + total.mass.23 + wspg.21 +wspg.22+ wspg.23 + c3pg.21+  c3pg.22 + c3pg.23 + forb.21 +  forb.22  + forb.23 + plot_layout(ncol = 3, axis_titles = 'collect', guides = 'collect')


total.mass.21+ wspg.21  +total.mass.22 +wspg.22+ total.mass.23 + wspg.23 + plot_layout(ncol = 2, axis_titles = 'collect', guides = 'collect')


total.mass.21 +total.mass.22 + total.mass.23 + plot_layout(ncol = 1, axis_titles = 'collect', guides = 'collect')


wspg.21 +wspg.22+ wspg.23 + plot_layout(ncol = 1, axis_titles = 'collect', guides = 'collect')


##plot maps for 2021-23
maps2123 <- ppt.vars.sf%>%
  filter(water.yr %in% c(2021:2023))%>%
  ggplot()+
  geom_sf(aes(fill = gs.ppt))+
  scale_fill_viridis_c(option = 'turbo', direction = -1)+
  facet_wrap(~water.yr, ncol = 1)+
  labs(title = 'Growing Season Precipitation', fill = 'GS PPT')+
  #theme()
  theme_void(base_size = 13)+
  theme(text=element_text(family="sans"), plot.title = element_text(hjust = 0.5))

layout <- c('abg
            cdg
            efg' )

fig_case_study <- total.mass.21+ wspg.21  +total.mass.22 +wspg.22+ total.mass.23 + wspg.23 + maps2123 + plot_layout(design = layout, axis_titles = 'collect', guides = 'collect')


ggsave('case_study.png', fig_case_study, width = 10.5, height = 6)



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
    
    