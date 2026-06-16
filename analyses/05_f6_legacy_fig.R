#Generate legacy figure
#

source('analyses/')

library(tidyverse)
library(effects)
library(emmeans)

#warm season 
aic.ws.df 


gs.ws1 <- effects::predictorEffect('gs.ppt', m.ws.gs.cs.gs1.gs2  , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = gs.ppt)%>%
  mutate(var = 'gs.ppt')%>%
  merge(sig.ws.gs, by = 'ecosite')



cs.ws1 <- effects::predictorEffect('cs.ppt',m.ws.gs.cs.gs1.gs2  , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = cs.ppt)%>%
  mutate(var = 'cs.ppt')%>%
  merge(sig.ws.cs, by = 'ecosite')



lgs.ws1 <- effects::predictorEffect('lag.gs',m.ws.gs.cs.gs1.gs2  , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = lag.gs)%>%
  mutate(var = 'lag.gs')%>%
  merge(sig.ws.lgs, by = 'ecosite')



lgs2.ws1 <- effects::predictorEffect('lag2.gs', m.ws.gs.cs.gs1.gs2  , focal.levels = 10)%>%
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
  labs(x = 'Precipitation (mm)', y = bquote('Predicted ANPP (kg ha'^-1*')'), title = 'Warm Season', fill = 'Ecosite', color = 'Ecosite')+
  guides(linetype = 'none', fill = 'none', color = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c( 'lag2.gs','lag.gs', 'cs.ppt', 'gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))+
  theme(axis.text.x= element_text(angle = 45, hjust = 1),
        legend.position = 'top')

var_names


#plot best model effect plot
#m.ws.gs.cs.gs1.gs2 = best warm season model
require(sjPlot)
require(broom)


## cool season 
aic.c3.df 
#m.c3.gs.cs.f.gs1.f1


gs.c31 <- effects::predictorEffect('gs.ppt', m.c3.gs.cs.f.f1   , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = gs.ppt)%>%
  mutate(var = 'gs.ppt')%>%
  merge(sig.c3.gs, by = 'ecosite')


cs.c31 <- effects::predictorEffect('cs.ppt',m.c3.gs.cs.f.f1   , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = cs.ppt)%>%
  mutate(var = 'cs.ppt')%>%
  merge(sig.c3.cs, by = 'ecosite')


f.c31 <- effects::predictorEffect('f.ppt',m.c3.gs.cs.f.f1   , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = f.ppt)%>%
  mutate(var = 'f.ppt')%>%
  merge(sig.c3.f, by = 'ecosite')




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
  labs(x = 'Precipitation (mm)', y = bquote('Predicted ANPP (kg ha'^-1*')'), title = 'Cool Season', fill = 'Ecosite', color = 'Ecosite')+
  geom_line(size = 2)+
  guides(linetype = 'none', fill = 'none', color = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c('lag.f',  'f.ppt', 'cs.ppt','gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))+
  theme(axis.text.x= element_text(angle = 45, hjust = 1),
        legend.position = 'top')





## forb 
aic.fo.df 



gs.fo1 <- effects::predictorEffect('gs.ppt', m.fo.gs.cs.f.f1    , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = gs.ppt)%>%
  mutate(var = 'gs.ppt')%>%
  merge(sig.fo.gs, by = 'ecosite')



cs.fo1 <- effects::predictorEffect('cs.ppt',m.fo.gs.cs.f.f1    , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = cs.ppt)%>%
  mutate(var = 'cs.ppt')%>%
  merge(sig.fo.cs, by = 'ecosite')



f.fo1 <- effects::predictorEffect('f.ppt', m.fo.gs.cs.f.f1    , focal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = f.ppt)%>%
  mutate(var = 'f.ppt')%>%
  merge(sig.fo.f, by = 'ecosite')



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
  guides(linetype = 'none', fill = 'none', color = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c('lag.f',  'f.ppt', 'cs.ppt','gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))+
  theme(axis.text.x= element_text(angle = 45, hjust = 1),
        legend.position = 'top')




### total mass: gs2

gs.tm1 <- effects::predictorEffect('gs.ppt', m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = gs.ppt)%>%
  mutate(var = 'gs.ppt')%>%
  merge(sig.tm.gs, by = 'ecosite')



cs.tm1 <- effects::predictorEffect('cs.ppt',m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = cs.ppt)%>%
  mutate(var = 'cs.ppt')%>%
  merge(sig.tm.cs, by = 'ecosite')


f.tm1 <- effects::predictorEffect('f.ppt', m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = f.ppt)%>%
  mutate(var = 'f.ppt')%>%
  merge(sig.tm.f, by = 'ecosite')


lgs.tm1 <- effects::predictorEffect('lag.gs',m.tm.gs.cs.f.gs1.gs2  , tmcal.levels = 10)%>%
  as.data.frame()%>%
  rename(var.value = lag.gs)%>%
  mutate(var = 'lag.gs')%>%
  merge(sig.tm.lgs, by = 'ecosite')


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
  guides(linetype = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c('lag2.gs',  'lag.gs', 'f.ppt','cs.ppt',  'gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))+
  theme(axis.text.x= element_text(angle = 45, hjust = 1),
        legend.position = 'top')

layout <- c(
  'a
b
b
b
c
c
c
d
d
d
e
e
e')
###combine all plots for figure 4
fig.6 <-  guide_area() + tm.slp.plt + ws.slp.plt  + cs.slp.plt + fo.slp.plt + 
  plot_layout(design = layout, guides = 'collect', axis_titles = 'collect') +
  plot_annotation(tag_levels = 'A')

fig.6

ggsave( 'figures/figure_6.pdf',fig.6, width = 18, height = 22, units = 'cm')
