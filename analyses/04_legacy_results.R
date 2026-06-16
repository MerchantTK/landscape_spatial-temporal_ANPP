#CPER Legacy results 
#updated 1.26.26
library(tidyverse)
library(emmeans)
library(effects)
library(lme4)
library(patchwork)
#read in data

#full dataframe
veg.precip.df <- read.csv('data/legacy_full_df.csv')

#define fig palette
precip.palette <- c(  '#F78154','#389DE5', '#5DDCAC' ,'#9368B7')


`%notin%` <- Negate(`%in%`) 



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
lattice::qqmath(resid(m.ws.gs.cs.gs1.gs2 ))
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

#check log transform
m.log.c3.gs.cs.f.f1 <- lmer(log(c3pg) ~ ecosite*gs.ppt + ecosite*cs.ppt + ecosite*f.ppt + ecosite*lag.f+ (1|pairblock) + (1|water.yr),
                     data = veg.precip.df)
plot(m.log.c3.gs.cs.f.f1   )
summary(m.log.c3.gs.cs.f.f1)
MuMIn::r.squaredGLMM(m.log.c3.gs.cs.f.f1)
MuMIn::r.squaredGLMM(m.c3.gs.cs.f.f1)
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


sig.ws.lgs2 <- emmeans::emtrends( m.ws.gs.cs.gs1.gs2, ~ecosite, var = 'lag2.gs') %>%
  as.data.frame()%>%
  mutate(sig = ifelse(lower.CL < 0 & upper.CL > 0, 'n', 'y'))%>%
  dplyr::select(ecosite, sig)

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
  guides(linetype = 'none', fill = 'none', color = 'none')+
  scale_color_manual(values = precip.palette)+
  scale_fill_manual(values = precip.palette)+
  theme_bw(base_size = 11)+
  facet_wrap(~forcats::fct_relevel(var, c('lag.f',  'f.ppt', 'cs.ppt','gs.ppt')), nrow = 1, scales = 'free_x', labeller = as_labeller(var_names))+
  theme(axis.text.x= element_text(angle = 45, hjust = 1),
        legend.position = 'top')




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

