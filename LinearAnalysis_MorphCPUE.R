library(tidyverse)#to make tidy data
library(lmerTest)#to get pvalues until I sort out which df are appropriate
library(performance)#to use check_model
library(see)#to use check_model
library(optimx)
library(emmeans)

##########################
###   CPUE   #############
##########################

##read in data
rg_cpue <- read.csv("CPUE_BaitedOnly.csv")%>%
  dplyr::select(-X)%>%
  mutate(Month = factor(Month, 
                        levels = c("May","June","July","August","September","October")),
         Site=factor(Site),
         SiteType=factor(Contamination))

#set constant for log10 transformation
t <- rg_cpue[rg_cpue$cpue_Total!=0,] 
c <- min(t$cpue_Total)/2

#center year at the first year goby were collected from contaminated sites
rg_cpue$cYear <- rg_cpue$Year - 2006

#run linear model
cpue_lmer <- lmer(log10(cpue_Total+c)~Year+I(Year^2)+SiteType+Month+(1|Site),data=rg_cpue)
check_model(cpue_lmer)
summary(cpue_lmer)

#cpue_trends <- emtrends(cpue_lmer, var="cYear")
#summary(cpue_trends)

######################
###   Morphology #####
######################

#read in the data
rg_morph <- read.csv("RGMorph_BaitedOnly.csv")%>%
  #  select(-X)%>%
  mutate(Sex = factor(Sex),Site = factor(Site), SiteType=factor(SiteType),
         Repro_Strategy = factor(Repro_Strategy), 
         Month = factor(Month, levels = c("May","June","July","August","September","October")))


#SL#
SL_lmer <- lmer(SL_cm ~ Year + SiteType + Sex + Month + (1 | Site), 
                data = rg_morph,
                control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
                #changed the default optimizer b/c of convergence warnings
                )
check_model(SL_lmer)
summary(SL_lmer)

#Body Mass#
Mass_lmer <- lmer(log10(Mass_g) ~ Year + SiteType + Sex + Month + (1 | Site), 
                data = rg_morph,
                control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
                )
check_model(Mass_lmer)
summary(Mass_lmer)

#Length-Weight relationship

ggplot(data = rg_morph,aes(x = SL_cm, y = Mass_g))+
  geom_point()#power of 3 gives most linear relationship
LW <- lm(Mass_g~I(SL_cm^3),data=rg_morph)
summary(LW)

#Body Condition#
Condition_lmer <- lmer(Fulton_Condition ~ Year + SiteType + Sex + Month + (1 | Site), 
                  data = rg_morph,
                  control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
)
check_model(Condition_lmer)
hist(resid(Condition_lmer)) #some extreme tails... b/c of a few rare observations that are very very high or very low
summary(Condition_lmer)

#GSI#
gsi.noNA  <- rg_morph[!is.na(rg_morph$GSI),]
gsi.no0 <- gsi.noNA[gsi.noNA$GSI!=0,]
gsi.c <- min(gsi.no0$GSI)/2 #making a constant for log10 transformation of (min value of GSI)/2

#exclude years prior to 2004 b/c GSI was rarely collected in the early years
rg_morph2004to2022 <- rg_morph[rg_morph$Year > 2003,]

rg_morph2004to2022$GSI_2 <- as.numeric(rg_morph2004to2022$Gonad_g)/(rg_morph2004to2022$Mass_g - as.numeric(rg_morph2004to2022$Gonad_g))

GSI_lmer <- lmer(log10(GSI+gsi.c) ~ Year + SiteType + Sex + Month + (1 | Site), 
                  data = rg_morph2004to2022,
                  control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
  )
check_model(GSI_lmer)
hist(resid(GSI_lmer))
summary(GSI_lmer)

GSI_lmer_2 <- lmer(log10(100*(GSI_2+gsi.c)) ~ cYear + SiteType + Sex + Month + (1 | Site), 
                          data = rg_morph2004to2022,
                          control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
  )

summary(GSI_lmer_2)

#HSI#
hsi.noNA  <- rg_morph[!is.na(rg_morph$HSI),]
hsi.no0 <- hsi.noNA[hsi.noNA$HSI!=0,]
hsi.c <- min(hsi.no0$HSI)/2 #making a constant for log10 transformation of (min value of HSI)/2

rg_HSI_no0_2004to2022 <- rg_morph2004to2022[rg_morph2004to2022$HSI > 0.000,]

HSI_lmer <- lmer(sqrt(HSI+hsi.c) ~ cYear + SiteType + Sex + Month + (1 | Site), 
                  data = rg_HSI_no0_2004to2022,
                  control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
        )
check_model(HSI_lmer)
hist(resid(HSI_lmer))
summary(HSI_lmer)

#Sigal also wanted a model where we include the 82 fish that were excluded b/c liver mass of 0
#and instead convert them to liver mass values of 0.001 (the smallest weight our scale
#could measure, so doing that below)

rg_morph2004to2022$Liver_test <- ifelse(rg_morph2004to2022$Liver_g < 0.00001, 0.001, rg_morph2004to2022$Liver_g)
rg_morph2004to2022$HSI_test <- rg_morph2004to2022$Liver_test/rg_morph2004to2022$Mass_g

HSI_lmer_test <- lmer(sqrt(HSI_test+hsi.c) ~ Year + SiteType + Sex + Month + (1 | Site), 
                 data = rg_morph2004to2022,
                 control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
  )

check_model(HSI_lmer_test)
hist(resid(HSI_lmer_test))
summary(HSI_lmer_test)


#not to use
rg_HSI_no0_2004to2022$HSI_2 <- rg_HSI_no0_2004to2022$Liver_g/(rg_HSI_no0_2004to2022$Mass_g-rg_HSI_no0_2004to2022$Liver_g)
HSI_lmer_2 <- lmer(sqrt(HSI_2+hsi.c) ~ Year + SiteType + Sex + Month + (1 | Site), 
                 data = rg_HSI_no0_2004to2022,
                 control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
  )

summary(HSI_lmer_2)

rg_morph2004to2022$HSI_2 <- rg_morph2004to2022$Liver_test/(rg_morph2004to2022$Mass_g-rg_morph2004to2022$Liver_test)
HSI_lmer_test_2 <- lmer(sqrt(HSI_2+hsi.c) ~ Year + SiteType + Sex + Month + (1 | Site), 
                   data = rg_morph2004to2022,
                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
)

summary(HSI_lmer_test_2)
