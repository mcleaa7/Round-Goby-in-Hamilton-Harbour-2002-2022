library(tidyverse)#to make tidy data
library(lmerTest)#to get pvalues until I sort out which df are appropriate
library(performance)#to use check_model
## BMB: see needs to be *installed*, but not *loaded* (auto-loads if available)
library(see)#to use check_model
library(tseries)#to test for autocorrelation
library(glmmTMB) #use glmtmb to add autocorrelation values
library(nlme) #use nlme to add autocorrelation values - but ran into bugs
library(bbmle) #Aic comparisons
library(broom.mixed) #for coef plots
library(ggplot2)

#some code to help for later
get_coefs <- function(mlist) { #function to pull coeffs for comparison of glmmTMB & lmer models
  coefs <- mlist |>
    purrr::map_dfr(\(x) broom.mixed::tidy(x,
                                          effects = "fixed",
                                          conf.int = TRUE),
                   .id = "model") |>
    dplyr::mutate(across(term, \(x) forcats::fct_inorder(factor(x))))
}

fit_glmmTMB <- function(lmer_model, data) {#function to make lmer models glmmtmb models and update for possible temporal autocorrelation effects
  
  form <- formula(lmer_model)
  base_model <- glmmTMB(form,
                        data = data,
                        REML=TRUE)  ## REML=FALSE is the default ...
  ## make sure the coefficients are all (almost) the same
  stopifnot(  ## error if not true
    ## test if equivalent up to a factor of 1e-6
    all.equal(fixef(base_model)$cond, fixef(lmer_model),
              tolerance = 1e-6)
  )
  ## year-to-year autocorr
  ar1year_model <- update(base_model, . ~ . + ar1(0 + factor(Year) | Site))
  ## month-to-month autocorr within year/site combinations
  ar1month_model <- update(base_model, . ~ . + ar1(0 + Month | factor(Year):Site))
  ## both
#  ar1both_model <- update(base_model, . ~ . + ar1(0 + Month | factor(Year):Site) +
 #                           ar1(0 + factor(Year) | Site))
 # return(lme4::namedList(base_model, ar1year_model, ar1month_model, ar1both_model))
  return(lme4::namedList(base_model, ar1year_model, ar1month_model))
}


theme_set(theme_bw(base_size=24))#sets theme for ggplot2

################
################
##CPUE##########
################
################

##read in data
rg_cpue <- read.csv("CPUE_BaitedOnly.csv")%>%
  mutate(Month = factor(Month, 
                        levels = c("May","June","July","August","September","October")),
         Site=factor(Site),
         SiteType=factor(Contamination)) %>%
  dplyr::select(Year, Site, Month, SiteType, cpue_Total)

##set constant for log10 transformation
tot <- rg_cpue[rg_cpue$cpue_Total!=0,] 
con <- min(tot$cpue_Total)/2

#run linear model
cpue_lmer <- lmer(log10(cpue_Total+con)~Year+I(Year^2)+
                    SiteType+Month+(1|Site),data=rg_cpue, na.action = na.omit)
check_model(cpue_lmer)
summary(cpue_lmer)
#scaling issues

CPUE_glmmTMB <- glmmTMB(log10(cpue_Total+con) ~ Year + I(Year^2) + SiteType + 
                          Month + (1|Site),
                        data = rg_cpue, REML=TRUE, na.action=na.omit)
#convergence issues

cpue2 <- update(CPUE_glmmTMB, data = mutate(rg_cpue, Year = Year -
                                              min(Year)))#update model with re-scale of year (starting at 2002 seems to be issue)
cpue3 <- update(cpue2, . ~ . + ar1(0 + factor(Year) | Site))
cpue4 <- update(cpue2, . ~ . + ar1(0 + Month | factor(Year):Site))

AICtab(cpue2, cpue3, cpue4, weights=TRUE, mnames = c("basic", "AR1_year", "AR1_month"))
ggcoef_cpue <- ggplot(cpue_coefs, aes(model, estimate, colour = model)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  facet_wrap(~term, scale = "free") +
  geom_hline(yintercept = 0, lty = 2)
print(ggcoef_cpue)

################
################
##Morphology####
################
################

#read in the data
rg_morph <- read.csv("RGMorph_BaitedOnly.csv")%>%
  #  select(-X)%>%
  mutate(Sex = factor(Sex),Site = factor(Site), SiteType=factor(SiteType),
         Repro_Strategy = factor(Repro_Strategy), 
         Month = factor(Month, levels = c("May","June","July","August","September","October")))

######
#SL
######

#run previous model
SL_lmer <- lmer(SL_cm ~ Year + SiteType + Sex + Month + (1 | Site), 
                  data = rg_morph)
#check_model(SL_lmer) ## BMB: hmm, why so slow? More data ...

SL_glmmTMB_models <- fit_glmmTMB(SL_lmer, data = rg_morph)
SL_glmmTMB_month <- glmmTMB(SL_cm ~ Year + SiteType + Sex + Month + (1|Site) + ar1(0 + Month | factor(Year):Site),
                            data=rg_morph, REML = TRUE)

AICtab(mlist = SL_glmmTMB_models,weights=TRUE)
ggcoef_cpue %+% get_coefs(SL_glmmTMB_models)

###############
###L/W analysis added at request of Associate Editor for BINV
###############
  
#run previous model
mass_SL_lmer <- lmer(log10(Mass_g) ~ SL_cm * Sex + Year + Month + SiteType + (1|Site),
                       data = rg_morph)
## check_model(mass_SL_lmer)

Mass_SL_glmmTMB_models <- fit_glmmTMB(mass_SL_lmer, data = rg_morph)
Mass_SL_month <- glmmTMB(log10(Mass_g) ~ SL_cm * Sex + Year + SiteType  + Month + 
                           (1|Site) + ar1(0 + Month | factor(Year):Site),
        data=rg_morph, REML = TRUE)

AICtab(mlist = Mass_SL_glmmTMB_models,weights=TRUE)
ggcoef_cpue %+% get_coefs(Mass_SL_glmmTMB_models)


##############
###GSI####
############
#prep data
gsi.noNA  <- rg_morph[!is.na(rg_morph$GSI),]
gsi.no0 <- gsi.noNA[gsi.noNA$GSI!=0,]
gsi.c <- min(gsi.no0$GSI)/2 #making a constant for log10 transformation of (min value of GSI)/2

#exclude years prior to 2004 b/c GSI was rarely collected in the early years
rg_morph2004to2022 <- rg_morph[rg_morph$Year > 2003,]

#run previous model
GSI_lmer <- lmer(log10(GSI+gsi.c) ~ Year + SiteType + Sex + Month + (1 | Site), 
                 data = rg_morph2004to2022,
                 na.action=na.omit)
check_model(GSI_lmer)

GSI_glmmTMB_models <- fit_glmmTMB(GSI_lmer, data = rg_morph2004to2022)
GSI_month <- glmmTMB(log10(GSI+gsi.c) ~ Year + SiteType + Sex + Month + 
                       (1|Site) + ar1(0 + Month | factor(Year):Site),
        data=rg_morph2004to2022, REML = TRUE)

AICtab(mlist = GSI_glmmTMB_models, weights=TRUE)
ggcoef_cpue %+% get_coefs(GSI_glmmTMB_models)


##########
##HSI####
#########

#prepare the data
hsi.noNA  <- rg_morph[!is.na(rg_morph$HSI),]
hsi.no0 <- hsi.noNA[hsi.noNA$HSI!=0,]
hsi.c <- min(hsi.no0$HSI)/2 #making a constant for log10 transformation of (min value of HSI)/2

#Assigning value of 0.001 to values less than 0.00001 (aka, zero)
rg_morph2004to2022$Liver_test <- ifelse(rg_morph2004to2022$Liver_g < 0.00001, 0.001, rg_morph2004to2022$Liver_g)
rg_morph2004to2022$HSI_test <- rg_morph2004to2022$Liver_test/rg_morph2004to2022$Mass_g

HSI_lmer <- lmer(sqrt(HSI_test+hsi.c) ~ Year + SiteType + Sex + Month + (1 | Site), 
                 data = rg_morph2004to2022)
check_model(HSI_lmer)

HSI_glmmTMB_models <- fit_glmmTMB(HSI_lmer, data = rg_morph2004to2022)
HSI_month <- glmmTMB(sqrt(HSI_test+hsi.c) ~ Year + SiteType + Sex + Month + 
                       (1|Site) + ar1(0 + Month | factor(Year):Site),
        data=rg_morph2004to2022, REML = TRUE)

AICtab(mlist = HSI_glmmTMB_models,weights=TRUE)
ggcoef_cpue %+% get_coefs(HSI_glmmTMB_models)
