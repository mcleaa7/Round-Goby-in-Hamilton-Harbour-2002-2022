library(ggplot2)
library(tidyverse)
library(lmerTest)
library(glmmTMB)
library(performance)
library(see)
library(mgcv)

dd2 <- read.csv("ProportionGuarders_BaitedOnly.csv") %>%
 # select(-X)%>%
  mutate(Site=factor(Site),
         SiteType=factor(SiteType, levels=c("Less Contaminated","More Contaminated")),
         Month = factor(Month, levels = c("April","May","June","July","August","September","October","November"))
  )

#We want to exclude April and Nov b/c we do not consistently collect data in those years
#dd2$Month <- factor(dd2$Month, levels = c("April","May","June","July","August","September","October","November"))
newdata <- dd2[ which(dd2$Month!='April'),]
newdata <- newdata[ which(newdata$Month!="November"),]
newdata <- droplevels(newdata)

newdata %>%
  summarise(PMTotal = sum(PM),
         SMTotal = sum(SM),
         TotalTotal = sum(total),
         SMPMPercent = ((PMTotal+SMTotal)/TotalTotal)*100
  )

ggplot(newdata, aes(x=Year,y=proportion,colour=SiteType))+
  geom_point(aes(size=total),position="jitter")+
  geom_smooth()


#relevel factor to make high contamination the reference level
newdata <- within(newdata, SiteType <- relevel(SiteType, ref = "More Contaminated"))

#assume binomial error distribution, weights for total # of males

m2.glmer <- glmmTMB(proportion ~ Year*SiteType + (1|Site), 
        data = newdata, family = betabinomial , weights = total)
check_model(m2.glmer)
qqnorm(residuals(m2.glmer,"pearson"), pch = 1, frame = FALSE)
qqline(residuals(m2.glmer,"pearson"), col = "darkgreen", lwd = 2)

summary(m2.glmer)

ggplot(newdata,aes(x=Year,y=resid(m2.glmer),colour=SiteType))+geom_point()+geom_smooth()

