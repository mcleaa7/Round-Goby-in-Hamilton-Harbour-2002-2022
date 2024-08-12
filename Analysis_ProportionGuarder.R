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

newdata$cyear <- newdata$Year - 2006 #centre year on 2006, when we first started
#collecting fish from the contaminated sites
newdata$obs <- factor(seq(nrow(newdata)))

ggplot(newdata, aes(x=Year,y=proportion,colour=SiteType))+
  geom_point(aes(size=total),position="jitter")+
  geom_smooth()

#so here's the model Rachel ran (or, well, one of them with interactions)
#deleted - can add back if nec

#binomial 
#m1cent.glmer <- glmer(proportion ~ cyear*SiteType + (1|Site) +(1|obs), 
#                      data = newdata, family="binomial", weights=total)

#qqnorm(resid(m1cent.glmer), pch = 1, frame = FALSE)
#qqline(resid(m1cent.glmer), col = "darkgreen", lwd = 2)
#hist(resid(m1cent.glmer))
#plot(m1cent.glmer)
#check_model(m1cent.glmer)
#assume binomial error distribution, weights for total # of males

m2.glmer <- glmmTMB(proportion ~ Year*SiteType + (1|Site), 
        data = newdata, family = betabinomial , weights = total)
check_model(m2.glmer)
qqnorm(residuals(m2.glmer,"pearson"), pch = 1, frame = FALSE)
qqline(residuals(m2.glmer,"pearson"), col = "darkgreen", lwd = 2)

summary(m2.glmer)

ggplot(newdata,aes(x=Year,y=resid(m2.glmer),colour=SiteType))+geom_point()+geom_smooth()

#run GAM following modelGS from Pedersen et al. (2019)
#The diagnostics show many issues with this approach so I didn't pursue further. 
#Using the glmmTMB model from above for reporting in results
#ART.gam <- gam(proportion ~ s(Year, Site, k=12, bs="fs", m=2) + 
#                 s(Year, by=SiteType, k=12, m=2) + 
#                 s(as.numeric(Month), k=6) + 
#                 SiteType,
#               data=newdata, method="REML",family=betar(link="logit")) 

#hist(resid(ART.gam))
#check_model(ART.gam)
##so not working. Maybe skip GAMs for this one. 
#gam.check(ART.gam) 
#par(mfrow=c(2,2))
#plot(ART.gam)
#par(mfrow=c(1,1))
#plot(newdata$Month,resid(ART.gam),notch=TRUE)
#abline(h=1)

