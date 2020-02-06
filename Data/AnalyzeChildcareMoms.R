library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)
setwd("~/GitHub/welfare-meta-study/Data/")
D <- read.csv("ChildCareMoms.csv") %>% mutate(Subs = Subs/Emp, OOP = CPI*OOP/(Emp/100*UsePaid/100)) #, UsePaid = UsePaid/Emp)
ggplot(D,aes(x=Subs,y=UsePaid)) + geom_point() + geom_smooth(method="lm")

D0 <- D %>% filter(Arm==0) %>% melt(c("Site","Year","Arm","AgeMin","AgeMax")) %>% select(-Arm)
D1 <- D %>% filter(Arm>0) %>% melt(c("Site","Year","Arm","AgeMin","AgeMax"))
TE <- merge(D0,D1,by=c("Site","Year","AgeMin","AgeMax","variable")) %>% mutate(TE = value.y-value.x)

# so next we can regress and test this regression

ggplot(D,aes(x=Subs,y=OOP)) + geom_point() + geom_smooth(method="lm",se=FALSE)
mod1 <- lm(OOP ~ Subs,D)
summary(mod1)
mod2 <- lm(log(OOP) ~ log(Subs),D)
summary(mod2)

ggplot(D,aes(x=log(Subs),y=log(OOP))) + geom_point() + geom_smooth(method="lm",se=FALSE)
D$Price = exp(predict(mod2))
D$Price2 = predict(mod1)

write.csv(D,"ChildCareMoms_estimated.csv")


mod3 <- lm(UsePaid ~ Price2,D)
mod4 <- lm(log(UsePaid) ~ log(Price2),D)
mod5 <- lm(log(UsePaid) ~ log(OOP),D)
summary(mod3)
summary(mod4)
summary(mod5)
#ok
ggplot(D,aes(x=log(OOP),y=log(UsePaid))) + geom_point() + geom_smooth(method="lm",se=FALSE)

mod3 <- lm(UsePaid ~ Subs,D)
mod4 <- lm(log(UsePaid) ~ log(Subs),D)
summary(mod3)
summary(mod4)
