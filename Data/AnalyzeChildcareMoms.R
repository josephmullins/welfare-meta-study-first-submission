library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)
setwd("~/GitHub/welfare-meta-study/Data/")
D <- read.csv("ChildCareMoms.csv") %>% mutate(Subs = Subs/Emp, OOP = OOP/UsePaid*100, UsePaid = UsePaid/Emp)
ggplot(D,aes(x=Subs,y=UsePaid)) + geom_point() + geom_smooth(method="lm")

D0 <- D %>% filter(Arm==0) %>% melt(c("Site","Year","Arm","AgeMin","AgeMax")) %>% select(-Arm)
D1 <- D %>% filter(Arm>0) %>% melt(c("Site","Year","Arm","AgeMin","AgeMax"))
TE <- merge(D0,D1,by=c("Site","Year","AgeMin","AgeMax","variable")) %>% mutate(TE = value.y-value.x)

# so next we can regress and test this regression
