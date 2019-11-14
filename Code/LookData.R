setwd("~/GitHub/welfare-meta-study/Code/")
library(dplyr)
library(ggplot2)

D <- read.csv("../Data/QuarterlyData.csv")
D$Date <- D$Year + D$Quarter*0.25
D$Wage <- D$Earnings/D$LFP*100/(12*30)

D <- D[1:180,]

g1 = ggplot(D,aes(x=Date,y=LFP,linetype=Treatment,color=Site)) + geom_line()
g2 = ggplot(D,aes(x=Date,y=Wage,linetype=Treatment,color=Site)) + geom_line()
g3 = ggplot(D,aes(x=Date,y=Welfare,linetype=Treatment,color=Site)) + geom_line()
g4 = ggplot(D,aes(x=Date,y=Receipt,linetype=Treatment,color=Site)) + geom_line()
g7 = ggplot(D,aes(x=Date,y=FoodStamps,linetype=Treatment,color=Site)) + geom_line()

d1 <- D[D$Treatment=="Control",]
d2 <- D[D$Treatment=="Treatment",]
d1$TE <- d2$Wage - d1$Wage
d1$TE_emp <- d2$LFP - d1$LFP

d1 <- d1 %>% group_by(Site,Treatment) %>% mutate(Time = Date - Year[1])
g5 = ggplot(d1,aes(x=Time,y=TE,color=Site)) + geom_line()
g6 = ggplot(d1,aes(x=Time,y=TE_emp,color=Site)) + geom_line()

ggsave("Plot_LFP.png",g1)
ggsave("Plot_Wage.png",g2)
ggsave("Plot_Welfare.png",g3)
ggsave("Plot_Receipt.png",g4)


D <- read.csv("../Data/ChildOutcomes.csv")

d0 <- D[D$Treatment=="C",]
d1 <- D[D$Treatment!="C",]
d0 <- d0 %>% melt(c("Site","Treatment","YearOfMeasurement","Age.at.Measurement"))
d0 <- rename_(d0,value0 = "value") %>% select(-"Treatment")
d1 <- d1 %>% melt(c("Site","Treatment","YearOfMeasurement","Age.at.Measurement"))
d2 <- merge(d0,d1)

d2$TE <- d2$value - d2$value0
cols = c("Site","YearOfMeasurement","Age.at.Measurement","Treatment")
d2 <- select(d2,-c("value0","value"))


