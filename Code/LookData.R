#setwd("~/GitHub/welfare-meta-study/Code/")
setwd("~/welfare-meta-study/Code/")
library(dplyr)
library(ggplot2)

colors = c("#221E1D","#63AA9C","#E9633B") # color scheme
cback = "#fff5EE"  ###"#f5f4ef"
th_ = theme_minimal() + theme(rect = element_rect(fill = cback),panel.background = element_rect(fill = cback, color=cback),legend.position="left",axis.title.y=element_text(vjust=1.2),legend.background = element_rect(fill=cback,color=cback),panel.border = element_blank())


D <- read.csv("../Data/QuarterlyDataGraphs.csv")
D$Date <- D$Year + D$Quarter*0.25
D$Wage <- D$Earnings/D$LFP*100/(12*30)

D <- D[1:156,]

g1_control = ggplot(D,aes(x=Date,y=LFP,linetype=Treatment,color=Site,alpha=Treatment)) + geom_line(size=1.5) + scale_alpha_manual(values=c(1,0.2,0.2))

g1_treat = ggplot(D,aes(x=Date,y=LFP,linetype=Treatment,color=Site,alpha=Treatment)) + geom_line(size=1.5) + scale_alpha_manual(values=c(1,0.3,1))

g1_full = ggplot(D,aes(x=Date,y=LFP,linetype=Treatment,color=Site,alpha=Treatment)) + geom_line(size=1.5) + scale_alpha_manual(values=c(1,1,1))

tikz(file = "~/Dropbox/Research Projects/WelfareMetaAnalysis/Figures/LFPTrendControl.tex",width=4.,height = 3)
print(g1_control+ th_)
dev.off()

tikz(file = "~/Dropbox/Research Projects/WelfareMetaAnalysis/Figures/LFPTrendTreat.tex",width=4.,height = 3)
print(g1_treat + th_)
dev.off()

g2_control = ggplot(D,aes(x=Date,y=Welfare,linetype=Treatment,color=Site,alpha=Treatment)) + geom_line(size=1.5) + scale_alpha_manual(values=c(1,0.2,0.2))

g2_treat = ggplot(D,aes(x=Date,y=Welfare,linetype=Treatment,color=Site,alpha=Treatment)) + geom_line(size=1.5) + scale_alpha_manual(values=c(1,0.3,1))

g2_full = ggplot(D,aes(x=Date,y=Welfare,linetype=Treatment,color=Site,alpha=Treatment)) + geom_line(size=1.5) + scale_alpha_manual(values=c(1,1,1))

tikz(file = "~/Dropbox/Research Projects/WelfareMetaAnalysis/Figures/WelfareTrendControl.tex",width=4.,height = 3)
print(g2_control+ th_)
dev.off()

tikz(file = "~/Dropbox/Research Projects/WelfareMetaAnalysis/Figures/WelfareTrendTreat.tex",width=4.,height = 3)
print(g2_treat + th_)
dev.off()

g3 = ggplot(D,aes(x=Date,y=Wage,linetype=Treatment,color=Site,alpha)) + geom_line(size=1.5) + scale_alpha_manual(values=c(1,0.3,1))

tikz(file = "~/Dropbox/Research Projects/WelfareMetaAnalysis/Figures/WageTrendTreat.tex",width=4.,height = 3)
print(g3 + th_ + ylab("Wage (\\$/hr)"))
dev.off()


# g2 = ggplot(D,aes(x=Date,y=Wage,linetype=Treatment,color=Site)) + geom_line()
# g3 = ggplot(D,aes(x=Date,y=Welfare,linetype=Treatment,color=Site)) + geom_line()
# g4 = ggplot(D,aes(x=Date,y=Receipt,linetype=Treatment,color=Site)) + geom_line()
# g7 = ggplot(D,aes(x=Date,y=FoodStamps,linetype=Treatment,color=Site)) + geom_line()
# 
# d1 <- D[D$Treatment=="Control",]
# d2 <- D[D$Treatment=="Treatment",]
# d1$TE <- d2$Wage - d1$Wage
# d1$TE_emp <- d2$LFP - d1$LFP
# 
# d1 <- d1 %>% group_by(Site,Treatment) %>% mutate(Time = Date - Year[1])
# g5 = ggplot(d1,aes(x=Time,y=TE,color=Site)) + geom_line()
# g6 = ggplot(d1,aes(x=Time,y=TE_emp,color=Site)) + geom_line()
# 
# ggsave("Plot_LFP.png",g1)
# ggsave("Plot_Wage.png",g2)
# ggsave("Plot_Welfare.png",g3)
# ggsave("Plot_Receipt.png",g4)


D <- read.csv("../Data/ChildTreatmentEffects.csv") %>%
  select(Site,N_control,Achievement,AchieveBelowAverage,PB,BPI,Suspend) %>%
  rename("Behavioral Problems" = BPI, "Positive Behaviors" = PB) %>%
  melt(c("Site","N_control","Achievement"))

g4 = ggplot(D,aes(x=Achievement,y=value,size=0.1*sqrt(N_control))) + geom_point(color=colors[3]) + facet_wrap( ~ variable) + geom_smooth(method="lm",se=FALSE,color=colors[2]) + th_ + theme(legend.position = "none") + ylab(NULL)

tikz(file = "~/Dropbox/Research Projects/WelfareMetaAnalysis/Figures/TEScatter.tex",width=3.5,height = 3)
print(g4)
dev.off()


#ggplot(D,aes(x=Achievement,y=value,size=sqrt(N_control),color=variable,fill=variable)) + geom_point() + geom_smooth(method="lm",se=FALSE)

D <- read.csv("~/welfare-meta-study/Code/ModelFit.csv")
D$Date <- D$Year + D$Quarter*0.25
D$Treatment <- ordered(D$Treatment,levels=c(0,1,2),labels=c("Control","Treatment","Incentives"))

g5 = ggplot(D,aes(x=Date,y=LFP,color=Treatment,linetype=Case)) + geom_line(size=1.2) + facet_wrap(. ~ Site,scales="free_x") + scale_color_manual(values=colors) + th_ + theme(legend.position = "right")

tikz(file = "~/Dropbox/Research Projects/WelfareMetaAnalysis/Figures/ModelFit.tex",width=5.,height = 3.)
print(g5)
dev.off()

g6 = ggplot(D,aes(x=Date,y=Participation,color=Treatment,linetype=Case)) + geom_line(size=1.2) + facet_wrap(. ~ Site,scales="free") + scale_color_manual(values=colors) + th_ + theme(legend.position = "right")

tikz(file = "~/Dropbox/Research Projects/WelfareMetaAnalysis/Figures/ModelFitWelf.tex",width=5.,height = 3.)
print(g6)
dev.off()



