### Analysis to accompany "The role of early acoustic experience in song discrimination" Hudson, Creanza and Shizuka 2020

require(pacman)
p_load(dplyr,MASS,sjPlot,lme4)


# Load Data ---------------------------------------------------------------


all2017pbdata <- read.csv("~/RoleOfAcousticExperience_data.csv") #Read in RoleOfAcousticExperience_data.csv from your current working directory 


# Basic species recognition test ------------------------------------------

#Likelihood ratio test: does including species of playback improve model fit?

#Global model
lm1 <- lm(Track~pb+WCSPprop+Pre+as.numeric(Feather)+Sex+Clutch,data=all2017pbdata)
#Same model but with playback species (pb) removed
lmred <- lm(Track~WCSPprop+Pre+as.numeric(Feather)+Sex+Clutch,data=all2017pbdata) 
#The models are significantly different
anova(lm1,lmred)

#Other stats
median(all2017pbdata$Track[all2017pbdata$pb == "WCSP"])
range(all2017pbdata$Track[all2017pbdata$pb == "WCSP"])
median(all2017pbdata$Track[all2017pbdata$pb == "GCSP"])
range(all2017pbdata$Track[all2017pbdata$pb == "GCSP"])


# Figure 2: Chirp responses by species of playback song -------------------


xmins=c(0.75,1.75)
xmaxs=c(1.25,2.25)
ggplot(all2017pbdata, aes(x=pb,y=(Track/2)-Pre))  +stat_summary(fun.y=mean, colour="black", geom="errorbarh", xmin=xmins,xmax=xmaxs,lwd=1.25,show.legend = FALSE, aes(height=1)) + geom_point(position=position_jitter(width=.1,height=0.08),alpha=0.5,aes(x=pb,size=3,colour=pb), stroke=.4,show.legend = FALSE)+labs(x="Playback song",y="Standardized chirp response")+geom_hline(yintercept=0,lty=3)+theme(plot.margin=unit(c(5.5, 5.5, 8, 5.5), "points"))+annotate("text",1.5,68,label="*",size=10)+ylim(c(-10,75))



# Amplitude data ----------------------------------------------------------

# Read in data 
ampout <- read.csv(file ="~/AmplitudeOutputs.csv",header = T,stringsAsFactors = F) #Read in amplitude data "AmplitudeOutputs.csv"
ampout$logamp <- log(ampout$AvgAmp)

#Figure 3a
p <- ggplot(ampout,aes(x=SppFactor,y=logamp,col=SppFactor))+geom_violin(show.legend = FALSE)+scale_color_manual(values=c("#F8766D","#00BFC4"))+stat_summary(fun.y=median, geom="errorbarh", height=0, xmin=xmins,xmax=xmaxs,lwd=1.25,show.legend = FALSE)
p+labs(x="Species",y="Log amplitude of each song")

GCSPamps <- ampout[ampout$Species==" GCSP ",]
WCSPamps <- ampout[ampout$Species==" WCSP ",]
wilcox.test(GCSPamps$NewAmp,WCSPamps$NewAmp)

#Make separate data set with only WCSP playbacks
WCpbs1=subset(all2017pbdata,all2017pbdata$pb == "WCSP")

WCpbs <- cbind(aggregate(WCpbs1$Track~WCpbs1$Nest, FUN=mean),aggregate(WCpbs1$WCSPprop~WCpbs1$Nest, FUN=mean)[2],aggregate(WCpbs1$WCrate~WCpbs1$Nest, FUN=mean)[2],aggregate(WCpbs1$WCAmp~WCpbs1$Nest, FUN=mean, na.action=na.pass)[2], aggregate(WCpbs1$PropAmp~WCpbs1$Nest, FUN=mean, na.action=na.pass)[2]) 

colnames(WCpbs) <-c("Nest","Resp","WCSPprop","WCrate","PropAmp","WCAmp")





# Average nest response as a function of WCSP exposure --------------

summary(lm(WCpbs$Resp~WCpbs$WCrate))

#Figure 3b
plot(jitter(WCpbs$WCSPprop),jitter(WCpbs$Resp), ylab = "Average nest response to white-crowned sparrow song",xlab="Proportion of white-crowned sparrow song songs heard at nest \n relative to total Zonotrichia song",las=1,cex=1.5,pch=21,bg="lightseagreen")


##Mixed effects model (includes nest ID) with quasipoisson distribution and penalized quasi likelihood

pqlmod <- glmmPQL(Track~pb+WCSPprop+pb*WCSPprop+Pre+as.numeric(Feather)+Sex+Clutch,data=all2017pbdata,na.action = "na.fail", random = ~1|Nest,family = quasipoisson)
summary(pqlmod)
VarCorr(pqlmod) 

tab_model(pqlmod) #Table 1


pqlmod_rate <- glmmPQL(Track~pb+WCrate+pb*WCrate+Pre+as.numeric(Feather)+Sex+Clutch,data=all2017pbdata,na.action = "na.fail", random = ~1|Nest,family = quasipoisson)
summary(pqlmod_rate)
VarCorr(pqlmod_rate) 

tab_model(pqlmod_rate) #Table s3



# Average nest response as a function of relative amplituce of WC song --------



newdat1<-as.data.frame(WCpbs) 

newdat<-aggregate(WCpbs,by=list(newdat1$Resp,newdat1$WCAmp),length)[,1:3]
newdat$logratio_r <-log(round(newdat$Group.2,digits =4)) #Round one very very small value to zero to make visualizing easier


newdat$col <- "lightseagreen"

newdat$logratio_r <- replace(newdat$logratio_r,newdat$logratio_r=="-Inf","-7")
newdat$col <- replace(newdat$col,newdat$logratio_r=="-7","cadetblue2")

#Figure 3c
symbols(as.numeric(newdat$logratio_r),newdat$Group.1,circles=newdat$Nest/10,inches=FALSE,col=newdat$col,xlab="Log transformed white-crowned sparrow song amplitude\n relative to golden-crowned sparrow song amplitude (corrected for background noise)",ylab="Average nest-wide response to white-crowned sparrow song", bg=newdat$col,las=1)
legend("topright",legend = c("nests with WCSP song","nests with zero WCSP song"),fill=  c("lightseagreen","cadetblue2"),inset=c(0.4,0.02), cex=.75)
