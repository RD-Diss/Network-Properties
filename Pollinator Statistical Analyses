##Loading required packages 
install.packages('bipartite')
install.packages("visreg")
install.packages("MuMIn")
library(visreg)
library('bipartite')
library(raster)
library(sp)
library(igraph)
library(car)
library(MuMIn)

library(broom)

##Testing Statistical Analysis


StatTestPol <- read.csv("All_Pollinator_Data.csv" , header = T,)

StatTestPol

## Transformation methods to ensure normality when needed.

AverageLinksLog <- log10(StatTestPol[,c('PolTotAV')])
AverageLinksLog

ConnectanceLog <-  log10(StatTestPol[,c('PolTotC')])
ConnectanceLog

Con <- StatTestPol[,c('PolTotC')]
SQcon <- sqrt(Con)

Modularity <-  StatTestPol[,c('PolTotMod')]

Random_Rob <-  StatTestPol[,c('PolTotRrob')]

High_LeastAbundant <-  StatTestPol[,c('PolTotrobHAB')]

Low_LeastAbundant <-  StatTestPol[,c('PolTotrobLAB')]

High_MostAbundantF <-  StatTestPol[,c('PolTotrobHBTL')]

Low_MostAbundantF <-  StatTestPol[,c('PolTotrobLBTL')]

Gen <- StatTestPol[,c("PolTotGen")]

Vul<- StatTestPol[,c("PolTotVul")]  

Nested <- StatTestPol[,c("PolTotNest")]

Spec <- StatTestPol[,c("PolTotS")]

LogCon <- log10(StatTestPol$PolTotC)
LogCon

Annual_Precipitation<- StatTestPol$APrec
Precipitation_Seasonality <- StatTestPol$PSeas
Temperature_Seasonality <- StatTestPol$TSeas
Annual_Temperature <- StatTestPol$ATemp
Lat <- StatTestPol$Lat
Long <- StatTestPol$Long

Overlap <- StatTestPol$PolTotOver

##Connectance Precipitation

lineartest1 <- lm(LogCon~Annual_Precipitation+Precipitation_Seasonality) #What effect does increasing precipitation etc have on connectance

shapiro.test(lineartest1$residuals)

summary(lineartest1)

avPlots(lineartest1, ylab= "Log10 Connectance", main = "Connectance v Precipitation Measures")


##Connectance Temperature

lineartest2 <- lm(Con~Temperature_Seasonality+Annual_Temperature) #What effect does increasing precipitation etc have on connectance
shapiro.test (lineartest2$residuals) 
summary(lineartest2)

##Connectance Geo Position

lineartest3 <- lm(LogCon~Long)
shapiro.test (lineartest3$residuals)
summary(lineartest3)

visreg(lineartest3, ylab = "Log10 Connectance", main = "Connectance v Longitude")


##Modularity Precipitation ****************

lineartest4 <- lm(Modularity~Annual_Precipitation+Precipitation_Seasonality)
shapiro.test (lineartest4$residuals)
summary(lineartest4)

avPlots(lineartest4, ylab= "Modularity", main = "Modularity v Precipitation Measures in Pollinator Networks")

##Modularity Temperature

lineartest5 <- lm(Modularity~Annual_Temperature+Temperature_Seasonality)
shapiro.test (lineartest5$residuals)
summary(lineartest5)

avPlots(lineartest5, ylab= "Modularity", main = "Modularity v Temperature Measures")

##Modularity Geo Position

lineartest6 <- lm(Modularity~Long)
shapiro.test (lineartest6$residuals)
summary(lineartest6)

visreg(lineartest3, ylab = "Modularity", main = "Modularity v Longitude")


## Random Robustness Connectance

lineartest10 <- lm(Random_Rob~Con)
shapiro.test (lineartest10$residuals)
summary(lineartest10)
visreg(lineartest10,ylab= "Log10 Random Robustness", xlab="Connectance", main= "Random Robustness v Connectance")

## Least Abundant first Robustness (Higher trophic) v Connectance &&&&&&&&&&&&&&&&

lineartest11 <- lm(High_LeastAbundant~Con)
shapiro.test (lineartest11$residuals)
summary(lineartest11)
visreg(lineartest11,ylab= "Least Abundant First Robustness", xlab="Connectance", main= "Least Abundant Species First Robustness, Higher Trophic levels v Connectance")

## Least Abundant first Robustness (Lower Trophic) v Connectance

lineartest12 <- lm(Low_LeastAbundant~Con)
shapiro.test (lineartest12$residuals)
summary(lineartest12)
visreg(lineartest12,ylab= "Least Abundant First Robustness", xlab="Connectance", main= "Least Abundant Species First Robustness, Lower Trophic levels v Connectance")

## Most abundant first robustness (Higher trophic) v Connectance *****************

lineartest13 <- lm(High_MostAbundantF~Con)
shapiro.test (lineartest13$residuals)
summary(lineartest13)
visreg(lineartest13,ylab= "Most Abundant First Robustness", xlab="Connectance", main= "Most Abundant Species First Robustness, Higher Trophic levels v Connectance")

## Most abundant first robustness (Lower trophic) v Connectance

lineartest14 <- lm(Low_MostAbundantF~Con)
shapiro.test (lineartest14$residuals)
summary(lineartest14)
visreg(lineartest14,ylab= "Most Abundant First Robustness", xlab="Connectance", main= "Most Abundant Species First Robustness, Lower Trophic levels v Connectance")

## Random Robustness Modularity

lineartest16 <- lm(Random_Rob~Modularity)
shapiro.test (lineartest16$residuals)
summary(lineartest16)
visreg(lineartest16)


## Least Abundant first Robustness (Higher trophic) v modularity  ********************

logHigh_LeastAbundant <- log10(High_LeastAbundant)
logHigh_LeastAbundant <- c(logHigh_LeastAbundant)

lineartest17 <- lm(logHigh_LeastAbundant~Modularity)

shapiro.test (lineartest17$residuals)

summary(lineartest17)

visreg(lineartest17, ylab= "Least Linked First Robustness", xlab="Modularity", main= "Least Linked Species First Robustness, Higher Trophic levels v Modularity")


## Least Abundant first Robustness (Lower trophic) v modularity ********************

lineartest18 <- lm(Low_LeastAbundant~Modularity)

shapiro.test (lineartest18$residuals)

summary(lineartest18)

visreg(lineartest18, ylab= "Least Linked First Robustness", xlab="Modularity", main= "Least Linked Species First Robustness, Lower Trophic levels v Modularity 
       in Pollinator Networks")

## Most Abundant first Robustness (Higher Trophic) v modularity

lineartest19 <- lm(High_MostAbundantF~Modularity)

shapiro.test (lineartest19$residuals)

summary(lineartest19)

visreg(lineartest19, ylab= "Most Abundant First Robustness", xlab="Modularity", main= "Most Abundant Species First Robustness, Higher Trophic levels v Modularity")

## Most Abundant first Robustness (Lower Trophic) v modularity

lineartest20 <- lm(Low_MostAbundantF~Modularity)

shapiro.test (lineartest20$residuals)

summary(lineartest20)

visreg(lineartest20, ylab= "Most Abundant First Robustness", xlab="Modularity", main= "Most Abundant Species First Robustness, Lower Trophic levels v Modularity")


## Nestedness v temperature &&&

lineartest21 <- lm(Nested~Temperature_Seasonality+Annual_Temperature)
shapiro.test (lineartest21$residuals)

summary (lineartest21)

avPlots(lineartest21, ylab= "Nestedness", main = "Nestedness v Temperature Measures")


## Nestedness v Precipitation &&&

lineartest22 <- lm(Nested~Precipitation_Seasonality+Annual_Precipitation)
shapiro.test (lineartest22$residuals)

summary (lineartest22)

avPlots(lineartest22, ylab= "Nestedness", main = "Nestedness v Precipitation Measures")

## Nestedness v latitude

lineartest23 <- lm(Nested~Lat)
shapiro.test (lineartest23$residuals)

summary (lineartest23)

visreg(lineartest23, ylab= "Nestedness", xlab="Latitude", main = "Nestedness v Precipitation Measures")


## Generality v Temperature

lineartest24 <- lm(Gen~Annual_Temperature+Temperature_Seasonality)
shapiro.test (lineartest24$residuals)

summary (lineartest24)

avPlots(lineartest24, ylab= "Genarality", main = "Generality v Temperature Measures") 

## Gen v Prec

gen10 <- log10(Gen)
lineartest25 <- lm(gen10~Precipitation_Seasonality+Annual_Precipitation)
shapiro.test (lineartest25$residuals)

summary (lineartest25)

avPlots(lineartest25, ylab= "Generality", main = "Log10 Generality v Precipitation Measures")

## Vul v Temp

lineartest26 <- lm(Vul~Annual_Temperature+Temperature_Seasonality)
shapiro.test (lineartest26$residuals)

summary (lineartest26)

avPlots(lineartest26, ylab= "Vulnerability", main = "Vulnerability v Temperature Measures") 

##Vulnerability v Precipitation


lineartest27 <- lm(Vul~Annual_Precipitation+Precipitation_Seasonality)
shapiro.test (lineartest27$residuals)

summary (lineartest27)

avPlots(lineartest27, ylab= "Vulnerability", main = "Vulnerability v Precipitation Measures") 


## Consumer Diet Overlap v Temp 


lineartest29 <- lm(Overlap~Annual_Temperature+Temperature_Seasonality)
shapiro.test (lineartest29$residuals)
summary (lineartest29)
avPlots(lineartest29, ylab= "Log10 Diet Overlap", main = "Consumer Diet Overlap v Temperature Measures") 

## Consumer Diet Overlap V Precipitation &&&

Overlap10 <- log10(Overlap)
lineartest30 <- lm(Overlap~Annual_Precipitation+Precipitation_Seasonality)
shapiro.test (lineartest30$residuals)
summary (lineartest30)

avPlots(lineartest30, ylab= "Diet Overlap", main = "Consumer Diet Overlap v Precipitation Measures") 

## Nestedness AIC

NestSelect <- lm(Nested~Temperature_Seasonality+Annual_Temperature+Precipitation_Seasonality+Annual_Precipitation, na.action = "na.fail")
dredge(NestSelect)

## Mod AIC

ModSelect <- lm(Modularity~Annual_Precipitation+Precipitation_Seasonality+Temperature_Seasonality+Annual_Temperature, na.action = "na.fail")
dredge(ModSelect)

ModSelectssss <-lm(Modularity~Annual_Precipitation)
summary (ModSelectssss)

## GEN AIC

GenSelect <- lm(Gen~Annual_Precipitation+Precipitation_Seasonality+Temperature_Seasonality+Annual_Temperature, na.action = "na.fail")
dredge(GenSelect)


## VUL AIC

VulSelect <- lm(Vul~Annual_Precipitation+Precipitation_Seasonality+Temperature_Seasonality+Annual_Temperature, na.action = "na.fail")
dredge(VulSelect)

## Overlap AIC

OvSelect <- lm(Overlap~Annual_Precipitation+Precipitation_Seasonality+Temperature_Seasonality+Annual_Temperature, na.action = "na.fail")
dredge(OvSelect)
