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

##Testing Statistical Analysis


StatTestPar <- read.csv("All_Parasite_Data.csv" , header = T,)

StatTestPar

#######

Con <- StatTestPar[,c('ParTotC')]
SQcon <- sqrt(Con)

Modularity <-  StatTestPar[,c('ParTotMod')]

Random_Rob <-  StatTestPar[,c('ParTotRrob')]

High_LeastAbundant <-  StatTestPar[,c('ParTotrobHAB')]

Low_LeastAbundant <-  StatTestPar[,c('ParTotrobLAB')]

High_MostAbundantF <-  StatTestPar[,c('ParTotrobHBTL')]

Low_MostAbundantF <-  StatTestPar[,c('ParTotrobLBTL')]

Gen <- StatTestPar[,c("ParTotGen")]

Vul<- StatTestPar[,c("ParTotVul")]  

Nested <- StatTestPar[,c("ParTotNest")]

Spec <- StatTestPar[,c("ParTotS")]

LogCon <- log10(StatTestPar$ParTotC)
LogCon

Overlap <- StatTestPar$ParTotOverlap

Annual_Precipitation<- StatTestPar$APrec
Precipitation_Seasonality <- StatTestPar$PSeas
Temperature_Seasonality <- StatTestPar$TSeas
Annual_Temperature <- StatTestPar$ATemp
Lat <- StatTestPar$Lat
Long <- StatTestPar$Long


##Connectance Precipitation

lineartest1 <- lm(Con~Annual_Precipitation+Precipitation_Seasonality) #What effect does increasing precipitation etc have on connectance

shapiro.test(lineartest1$residuals)

summary(lineartest1)

avPlots(lineartest1, ylab= "Connectance", main = "Connectance v Precipitation Measures")

##Connectance Temperature **********

lineartest2 <- lm(Con~Temperature_Seasonality+Annual_Temperature) #What effect does increasing precipitation etc have on connectance
shapiro.test (lineartest2$residuals) 
summary(lineartest2)

avPlots(lineartest2, ylab= "Connectance", main = "Connectance v Temperature Measures")

##Connectance Geo Position *******

lineartest3 <- lm(LogCon~Lat)
shapiro.test (lineartest3$residuals)
summary(lineartest3)

visreg(lineartest3, ylab = "Log10 Connectance", main = "Connectance v Latitude")


##Modularity Precipitation ****************

lineartest4 <- lm(Modularity~Annual_Precipitation+Precipitation_Seasonality)
shapiro.test (lineartest4$residuals)
summary(lineartest4)

avPlots(lineartest4, ylab= "Modularity", main = "Modularity v Precipitation Measures")

##Modularity Temperature **********

lineartest5 <- lm(Modularity~Annual_Temperature+Temperature_Seasonality)
shapiro.test (lineartest5$residuals)
summary(lineartest5)

avPlots(lineartest5, ylab= "Modularity", main = "Modularity v Temperature Measures")

##Modularity Geo Position *******

lineartest6 <- lm(Modularity~Lat)
shapiro.test (lineartest6$residuals)
summary(lineartest6)

visreg(lineartest3, ylab = "Modularity", main = "Modularity v Latitude")

## Robustness Lat


lineartest7 <- lm(Random_Rob~Lat)
shapiro.test (lineartest7$residuals)
summary(lineartest7)
visreg(lineartest7)


## Random Robustness Connectance

lineartest10 <- lm(Random_Rob~Con)
shapiro.test (lineartest10$residuals)
summary(lineartest10)
visreg(lineartest10,ylab= "Random Robustness", xlab="Connectance", main= "Random Robustness v Connectance")

## Least Abundant first Robustness (Higher trophic) v Connectance &&&&&&&&&&&&&&&&

lineartest11 <- lm(High_LeastAbundant~Con)
shapiro.test (lineartest11$residuals)
summary(lineartest11)
visreg(lineartest11,ylab= "Least Abundant First Robustness", xlab="Connectance", main= "Least Abundant Species First Robustness, Higher Trophic levels v Connectance")

## Least Abundant first Robustness (Lower Trophic) v Connectance **************

lineartest12 <- lm(Low_LeastAbundant~Con) 
shapiro.test (lineartest12$residuals)
summary(lineartest12)
visreg(lineartest12,ylab= "Least Abundant First Robustness", xlab="Connectance", main= "Least Abundant Species First Robustness, Lower Trophic levels v Connectance")

## Most abundant first robustness (Higher trophic) v Connectance *****************

lineartest13 <- lm(High_MostAbundantF~Con)
shapiro.test (lineartest13$residuals)
summary(lineartest13)
visreg(lineartest13,ylab= "Most Abundant First Robustness", xlab="Connectance", main= "Most Abundant Species First Robustness, Higher Trophic levels v Connectance")

## Most abundant first robustness (Lower trophic) v Connectance ****************

lineartest14 <- lm(Low_MostAbundantF~Con)
shapiro.test (lineartest14$residuals)
summary(lineartest14)
visreg(lineartest14,ylab= "Most Abundant First Robustness", xlab="Connectance", main= "Most Abundant Species First Robustness, Lower Trophic levels v Connectance")


## Random Robustness Modularity

lineartest16 <- lm(Random_Rob~Modularity)
shapiro.test (lineartest16$residuals)
summary(lineartest16)
visreg(lineartest16, ylab= "Random Robustness", main = "Random Robustness v Modularity")


## Least Abundant first Robustness (Higher trophic) v modularity ************


lineartest17 <- lm(High_LeastAbundant~Modularity)

shapiro.test (lineartest17$residuals)

summary(lineartest17)

visreg (lineartest17, ylab="Least Abundant First Robustness", main= "Least Abundant Species First Robustness, Higher Trophic Levels v Modularity")


## Least Abundant first Robustness (Lower trophic) v modularity ********************

lineartest18 <- lm(Low_LeastAbundant~Modularity)

shapiro.test (lineartest18$residuals)

summary(lineartest18)

visreg(lineartest18, ylab= "Least Abundant First Robustness", xlab="Modularity", main= "Least Abundant Species First Robustness, Lower Trophic levels v Modularity")

## Most Abundant first Robustness (Higher Trophic) v modularity *************

lineartest19 <- lm(High_MostAbundantF~Modularity)

shapiro.test (lineartest19$residuals)

summary(lineartest19)

visreg(lineartest19, ylab= "Most Abundant First Robustness", xlab="Modularity", main= "Most Abundant Species First Robustness, Higher Trophic levels v Modularity")

## Most Abundant first Robustness (Lower Trophic) v modularity **************

lineartest20 <- lm(Low_MostAbundantF~Modularity)

shapiro.test (lineartest20$residuals)

summary(lineartest20)

visreg(lineartest20, ylab= "Most Linked First Robustness", xlab="Modularity", main= "Most Linked Species First Robustness, Lower Trophic levels v Modularity
       in Parasite Networks")

## Nestedness v temperature *************

lineartest21 <- lm(Nested~Temperature_Seasonality+Annual_Temperature)
shapiro.test (lineartest21$residuals)

summary (lineartest21)

avPlots(lineartest21, ylab= "Nestedness", main = "Nestedness v Temperature Measures in Parasite Networks")
avPlots(lineartest21)

## Nestedness v Precipitation 

lineartest22 <- lm(Nested~Precipitation_Seasonality+Annual_Precipitation)
shapiro.test (lineartest22$residuals)

summary (lineartest22)

avPlots(lineartest22, ylab= "Nestedness", main = "Nestedness v Precipitation Measures")

## Nestedness v latitude **************

lineartest23 <- lm(Nested~Lat)
shapiro.test (lineartest23$residuals)

summary (lineartest23)

visreg(lineartest23, ylab= "Nestedness", xlab="Latitude", main = "Nestedness v Latitude")


## Generality v Temperature &&&&&&&&&&&

gen10 <- log10(Gen)
lineartest24 <- lm(gen10~Annual_Temperature+Temperature_Seasonality)
shapiro.test (lineartest24$residuals)

summary (lineartest24)

avPlots(lineartest24, ylab= "Genarality", main = "Generality v Temperature Measures") 

## Gen v Prec *************

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

