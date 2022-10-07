### Statistical Analysis Testing

StatTestPol <- read.csv("All_Pollinator_Data.csv" , header = T,)

StatTestPol

LogCon <- log10(StatTestPol$PolTotC) ##Log10 to ensure normality
LogCon

Annual_Precipitation<- StatTestPol$APrec
Precipitation_Seasonality <- StatTestPol$PSeas
Temperature_Seasonality <- StatTestPol$TSeas
Annual_Temperature <- StatTestPol$ATemp
Modularity <-  StatTestPar[,c('ParTotMod')]

##############
lineartest1 <- lm(LogCon~Annual_Precipitation+Precipitation_Seasonality+Temperature_Seasonality+Annual_Temperature) 


shapiro.test(lineartest1$residuals)

summary(lineartest1)

##Modularity Precipitation ****************

lineartest4 <- lm(Modularity~Annual_Precipitation+Precipitation_Seasonality)
shapiro.test (lineartest4$residuals)
summary(lineartest4)

avPlots(lineartest4, ylab= "Modularity", main = "Modularity v Precipitation Measures")

##Modularity Precipitation ****************

lineartest4.5<- lm(Modularity~Annual_Precipitation)
shapiro.test (lineartest4.5$residuals)
summary(lineartest4.5)

visreg(lineartest4.5, ylab = "Modularity", main = "Modularity v PREC AN")
