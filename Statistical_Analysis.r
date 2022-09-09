### Statistical Analysis Testing

StatTestPol <- read.csv("All_Pollinator_Data.csv" , header = T,)

StatTestPol

LogCon <- log10(StatTestPol$PolTotC) ##Log10 to ensure normality
LogCon

Annual_Precipitation<- StatTestPol$APrec
Precipitation_Seasonality <- StatTestPol$PSeas
Temperature_Seasonality <- StatTestPol$TSeas
Annual_Temperature <- StatTestPol$ATemp

lineartest1 <- lm(LogCon~Annual_Precipitation+Precipitation_Seasonality+Temperature_Seasonality+Annual_Temperature) 


shapiro.test(lineartest1$residuals)

summary(lineartest1)
