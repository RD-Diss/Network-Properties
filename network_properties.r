library(raster)
library(sp)

###Polinator Data

r <- getData("worldclim",var="bio",res=2.5)

r <- r[[c(1,4,12,15)]]
names(r) <- c("ATemp","TSeas","APrec","PSeas")

Pollat <- c(68.35,56.066667,51.574994,46.553731,36.075391,17.916667,8.933333,5.583333,-29.616667,-33,-37.75336478,-42,-41.077069,-43.099531,20.7,33.4,37.016667,28.268611,-37.885354,81.816667,56.066667,-0.809883455,52.762395,42.2973,40.133333)

Pollon <- c(18.5,10.216667,-2.589902,-66.071245,-79.001843,-76.191667,-67.416667,-61.716667,30.133333,-69.283333,-58.28485516,-73.583333,-71.525206,171.720224,57.733333,131.5,-6.55,-16.605556,-57.841955,-71.3,9.266667,-91.12445681,1.575532,3.235217,-88.166667)

coordinates <- data.frame(x=Pollon,y=Pollat)

points <- SpatialPoints(coordinates, proj4string = r@crs)

values <- extract(r,points)

PolDat <- cbind.data.frame(coordinates(points),values)

PolDat

plot(r[[1]])
plot(points,add=T)

###Parasite data

r <- getData("worldclim",var="bio",res=2.5)

r <- r[[c(1,4,12,15)]]
names(r) <- c("ATemp","TSeas","APrec","PSeas")

HParlat <- c(40.183529,48.719961,52.233333,38.7775,51.166667,65.694476,49,50.283333,66.4,42.316667,45,46.166667,48.483333,57.766667,56.016667,55.016667,65.127638,39.661742,72.659588,47.734879,54.562508,51.851049,51.420062,34.10231964,16.73348654)

HParlon <- c(44.516602,19.493866,21.016667,68.423056,71.416667,148.710938,89,127.533333,129.166667,69.595833,80,74.333333,135.066667,40.933333,93.066667,82.933333,60.644531,46.761932,96.152344,82.036028,75.591431,95.765076,52.362815,-114.6822193,-88.98623133)

coordinates <- data.frame(x=HParlon,y=HParlat)

points <- SpatialPoints(coordinates, proj4string = r@crs)

values <- extract(r,points)

HPardat <- cbind.data.frame(coordinates(points),values)

HPardat

plot(r[[1]])
plot(points,add=T)

### Pollinator Network Analysis

PL004 <- read.csv("004_Bin.csv", header = F)


PL004 

names(PL004) <- paste("Pol", 1:84, sep = "")
row.names(PL004) <- paste("Pla", 1:12, sep = "")  
PL004  <- as.matrix(PL004)

n_Pla4 <- dim(PL004)[1]
n_Pla4

n_Pol4 <- dim(PL004)[2]
n_Pol4

S <- n_Pla4 + n_Pol4
S

L <- sum(PL004)
L

AV <- L/S
AV

CON <- L/(n_Pla4 * n_Pol4)
CON

### 006

dir()

PL006 <- read.csv("006_Bin.csv", header = F)


PL006 

names(PL006) <- paste("Pol", 1:61, sep = "")
row.names(PL006) <- paste("Pla", 1:17, sep = "")  
PL006  <- as.matrix(PL006)
PL006
n_Pla <- dim(PL006)[1]
n_Pla

n_Pol <- dim(PL006)[2]
n_Pol

S <- n_Pla + n_Pol
S

L <- sum(PL006)
L

AV <- L/S
AV

CON <- L/(n_Pla * n_Pol)
CON

##Binary Network Creation

dir()

PL006 <- read.csv("M_PL_006.csv", header = T)
PL006 

names(PL006) <- paste("Pol", 1:62, sep = "")
row.names(PL006) <- paste("Pla", 1:17, sep = "")  
PL006  <- as.matrix(PL006)
PL006

PL006[PL006 >= 1] <- 1
PL006
