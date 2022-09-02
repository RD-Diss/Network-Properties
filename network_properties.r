##Loading required packages 
library(raster)
library(sp)

###Polinator Data

##Extracting wordlclim data

r <- getData("worldclim",var="bio",res=2.5)

r <- r[[c(1,4,12,15)]]
names(r) <- c("ATemp","TSeas","APrec","PSeas")

##Importing network coordinates

Polmeta <- read.csv("Polreferences.csv" , header = T,)

Pollat <- c(Polmeta[,8])

Pollon <- c(Polmeta[,9])

coordinates <- data.frame(x=Pollon,y=Pollat)

points <- SpatialPoints(coordinates, proj4string = r@crs)

values <- extract(r,points)

PolDat <- cbind.data.frame(coordinates(points),values)

PolDat

plot(r[[1]])
plot(points,add=T)

###Parasite data

##Extracting wordlclim data

r <- getData("worldclim",var="bio",res=2.5)

r <- r[[c(1,4,12,15)]]
names(r) <- c("ATemp","TSeas","APrec","PSeas")

##Importing network coordinates

HParmeta <- read.csv("Parreferences.csv", header = T)

HParlat <- c(HParmeta[,9])
  
HParlon <- c(HParmeta[,10])
  

coordinates <- data.frame(x=HParlon,y=HParlat)

points <- SpatialPoints(coordinates, proj4string = r@crs)

values <- extract(r,points)

HPardat <- cbind.data.frame(coordinates(points),values)

HPardat

plot(r[[1]])
plot(points,add=T)

### Pollinator Network Analysis


##Load Network

net <- as.matrix(read.csv("M_PL_004.csv", header = T, row.names = 1))

## Binary Change 
net[net != 0] <- 1
net

S_r <- dim(net)[1]
S_c <- dim(net)[2]
S <- S_r + S_c
L <- sum(net)
C <- L/(S_r * S_c)

S
L
C

dir()

## Loop Test

Polnetworks <- list.files("C:/Users/danie/Desktop/Master's/Diss/Geographical Variability of Networks/Data/Web Of Life/Polinators")
Polnetworks

for (i in 1:length(Polnetworks)) {
  assign(paste0("network", i),
         as.matrix(read.csv(paste0("C:/Users/danie/Desktop/Master's/Diss/Geographical Variability of Networks/Data/Web Of Life/Polinators"), header = T, row.names = 1, Polnetworks [1])))
  
}

