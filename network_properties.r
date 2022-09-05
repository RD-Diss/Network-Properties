##Loading required packages 
install.packages('bipartite')
library('bipartite')
library(raster)
library(sp)
library(igraph)

###Pollinator Data

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

##Network Loop 

PolTotS <- c() #Pollinator all species number
PolTotL <- c() #Pollinator all links
PolTotC <- c() #Pollinator all connectance
PolTotAV <- c() #Pollinator all average number of links per species
PolTotMod <- c() #Pollinator all modularity
PolTotrob <- c() #Pollinator all robustness

Polnetworks <- list.files("C:\\Users\\danie\\Desktop\\Polinators")
Polnetworks

for (i in 1:length(Polnetworks)) {
  
  net <- as.matrix(read.csv(paste0("C:\\Users\\danie\\Desktop\\Polinators\\",Polnetworks[i]), header = T, row.names = 1))
  
  
  ## Binary Change 
  net[net != 0] <- 1
  net
  
  #Network Properties
  S_r <- dim(net)[1]
  S_c <- dim(net)[2]
  S <- S_r + S_c
  L <- sum(net)
  AV <- L/S
  C <- L/(S_r * S_c)
  
  
  PolTotS <- append(PolTotS,S)
  
  PolTotL <- append(PolTotL,L)

  PolTotC <- append(PolTotC,C)
  
  PolTotAV <- append(PolTotAV,AV)
  
  PolMod <- computeModules(net,forceLPA = T)@likelihood #Calculating Modularity
  
  PolTotMod <- append(PolTotMod,PolMod)
  
  Polex <- second.extinct(net,participant = "both", details = false) 
  
  Polrob <- robustness(Polex) #Calculating robustness
  PolTotrob <- append(PolTotrob, Polrob)
  
  
  }

PolTotS
PolTotL
PolTotAV
PolTotC

PolTotMod

PolTotrob
