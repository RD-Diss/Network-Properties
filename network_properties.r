##Loading required packages 
install.packages('bipartite')
library('bipartite')
library(raster)
library(sp)
library(igraph)

###Pollinator Data

##Extracting wordlclim data

Polr <- getData("worldclim",var="bio",res=2.5)

Polr <- Polr[[c(1,4,12,15)]]
names(Polr) <- c("ATemp","TSeas","APrec","PSeas")

##Importing network coordinates

Polmeta <- read.csv("Polreferences.csv" , header = T,)

Pollat <- c(Polmeta[,8])

Pollon <- c(Polmeta[,9])

coordinates <- data.frame(Long=Pollon,Lat=Pollat)

points <- SpatialPoints(coordinates, proj4string = Polr@crs)

values <- extract(Polr,points)

PolDat <- cbind.data.frame(coordinates(points),values)

PolDat

plot(Polr[[1]])
plot(points,add=T)

###Parasite data

##Extracting wordlclim data

Parr <- getData("worldclim",var="bio",res=2.5)

Parr <- Parr[[c(1,4,12,15)]]
names(Parr) <- c("ATemp","TSeas","APrec","PSeas")

##Importing network coordinates

HParmeta <- read.csv("Parreferences.csv", header = T)
HParmeta
HParlat <- c(HParmeta[,9])

HParlon <- c(HParmeta[,10])


coordinates <- data.frame(Long=HParlon,Lat=HParlat)

points <- SpatialPoints(coordinates, proj4string = Parr@crs)

values <- extract(Parr,points)

HParDat <- cbind.data.frame(coordinates(points),values)

HParDat

plot(Parr[[1]])
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


PolNetAn <- data.frame(PolTotS,PolTotL,PolTotAV,PolTotC,PolTotMod,PolTotrob)
PolNetAn

PolNetAll <-data.frame(PolNetAn,PolDat)
PolNetAll

### Parasite Network Analysis

##Network Loop 

ParTotS <- c() #Parasite all species number
ParTotL <- c() #Parasite all links
ParTotC <- c() #Parasite all connectance
ParTotAV <- c() #Parasite all average number of links per species
ParTotMod <- c() #Parasite all modularity
ParTotrob <- c() #Parasite all robustness

Parnetworks <- list.files("C:\\Users\\danie\\Desktop\\Parasites")
Parnetworks

for (i in 1:length(Parnetworks)) {
  
  net <- as.matrix(read.csv(paste0("C:\\Users\\danie\\Desktop\\Parasites\\",Parnetworks[i]), header = T, row.names = 1))
  
  
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
  
  
  ParTotS <- append(ParTotS,S)
  
  ParTotL <- append(ParTotL,L)
  
  ParTotC <- append(ParTotC,C)
  
  ParTotAV <- append(ParTotAV,AV)
  
  ParMod <- computeModules(net,forceLPA = T)@likelihood #Calculating Modularity
  
  ParTotMod <- append(ParTotMod,ParMod)
  
  Parex <- second.extinct(net,participant = "both", details = false) 
  
  Parrob <- robustness(Parex) #Calculating robustness
  ParTotrob <- append(ParTotrob, Parrob)
  
  
}

ParTotS
ParTotL
ParTotAV
ParTotC

ParTotMod

ParTotrob


ParNetAn <- data.frame(ParTotS,ParTotL,ParTotAV,ParTotC,ParTotMod,ParTotrob)
ParNetAn

ParNetAll <-data.frame(ParNetAn,HParDat)
ParNetAll

#modularity

test <- computeModules(net,forceLPA = T)@likelihood

plotModuleWeb(test)
test

#robustness

robtest <- second.extinct(net,participant = "both", details = false)
robustness(robtest)

