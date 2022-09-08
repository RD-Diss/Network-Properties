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
PolTotRrob <- c() #Pollinator all Random robustness
PolTotrobHAB <- c() #Pollinator all Higher Trophic Abundance robustness
PolTotrobLAB <- c() #Pollinator all Lower Trophic Abundance robustness
PolTotrobHBTL <- c() #Pollinator all Higher Trophic Generalist First robustness
PolTotrobLBTL <- c() #Pollinator all Lower Trophic Generalist First robustness

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
  
  Polex <- second.extinct(net,participant = "both", details = false) #Calculating Random robustness
  
  PolRrob <- robustness(Polex) 
  PolTotRrob <- append(PolTotRrob, PolRrob)
  
  PolHABex <- second.extinct(net,participant = "higher", method = "abun", details = false) #Calculating least abundant first robustness
  PolLABex<- second.extinct(net,participant = "lower", method = "abun", details = false)
  
  
  PolrobHAB <- robustness(PolHABex)
  PolTotrobHAB <- append(PolTotrobHAB,PolrobHAB)
  
  PolrobLAB <- robustness(PolLABex)
  PolTotrobLAB <- append(PolTotrobLAB, PolrobLAB)
  
  PolHBTLex <- second.extinct(net,participant = "higher", method = "degree", details = false) #Calculating most abundant to least robustness
  PolLBTLex <- second.extinct(net,participant = "lower", method = "degree", details = false)
  
  
  PolrobHBTL <- robustness(PolHBTLex)
  PolTotrobHBTL <- append(PolTotrobHBTL, PolrobHBTL)
  
  PolrobLBTL <- robustness(PolLBTLex)
  PolTotrobLBTL <- append(PolTotrobLBTL, PolrobLBTL)
  
  
}

PolTotS
PolTotL
PolTotAV
PolTotC

PolTotMod

PolTotRrob
PolTotrobHAB  
PolTotrobLAB  
PolTotrobHBTL 
PolTotrobLBTL


PolNetAn <- data.frame(PolTotS,PolTotL,PolTotAV,PolTotC,PolTotMod,PolTotRrob,PolTotrobHAB,PolTotrobLAB,PolTotrobHBTL,PolTotrobLBTL)
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
ParTotRrob <- c() #Parasite all Random robustness
ParTotrobHAB <- c() #Parasite all Higher Trophic Abundance robustness
ParTotrobLAB <- c() #Parasite all Lower Trophic Abundance robustness
ParTotrobHBTL <- c() #Parasite all Higher Trophic Generalist First robustness
ParTotrobLBTL <- c() #Parasite all Lower Trophic Generalist First robustness

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
  
  Parex <- second.extinct(net,participant = "both", details = false) #Calculating Random robustness
  
  ParRrob <- robustness(Parex) 
  ParTotRrob <- append(ParTotRrob, ParRrob)
  
  ParHABex <- second.extinct(net,participant = "higher", method = "abun", details = false) #Calculating least abundant first robustness
  ParLABex<- second.extinct(net,participant = "lower", method = "abun", details = false)
  
  
  ParrobHAB <- robustness(ParHABex)
  ParTotrobHAB <- append(ParTotrobHAB,ParrobHAB)
  
  ParrobLAB <- robustness(ParLABex)
  ParTotrobLAB <- append(ParTotrobLAB, ParrobLAB)
  
  ParHBTLex <- second.extinct(net,participant = "higher", method = "degree", details = false) #Calculating most abundant to least robustness
  ParLBTLex <- second.extinct(net,participant = "lower", method = "degree", details = false)
  
  
  ParrobHBTL <- robustness(ParHBTLex)
  ParTotrobHBTL <- append(ParTotrobHBTL, ParrobHBTL)
  
  ParrobLBTL <- robustness(ParLBTLex)
  ParTotrobLBTL <- append(ParTotrobLBTL, ParrobLBTL)
  
  
}

ParTotS
ParTotL
ParTotAV
ParTotC

ParTotMod

ParTotRrob
ParTotrobHAB  
ParTotrobLAB  
ParTotrobHBTL 
ParTotrobLBTL


ParNetAn <- data.frame(ParTotS,ParTotL,ParTotAV,ParTotC,ParTotMod,ParTotRrob,ParTotrobHAB, ParTotrobLAB,ParTotrobHBTL,ParTotrobLBTL)
ParNetAn

ParNetAll <-data.frame(ParNetAn,HParDat)
ParNetAll

#modularity

test <- computeModules(net,forceLPA = T)@likelihood

plotModuleWeb(test)
test

#robustness

robtest <- second.extinct(net,participant = "both", details = false)


robtestHAB <- second.extinct(net,participant = "higher", method = "abun", details = false)
robtestLAB<- second.extinct(net,participant = "lower", method = "abun", details = false)

robtestHBTL <- second.extinct(net,participant = "higher", method = "degree", details = false)
robtestLBTL<- second.extinct(net,participant = "lower", method = "degree", details = false)

robustness(robtest)
robustness(robtestHAB)
robustness(robtestLAB)
robustness(robtestHBTL)
robustness(robtestLBTL)

## Exporting Data

write.csv(PolNetAll, "All_Pollinator_Data.csv")


write.csv(ParNetAll, "All_Parasite_Data.csv")
