##Loading required packages 
install.packages('bipartite')
install.packages("visreg")
library(visreg)
library('bipartite')
library(raster)
library(sp)
library(igraph)
library(car)

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
PolTotDietOver <- c() #Pollinator all Diet overlap
PolTotGen <- c() #Pollinator all Generality-measure of the degree of specialization in the web from the consumer perspective.
PolTotVul <- c() #Pollinator all Vulnerability- represents the degree of specialization found in the network from the resource perspective.
PolTotNest <- c() #Pollinator nestedness

Polnetworks <- list.files("C:\\Users\\danie\\Desktop\\Polinators")
Polnetworks

for (i in 1:length(Polnetworks)) {
  
  net <- as.matrix(read.csv(paste0("C:\\Users\\danie\\Desktop\\Polinators\\",Polnetworks[i]), header = T, row.names = 1))
  
  
  ## Binary Change 
  net[net != 0] <- 1
  net
  
  #Network Properties
  S_r <- dim(net)[1] ## Rows = Plant Species
  S_c <- dim(net)[2] ## Columns = Pollinator Species
  S <- S_r + S_c ## Total species number
  L <- sum(net) ## Number of links
  AV <- L/S ## Links per species
  C <- L/(S_r * S_c) ## Connectance
  #PolDietOv <- L/()
  #PolGen <- L/S_c
  
  
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
  
  #TPolDietOver <- append(TPolDietOver, PolDietOv)
  
  PolGen <- L/S_c
  PolTotGen <- append(PolTotGen, PolGen)
  
  PolVul <- L/S_r
  PolTotVul <- append(PolTotVul, PolVul)
  
  PolNest <- nested(net, method="NODF", rescale=FALSE, normalised=TRUE)
  PolTotNest <- append(PolTotNest,PolNest)
  
  
}


########## Diet Overlap

PolTotOverlapC <- c() # Pollinator all diet overlap connectance
PolTotOverlap <- c() # Pollinator all pollinator (consumer) diet overlap 

for (i in 1:length(Polnetworks)) {
  
  net <- as.matrix(read.csv(paste0("C:\\Users\\danie\\Desktop\\Polinators\\",Polnetworks[i]), header = T, row.names = 1))
  
  
  ## Binary Change 
  net[net != 0] <- 1
  net
  
  ## Construct consumer diet overlap graph
  
  net <- t(net) ##Transpose
  nrows <- NROW(net) #Count number of rows and columns
  ncols <- NCOL(net)
  out <- matrix(nrow = nrows, ncol = nrows) #Create an empty matrix with the number of rows
  
  for (i in 1:nrows){ #For every pollinator
    for (j in 1:nrows){ #For every other pollinator
      logical = any((net [i ,] & net[j,]) & net[j,]) ## the "any" checks if both vectors (rows) have 1s in the same location
      
      if (logical){ #If logical is true, produce 1 if not produce a 0
        out[i, j] <- 1
        out[j, i] <- 1
      } else {
        out[i, j] <- 0
        out[j, i] <- 0
      }
    }
  }
 
   
  #Diet overlap connectance
  OS_r <- dim(out)[1] ## Rows = Pollinator Species
  OS_r ## Total species number
  OL <- sum(out) ## Number of links
  
  OC <- OL/OS_r ## Connectance
  PolTotOverlapC <- append(PolTotOverlapC,OC)
  
  PolOverlap <- OC/(OS_c*(OS_c-1)/2) ## Pollinator (consumer) diet overlap connectance divided by the possible links between them
  PolTotOverlap <- append(PolTotOverlap,PolOverlap)
  
}

out

PolTotOverlap

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

PolTotDietOver
PolTotGen
PolTotVul
PolTotNest

PolNetAn <- data.frame(PolTotS,PolTotL,PolTotAV,PolTotC,PolTotMod,PolTotRrob,PolTotrobHAB,PolTotrobLAB,PolTotrobHBTL,PolTotrobLBTL, PolTotGen, PolTotVul, PolTotNest, PolTotOverlap)
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
ParTotDietOver <- c() #Parasite all Diet overlap
ParTotGen <- c() #Parasite all Generality-measure of the degree of specialization in the web from the consumer perspective.
ParTotVul <- c() #Parasite all Vulnerability- represents the degree of specialization found in the network from the resource perspective.
ParTotNest <- c() #Parasite nestedness

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
  
  ParGen <- L/S_c
  ParTotGen <- append(ParTotGen, ParGen)
  
  ParVul <- L/S_r
  ParTotVul <- append(ParTotVul, ParVul)
  
  ParNest <- nested(net, method="NODF", rescale=FALSE, normalised=TRUE)
  ParTotNest <- append(ParTotNest,ParNest)
}


########## Diet Overlap

ParTotOverlapC <- c() # Parasite all diet overlap connectance
ParTotOverlap <- c() # All Parasite (consumer) diet overlap 

for (i in 1:length(Parnetworks)) {
  
  net <- as.matrix(read.csv(paste0("C:\\Users\\danie\\Desktop\\Parasites\\",Parnetworks[i]), header = T, row.names = 1))
  
  
  ## Binary Change 
  net[net != 0] <- 1
  net
  
  ## Construct consumer diet overlap graph
  
  net <- t(net) ##Transpose
  nrows <- NROW(net) #Count number of rows and columns
  ncols <- NCOL(net)
  out <- matrix(nrow = nrows, ncol = nrows) #Create an empty matrix with the number of rows
  
  for (i in 1:nrows){ #For every parasite
    for (j in 1:nrows){ #For every other parasite
      logical = any((net [i ,] & net[j,]) & net[j,]) ## the "any" checks if both vectors (rows) have 1s in the same location
      
      if (logical){ #If logical is true, produce 1 if not produce a 0
        out[i, j] <- 1
        out[j, i] <- 1
      } else {
        out[i, j] <- 0
        out[j, i] <- 0
      }
    }
  }
  
  
  #Diet overlap connectance
  OS_r <- dim(out)[1] ## Rows = Parasite Species
  OS_r ## Total species number
  OL <- sum(out) ## Number of links
  
  OC <- OL/OS_r ## Connectance
  ParTotOverlapC <- append(ParTotOverlapC,OC)
  
  ParOverlap <- OC/(OS_c*(OS_c-1)/2) ## Parasite (consumer) diet overlap connectance divided by the possible links between them
  ParTotOverlap <- append(ParTotOverlap,ParOverlap)
  
}

out


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

ParTotGen 
ParTotVul
ParTotNest
ParTotOverlap

ParNetAn <- data.frame(ParTotS,ParTotL,ParTotAV,ParTotC,ParTotMod,ParTotRrob,ParTotrobHAB, ParTotrobLAB,ParTotrobHBTL,ParTotrobLBTL,ParTotGen, ParTotVul, ParTotNest, ParTotOverlap)
ParNetAn

ParNetAll <-data.frame(ParNetAn,HParDat)
ParNetAll



## Exporting Data

write.csv(PolNetAll, "All_Pollinator_Data.csv")


write.csv(ParNetAll, "All_Parasite_Data.csv")
