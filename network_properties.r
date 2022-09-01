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



##Load Network

net <- as.matrix(read.csv("M_PL_006.csv", header = T, row.names = 1))

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

## Loop Test

Polnetworks <- list.files("C:/Users/danie/Desktop/Master's/Diss/Geographical Variability of Networks/Data/Web Of Life/Polinators")
Polnetworks

for (i in 1:length(Polnetworks)) {
  assign(paste0("network", i),
         as.matrix(read.csv(paste0("C:/Users/danie/Desktop/Master's/Diss/Geographical Variability of Networks/Data/Web Of Life/Polinators"), header = T, row.names = 1, Polnetworks [1])))
  
}
