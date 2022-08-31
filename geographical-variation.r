
net <- as.matrix(read.csv("./data//M_PL_006.csv", header = T, row.names = 1))

## Binary Change Test
net[net != 0] <- 1
net

S_r <- dim(net)[1]
S_c <- dim(net)[2]
S <- S_r + S_c
L <- sum(net)
C <- L/(S_r * S_c)
