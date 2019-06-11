library(sp)
library(rgdal)

geodist <- read.table('data/Geodist_gabon_methane_samples_lat_long.txt')


colnames(geodist) <- c('Site', 'X', 'Y')
head(geodist)
tail(geodist)
coordinates(geodist) <- c('X', 'Y')
proj4string(geodist) <- CRS("+proj=longlat +datum=WGS84")  ## for example
res <- spTransform(geodist, CRS("+proj=utm +zone=32M ellps=WGS84"))
geodist <- as.data.frame(res)

 plot(geodist$X, geodist$Y)
(geodist$X <- geodist$X - min(geodist$X))
(geodist$Y <- geodist$Y - min(geodist$Y))
 plot(geodist$X, geodist$Y)
geodist$Site
write.table(geodist, 'output/geodist.tsv')
