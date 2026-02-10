source("acoustic_triangulation.R")

# Small manual example
# Decibel levels of 4 shots
dB <- data.frame(sht1 = c(NA, NA, 85, 
                          NA, NA, NA, 
                          NA, NA, NA),
                 sht2 = c(73, 86, NA, 
                          NA, NA, NA, 
                          NA, NA, NA),
                 sht3 = c(NA, NA, NA, 
                          NA, NA, 63, 
                          NA, 67, 76),
                 sht4 = c(NA, NA, NA, 
                          68, 80, NA, 
                          72, 75, NA))
# Coordinates of nine recorders
coord <- data.frame(expand.grid(x = 1:3,
                                y = 1:3))
res <- locate_sound(coord, dB)
plot(coord, asp=1, xlim=c(0, 4), ylim=c(0,4))
for(i in 1:ncol(dB)) text(coord$x, coord$y, dB[,i], col=i)
points(res$loc, pch=16, col=1:ncol(dB))


# Scalable example working from known shot locations
grd <- 8 # max grid points per side
sp <- 200 # grid point spacing
coord <- data.frame(expand.grid(x = sp*(0:grd),
                                y = sp*(0:grd)))
nsht <- 100 # number of shots
shts <- data.frame(x = runif(nsht, 0, sp*grd), # random shot locations
                   y = runif(nsht, 0, sp*grd))
# Calculate distance matrix and infer decibel levels at each recorder above a threshold
ij <- expand.grid(i=1:(grd+1)^2, j=1:nsht)
dist <- matrix(sqrt((coord$x[ij$i] - shts$x[ij$j])^2 + 
                      (coord$y[ij$i] - shts$y[ij$j])^2), 
               nrow = (grd+1)^2)
dB <- 160 - 20*log(dist) + rnorm(length(dist), sd=0.5)
dB[dB<60] <- NA
# Estimate shot locations
res <- locate_sound(coord, dB)
# Visualise
plot(coord, asp=1, xlim=sp*c(-1, grd+1), ylim=sp*c(-1, grd+1), cex=0.3)
for(i in 1:ncol(dB)) text(coord$x, coord$y, round(dB[,i]), col=i)
points(shts, pch=16, col=1:nsht)
points(res$loc, col=1:nsht)
table(apply(dB, 2, function(x) sum(!is.na(x)))) # number of recorders per shot
derr <- sqrt((shts$x-res$loc$x)^2 + (shts$y-res$loc$y)^2) # distance errors
mean(derr)
range(derr)
hist(derr)
