# Relative distance from sound source of dB (decibel) sound levels
# Expressed relative to dB maximum
# See https://www.omnicalculator.com/physics/distance-attenuation
rel_dist <- function(dB){
  10^((max(dB) - dB)/20)
}

# Mean of x weighted by w
wmean <- function(x, w){
  w[is.na(w)] <- 0
  sum(x * w) / sum(w)
}

# Loss function
# INPUT
#   prm: vector of parameters, source co-ordinates, first x then y, finally max discernible dB
#   coord: dataframe of recorder coordinates (columns x and y)
#   dB: dataframe of decibel sound levels at each recorder, one column per source, NA if not detected
# OUTPUT
#   Loss value (sum of squared differences between recorded and infered decibel levels)
lf <- function(prm, coord, dB){
  npnt <- nrow(coord)
  nsht <- ncol(dB)
  errs <- c(!all(c("x","y") %in% names(coord)),
            !nrow(dB) == npnt,
            !length(prm) == 2*nsht+1)
  msgs <- c("coord must have columns x and y",
            "coord and dB should have the same number of rows",
            "prm must hold 2*npnt+1 values")
  if(any(errs))
    stop(paste(msgs[errs], collapse = "; "))
  
  ij <- expand.grid(i=1:npnt, j=1:nsht)
  x <- prm[1:nsht]
  y <- prm[nsht+(1:nsht)]
  minDB <- tail(prm, 1)
  dist_xy <- matrix(sqrt((coord$x[ij$i] - x[ij$j])^2 + 
                           (coord$y[ij$i] - y[ij$j])^2),
                    nrow=nrow(coord))
  reldist_xy <- dist_xy/min(dist_xy)
  nas <- is.na(dB)
  minDistNA <- min(dist_xy[nas])
  i <- which(dist_xy == minDistNA)[1]
  reldB <- -20*log(reldist_xy)
  (db_infrd <- reldB - reldB[i] + minDB)
  sum((dB[!nas] - db_infrd[!nas])^2)
}

# Estimates the location of multiple sound sources on a grid of recorders
# INPUT
#   coord: dataframe of recorder coordinates (columns x and y)
#   dB: dataframe of decibel sound levels at each recorder, one column per source, NA if not detected
# OUTPUT
#   A list with elements:
#     loc: data frame of estimated x,y source locations
#     maxUnheard: the infered maximum unheard decibel level
#     loss: the minimised loss function value
locate_sound <- function(coord, dB){
  prm <- c(x=apply(dB, 2, function(w) wmean(coord$x, w)) * rnorm(ncol(dB), 1, sd=0.01),
           y=apply(dB, 2, function(w) wmean(coord$y, w)) * rnorm(ncol(dB), 1, sd=0.01),
           minDB=min(dB, na.rm=T))
  res <- optim(prm, lf, coord=coord, dB=dB, method="SANN")
  xy <- as.data.frame(matrix(head(res$par,-1), ncol=2))
  names(xy) <- c("x", "y")
  list(loc = xy,
       maxUnheard = tail(res$par, 1),
       loss = res$value)
}
