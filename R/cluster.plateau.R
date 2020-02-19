cluster.plateau <- function(nii.p, vol.p=1,
                            nii.sign=NULL, vol.sign=vol.p,
                            nii.mask=NULL, vol.mask=vol.p,
                            p.peak=0.001, p.plateau=0.05,
                            cluster.size, 
                            connectivity=26,
                            save.dir, file.name = NULL) {

which.min <- function(X) {
  x <- order(X)
  x <- which(x==1)
  return(x)
}

# find peak clusters
peak <- cluster.p(nii.p, vol.p, nii.sign, vol.p, nii.mask, vol.mask, p.peak, 1,
                  connectivity, save.nii=F)
# find plateau clusters
plateau <- cluster.p(nii.p, vol.p, nii.sign, vol.p, nii.mask, vol.mask, p.plateau,
                     cluster.size, connectivity, save.nii=F)

n.dir <- length(peak)
for (j in 1:n.dir) {
  peak.ls <- sort(unique(as.numeric(peak[[j]])))[-1]
  peak.n <- length(peak.ls)
  if (peak.n > 0) {
    peak.coords <- matrix(NA, nrow=peak.n, ncol=3)
    for (i in 1:peak.n) {
      temp.coords <- which(peak[[j]] == i, arr.ind=TRUE)
      peak.coords[i, ] <- round(colMeans(temp.coords))
    }
    
    plateau.ls <- sort(unique(as.numeric(plateau[[j]])))[-1]
    plateau.n <- length(plateau.ls)
    
    cluster <- peak[[j]] * 0
    for (i in 1:plateau.n) {
      plateau.coords <- which(plateau[[j]] == i, arr.ind=TRUE)
      temp <- find.matches(peak.coords, plateau.coords)
      which.peak <- which(temp$distance == 0)
      
      if (length(which.peak) == 1) {
        cluster[plateau.coords] <- which.peak
      } else if (length(which.peak) > 1) {
        peak.dist <- as.matrix(dist(rbind(peak.coords[which.peak, ], plateau.coords)))[1:length(which.peak), (length(which.peak)+1):(length(which.peak)+nrow(plateau.coords))]
        peak.nearest <- apply(peak.dist, 2, which.min)
        cluster[plateau.coords] <- which.peak[peak.nearest]
      }
    }
    
    fname <- unlist(strsplit(nii.p, "[/]"))
    fname <- fname[(length(fname))]
    fname <- unlist(strsplit(fname, "[.]"))
    fname <- paste(fname[-length(fname)], collapse=".")
    if (is.null(file.name)) {
      fname <- paste0(save.dir, "/", fname,
                      "_vol", vol.p,
                      "_cl", connectivity,
                      "_ppeak", p.peak,
                      "_pplateau", p.plateau,
                      "_sz", cluster.size)
    } else {
      fname <- paste0(save.dir, "/", file.name)
    }
    fname <- paste0(fname, ".", names(peak)[j], ".nii")
    
    if (!dir.exists(save.dir)) { dir.create(save.dir) }
    init.nii(file.name=fname, dims=nii.dims(nii.p)[1:3], 
             pixdim=unlist(nii.hdr(nii.p, "pixdim")),
             orient=nii.orient(nii.p))
    write.nii.volume(nii.file=fname, vol.num=1, values=cluster)
  }
}
}
