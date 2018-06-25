cluster.p <- function(nii.p, vol.p=1,
                      nii.sign=NULL, vol.sign=vol.p,
                      nii.mask=NULL, vol.mask=vol.p,
                      p.thresh=0.001,
                      cluster.size, 
                      connectivity=26,
                      save.dir, file.name = NULL) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  pmap <- read.nii.volume(nii.p, vol.p)
  
  if (!is.null(nii.mask)) {
    mask <- read.nii.volume(nii.mask, vol.mask)
    mask <- (mask != 0) * 1
    mask[mask==0] <- NA
    pmap <- pmap * mask
  }
  
  if (!is.null(nii.sign)) {
    sign.temp <- sign(read.nii.volume(nii.sign, vol.sign))
    sign.map <- list()
    sign.map[[1]] <- (sign.temp==1)*1
    sign.map[[2]] <- (sign.temp==-1)*1
    n.dir <- 2
  } else {
    n.dir = 1
    sign.map <- list()
    sign.map[[1]] <- (pmap>0)*1
  }
  
  # Threshold
  pmap.dir <- vector("list", n.dir)
  for (i in 1:n.dir) {
    pmap.dir[[i]] <- (pmap < p.thresh & !is.na(pmap) & sign.map[[i]]==1) * 1
  }
  
  img.dims <- nii.dims(nii.p)
  
  # Initialize connectivity neighbourhood ----
  if (is.null(connectivity)) { connectivity <- 26 }
  n <- 3
  dimorder <- 1:3
  cdim <- cumprod(c(1, img.dims[dimorder][-n]))
  neighborhood <- switch(
    as.character(connectivity),
    `6` = t(matrix(c(1,2,2,2,1,2,2,2,1,3,2,2,2,3,2,2,2,3), nrow=3)),
    `18` = t(matrix(c(1,1,2,1,2,1,1,2,2,1,2,3,1,3,2,2,1,1,2,1,2,2,1,3,2,2,1,2,2,2,
                      2,2,3,2,3,1,2,3,2,2,3,3,3,1,2,3,2,1,3,2,2,3,2,3,3,3,2), nrow=3)),
    `26` = t(matrix(c(1,1,1,1,1,2,1,1,3,1,2,1,1,2,2,1,2,3,1,3,1,1,3,2,1,3,3,2,1,1,
                      2,1,2,2,1,3,2,2,1,2,2,3,2,3,1,2,3,2,2,3,3,3,1,1,3,1,2,3,1,3,
                      3,2,1,3,2,2,3,2,3,3,3,1,3,3,2,3,3,3), nrow=3)),
    stop("Unrecognized connectivity"))
  center.pt <- t(matrix(c(2,2,2), nrow=3))
  offsets <- as.integer(colSums(t(neighborhood[ ,1:3, drop=FALSE]-1)*cdim) + 1L) -
    as.integer(colSums(t(center.pt[ ,1:3, drop=FALSE]-1)*cdim) + 1L)
  
  for (j in 1:n.dir) {
    # Initialize cluster volume ----
    connected <- array(0, dim=img.dims[1:3])
    
    num.clusters <- 0
    idx <- numeric()
    for (x in 1:img.dims[1]) {
      for (y in 1:img.dims[2]) {
        for (z in 1:img.dims[3]) {
          if (pmap.dir[[j]][x,y,z] == 1) {
            num.clusters <- num.clusters + 1
            current.pt <- t(as.matrix(c(x,y,z)))
            idx <- as.integer(colSums(t(current.pt[ ,1:3, drop=FALSE]-1)*cdim) + 1L)
            connected[idx] <- num.clusters
            while (length(idx)!=0) {
              pmap.dir[[j]][idx] <- 0
              neighbors = as.vector(apply(as.matrix(idx), 1, '+', offsets))
              neighbors = unique(neighbors[which(neighbors > 0)])
              idx = neighbors[which(pmap.dir[[j]][neighbors]!=0)]
              connected[idx] <- num.clusters
            }
          }
        }
      }
    }
    
    connected.size <- as.data.frame(table(connected))
    connected.size <- connected.size[-1, ]
    connected.size <- connected.size[order(-connected.size$Freq), ]
    connected.size <- connected.size[which(connected.size$Freq >= cluster.size), ]
    n.clusters <- nrow(connected.size)
    if (n.clusters != 0) {
      
      cluster.array <- array(0, dim=c(img.dims[1:3], n.clusters))
      
      fname <- unlist(strsplit(nii.p, "[/]"))
      fname <- fname[(length(fname))]
      fname <- unlist(strsplit(fname, "[.]"))
      fname <- paste(fname[-length(fname)], collapse=".")
      if (is.null(file.name)) {
        fname <- paste0(save.dir, "/", fname,
                        ".vol", vol.p,
                        ".cl", connectivity,
                        ".p", p.thresh,
                        ".sz", cluster.size)
      } else {
        fname <- paste0(save.dir, "/",
                        file.name)
      }
      if (j == 1) {
        fname <- paste0(fname, ".pos.nii")
      } else {
        fname <- paste0(fname, ".neg.nii")
      }
      
      if (!dir.exists(save.dir)) { dir.create(save.dir) }
      
      init.nii(file.name=fname, dims=c(img.dims[1:3], n.clusters),
               pixdim=NULL, orient=NULL)
      
      for (k in 1:n.clusters) {
        write.nii.volume(
          nii.file=fname,
          vol.num=k,
          values=array(as.numeric(connected==connected.size$connected[k]), dim=img.dims[1:3]))
      }
    } else {
      warning("No clusters matching criteria detected")
    }
  }
}