cluster.create<- function(data.4d, threshold, thresh.dir="gt", cluster.size, save.dir, connectivity=NULL, which.volume="all") {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  # Check 4D Input ---------------------------------------------------------------
  stopifnot(!missing(data.4d))
  if (is.character(data.4d)) {
    if (length(data.4d)==1 & dir.exists(data.4d)) {
      temp.4d <- data.4d
      if (length(temp.4d) == 1) { # Get filelist from specified directory
        fls <- list.files(path = temp.4d, pattern="*.nii$")
        data.4d <- paste(temp.4d, fls, sep="/")
      } else {
        data.4d <- temp.4d
      }
      for (i in 1:length(data.4d)) { # check if 4d files exist
        stopifnot(file.exists(data.4d[i]))
      }
    } else {
      for (i in 1:length(data.4d)) { # check if 4d files exist
        stopifnot(file.exists(data.4d[i]))
      }
    }
  }
  n.4d <- length(data.4d)
  
  # Check threshold input
  stopifnot(thresh.dir %in% c("ge", "gt", "le", "lt", "eq", "gtlt", "gele", "gtle", "gelt"))
  if (thresh.dir %in% c("gtlt", "gele", "gtle", "gelt") & length(threshold)==1) {
    threshold <- c(threshold,-threshold)
    sprintf("Thresholding masks at values %0.3f, %0.3f", threshold[1], threshold[2])
  }
  if (thresh.dir %in% c("gt", "ge", "lt", "le", "eq") & length(threshold)!=1) {
    sprintf("Thresholding masks at values %0.3f", threshold[1])
  }
  
  # Loop through files -------------------------------------------------------------------
  for (i in 1:n.4d) {
    # Get necessary NII file info) ----
    img.dims <- nii.dims(data.4d[i])
    n.vols <- img.dims[4] # Total number of volumes in file
    
    if (is.character(which.volume)) {
      if (which.volume == "all") {
        which.volume <- 1:n.vols
      } else {
        stop("Cannot parse volumes")
      }
    }
    stopifnot(all(which.volume %in% (1:n.vols)))
    
    # Loop through volumes ---------------------------------------------------------------
    for (j in which.volume) {
      volume <- read.nii.volume(data.4d[i], vol.num=j) # Load volume
      
      # Initialize connectivity neighbourhood ----
      if (is.null(connectivity)) {connectivity <- 26}
      n <- 3
      dimorder <- 1:3
      cdim <- cumprod(c(1, img.dims[dimorder][-n]))
      neighborhood <- switch(
        as.character(connectivity),
        `6` = t(matrix(c(1,2,2,2,1,2,2,2,1,3,2,2,2,3,2,2,2,3), nrow=3)),
        `18` = t(matrix(
          c(1,1,2,1,2,1,1,2,2,1,2,3,1,3,2,2,1,1,2,1,2,2,1,3,2,2,1,2,2,2,
            2,2,3,2,3,1,2,3,2,2,3,3,3,1,2,3,2,1,3,2,2,3,2,3,3,3,2), nrow=3)),
        `26` = t(matrix(
          c(1,1,1,1,1,2,1,1,3,1,2,1,1,2,2,1,2,3,1,3,1,1,3,2,1,3,3,2,1,1,
            2,1,2,2,1,3,2,2,1,2,2,3,2,3,1,2,3,2,2,3,3,3,1,1,3,1,2,3,1,3,
            3,2,1,3,2,2,3,2,3,3,3,1,3,3,2,3,3,3), nrow=3)),
        stop("Unrecognized connectivity"))
      center.pt <- t(matrix(c(2,2,2), nrow=3))
      offsets <- as.integer(colSums(t(neighborhood[ ,1:3, drop=FALSE]-1)*cdim) + 1L) -
        as.integer(colSums(t(center.pt[ ,1:3, drop=FALSE]-1)*cdim) + 1L)
      
      # Intensity threshold ----
      volume <- switch(thresh.dir,
        `gt`=array((volume > threshold[1])*1, dim=img.dims[1:3]),
        `ge`=array((volume >= threshold[1])*1, dim=img.dims[1:3]),
        `lt`=array((volume < threshold[1])*1, dim=img.dims[1:3]),
        `le`=array((volume <= threshold[1])*1, dim=img.dims[1:3]),
        `eq`=array((volume == threshold[1])*1, dim=img.dims[1:3]),
        `gtlt`=array(((volume > threshold[1])*1) * ((volume < threshold[2])*1), dim=img.dims[1:3]),
        `gele`=array(((volume >= threshold[1])*1) * ((volume <= threshold[2])*1), dim=img.dims[1:3]),
        `gtle`=array(((volume > threshold[1])*1) * ((volume <= threshold[2])*1), dim=img.dims[1:3]),
        `gelt`=array(((volume >= threshold[1])*1) * ((volume < threshold[2])*1), dim=img.dims[1:3]),
        error("Cannot parse threshold direction."))
      
      # Initialize cluster volume ----
      connected <- array(0, dim=img.dims[1:3])
      
      num.clusters <- 0
      idx <- numeric()
      for (x in 1:img.dims[1]) {
        for (y in 1:img.dims[2]) {
          for (z in 1:img.dims[3]) {
            if (volume[x,y,z] == 1) {
              num.clusters <- num.clusters + 1
              current.pt <- t(as.matrix(c(x,y,z)))
              idx <- as.integer(colSums(t(current.pt[ ,1:3, drop=FALSE]-1)*cdim) + 1L)
              connected[idx] <- num.clusters
              while (length(idx)!=0) {
                volume[idx] <- 0
                neighbors = as.vector(apply(as.matrix(idx), 1, '+', offsets))
                neighbors = unique(neighbors[which(neighbors > 0)])
                idx = neighbors[which(volume[neighbors]!=0)]
                connected[idx] <- num.clusters
              }
            }
          }
        }
      }
      
      connected.size <- as.data.frame(table(connected))
      connected.size <- connected.size[-1, ]
      if (nrow(connected.size) !=0) {
      connected.size <- connected.size[order(-connected.size$Freq), ]
      connected.size <- connected.size[which(connected.size$Freq >= cluster.size), ]
      
      cluster.array <- array(0, dim=c(img.dims[1:3], nrow(connected.size)))
      for (k in 1:nrow(connected.size)) {
        cluster.array[ , , , k] <- array(as.numeric(connected==connected.size$connected[k]), dim=img.dims[1:3])
      }
      
      # Write 4D output ----------------------------------------------------------------------
      fname <- unlist(strsplit(data.4d[i], "[/]"))
      fname <- fname[(length(fname))]
      fname <- unlist(strsplit(fname, "[.]"))
      fname <- paste(fname[-length(fname)], collapse=".")
      if (thresh.dir %in% c("gtlt", "gele", "gtle", "gelt")) {
        fname <- paste0(save.dir, "/", fname,
                        ".vol", j,
                        ".cl", connectivity,
                        ".th.", substr(thresh.dir,1,2), threshold[1],
                        ".", substr(thresh.dir,3,4), threshold[2],
                        ".sz", cluster.size, ".nii")
      }
      if (thresh.dir %in% c("gt", "ge", "lt", "le", "eq")) {
        fname <- paste0(save.dir, "/", fname,
                        ".vol", j,
                        ".cl", connectivity,
                        ".th.", thresh.dir, threshold[1],
                        ".sz", cluster.size, ".nii")
      }
      
      if (!dir.exists(save.dir)) { dir.create(save.dir) }
      
      init.nii(file.name = fname,
               dims = dim(cluster.array),
               pixdim = unlist(nii.hdr(data.4d[i], "pixdim")),
               orient = nii.orient(data.4d[i]))
      for (k in 1:dim(cluster.array)[4]) {
        write.nii.volume(fname, k, cluster.array[ , , ,k])
      }
      } else {
        warning("No clusters matching criteria detected")
      }
    }
  }
}


