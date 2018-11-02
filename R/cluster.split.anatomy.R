cluster.split.anatomy <- function(data.cluster, data.anatomy, which.cluster="all", cluster.sz) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  stopifnot(file.exists(data.cluster), file.exists(data.anatomy))
  img.dims <- nii.dims(data.cluster)[1:3]
  n.clusters <- nii.dims(data.cluster)[4]
  data.anatomy <- read.nii.volume(data.anatomy, 1)
  
  new.cluster <- vector("list", n.clusters)
  n.new <- numeric(n.clusters)
  for (i in 1:n.clusters) {
    cluster.temp <- read.nii.volume(data.cluster, i)
    cluster.temp <- data.anatomy * cluster.temp
    roi <- unique(as.vector(cluster.temp))
    roi <- roi[-which(roi==0)]
    if (length(roi) != 0) {
      new.cluster[[i]] <- array(0, dim=c(img.dims, length(roi)))
      sub.cluster.sz <- numeric(length(roi))
      for (j in 1:length(roi)) {
        sub.cluster.sz[j] <- sum(cluster.temp==roi[j])
        if (sub.cluster.sz[j] >= cluster.sz) {
          n.new[i] <- n.new[i] + 1
          new.cluster[[i]][ , , ,n.new[i]] <- cluster.temp == roi[j]
        }
      }
    } else {
      n.new[i] <- n.new[i] + 1
      new.cluster[[i]] <- array(read.nii.volume(data.cluster, i), dim=c(img.dims, 1))
    }
  }
  comb.cluster <- array(0, dim=c(img.dims, sum(n.new)))
  clust.count <- 0
  for (i in 1:n.clusters) {
    if (n.new[i] !=0) {
      for (j in 1:n.new[i]) {
        clust.count <- clust.count +1
        comb.cluster[ , , ,clust.count] <- new.cluster[[i]][ , , ,j]
      }
    }
  }
  
  init.nii(file.name = paste0(substr(data.cluster,1,(nchar(data.cluster)-3)), "anat.nii"),
           dims = dim(comb.cluster),
           pixdim = unlist(nii.hdr(data.cluster, "pixdim")),
           orient = nii.orient(data.cluster))
  for (i in 1:dim(comb.cluster)[4]) {
    write.nii.volume(nii.file = paste0(substr(data.cluster,1,(nchar(data.cluster)-3)), "anat.nii"),
                     vol.num = i,
                     values = comb.cluster[ , , ,i])
  }
}