cluster.coord <- function(cluster.nii,
                          cluster.vol = "all",
                          value.nii = NULL,
                          value.vol = NULL) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  #debug ----
  # rm(list=ls())
  # gc()
  # cluster.nii <- "C:/Users/tim.koscik/Documents/Data/Duplex.Decision/analyses/20170602/clusters/GoodANDBadEV.clusters.nii"
  # cluster.vol <- "all"
  # value.nii <- "C:/Users/tim.koscik/Documents/Data/Duplex.Decision/analyses/20170602/models/Decision-fMRI.EVAll.zvalue.nii"
  # value.vol <- 2
  # ----

  # Check specified volumes ----
  n.clusters <- 1:nii.dims(cluster.nii)[4]
  if (!is.numeric(cluster.vol)) {
    if (cluster.vol == "all") {
      cluster.vol <- n.clusters
    }
  }
  if (!all(cluster.vol %in% n.clusters)) {
    stop("Specified cluster volumes are invalid")
  }

  # Get transform ----
  tform <- nii.hdr(nii.file = cluster.nii,
                   field = c("srow_x", "srow_y", "srow_z"))
  tform <- rbind(tform$srow_x, tform$srow_y, tform$srow_z, c(0,0,0,1))


  # Setup output container ----
  coords <- list()
  coords$n.voxels <- numeric(length(cluster.vol))
  names(coords$n.voxels) <- paste0("Cluster", cluster.vol)

  coords$cog <- data.frame(
    world.x = numeric(length(cluster.vol)),
    world.y = numeric(length(cluster.vol)),
    world.z = numeric(length(cluster.vol)),
    mni.x = numeric(length(cluster.vol)),
    mni.y = numeric(length(cluster.vol)),
    mni.z = numeric(length(cluster.vol)))
  row.names(coords$cog) <- paste0("Cluster", cluster.vol)

  if (!is.null(value.nii)) {
    coords$peak <- data.frame(
      world.x = numeric(length(cluster.vol)),
      world.y = numeric(length(cluster.vol)),
      world.z = numeric(length(cluster.vol)),
      mni.x = numeric(length(cluster.vol)),
      mni.y = numeric(length(cluster.vol)),
      mni.z = numeric(length(cluster.vol)),
      value = numeric(length(cluster.vol)))
    row.names(coords$peak) <- paste0("Cluster", cluster.vol)
  }
  # Get World Coordinates ----
  for (i in 1:length(cluster.vol)) {

    cluster.mask <- read.nii.volume(cluster.nii, cluster.vol[i]) != 0
    cluster.coords <- which(cluster.mask, arr.ind = TRUE)
    cluster.values <- read.nii.volume(value.nii, value.vol)[cluster.mask]

    coords$n.voxels[i] <- sum(cluster.mask)

    coords$cog$world.x[i] <- round(weighted.mean(cluster.coords[ ,1], cluster.values))
    coords$cog$world.y[i] <- round(weighted.mean(cluster.coords[ ,2], cluster.values))
    coords$cog$world.z[i] <- round(weighted.mean(cluster.coords[ ,3], cluster.values))

    if (!is.null(value.nii)) {
      coords$peak$world.x[i] <- cluster.coords[which(abs(cluster.values) == max(abs(cluster.values)))[1],1]
      coords$peak$world.y[i] <- cluster.coords[which(abs(cluster.values) == max(abs(cluster.values)))[1],2]
      coords$peak$world.z[i] <- cluster.coords[which(abs(cluster.values) == max(abs(cluster.values)))[1],3]
      
      cluster.values <- read.nii.volume(value.nii, value.vol)
      coords$peak$value[i] <- cluster.values[coords$peak$world.x[i],
                                             coords$peak$world.y[i],
                                             coords$peak$world.z[i]]
    }
  }

  # Convert to MNI ----
  coords$cog[ ,4:6] <- t(tform %*% rbind(coords$cog$world.x-1,
                                         coords$cog$world.y-1,
                                         coords$cog$world.z-1,
                                         rep(1, length(cluster.vol))))[ ,-4]
  if (!is.null(value.nii)) {
    coords$peak[ ,4:6] <- t(tform %*% rbind(coords$peak$world.x-1,
                                            coords$peak$world.y-1,
                                            coords$peak$world.z-1,
                                            rep(1, length(cluster.vol))))[ ,-4]
  }

  # Output ----
  return(coords)
}
