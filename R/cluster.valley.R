cluster.valley <- function(data.3d, vol.3d,
                           mask, vol.mask,
                           min.size = 10,
                           tolerance = 0.95,
                           dir = c("native", "invert", "abs"),
                           save.dir,
                           file.prefix) {
  # debug ----
  rm(list=ls())
  library(nifti.io)
  data.3d <- "D:/data/duplex.20180219/m02/duplex.mask1.m02.coef.tvalue.nii"
  vol.3d <- 2
  mask <- "D:/data/duplex.20180219/m02/temp_mask.nii"
  vol.mask <- 1
  min.size = 10
  dir = c("native", "invert", "abs")
  source('D:/programs/nifti.cluster/R/cluster.3d.R')
  #----

  pixdim <- unlist(nii.hdr(data.3d, "pixdim"))
  orient <- nii.orient(data.3d)

# Load data --------------------------------------------------------------------
  img <- read.nii.volume(data.3d, vol.3d)
  mask <- read.nii.volume(mask, vol.mask)
  img.dims <- dim(img)

# remove NAs -------------------------------------------------------------------
  img[is.na(img)] <- 0
  mask[is.na(mask)] <- 0
  mask <- (mask > 0) * 1

# apply mask -------------------------------------------------------------------
  img <- img * mask

# invert mask if desired -------------------------------------------------------
  img <- switch(dir[1],
                `native` = img,
                `invert` = img * -1,
                `abs` = abs(img),
                otherwise=stop("unknown direction value"))

# Initial clustering to find small clusters and remove islands -----------------
  temp.cl <- cluster.3d(in.array=(img>0)*1, connectivity=6L, min.size=min.size)$clusters
  temp.cl <- (temp.cl > 0)*1
  img <- img * temp.cl

# Find valleys based on local gradients ----------------------------------------
  which.xyz <- which(img > 0, arr.ind = TRUE)
  # grad.pt.1 <- matrix(c(-1, 0, 0,
  #                       0,-1, 0,
  #                       0, 0,-1), ncol=3, byrow = TRUE)
  # grad.pt.2 <- matrix(c( 1, 0, 0,
  #                        0, 1, 0,
  #                        0, 0, 1), ncol=3, byrow = TRUE)
  grad.pt.1 <- matrix(c(-1, 0, 0,
                        0,-1, 0,
                        0, 0,-1,
                        0,-1, 1,
                        0,-1,-1,
                        -1,-1, 0,
                        1,-1, 0,
                        1,-1, 0,
                        1, 0,-1,
                        1, 0, 1), ncol=3, byrow = TRUE)
  grad.pt.2 <- matrix(c( 1, 0, 0,
                         0, 1, 0,
                         0, 0, 1,
                         0, 1,-1,
                         0, 1, 1,
                         1, 1, 0,
                         -1, 1, 0,
                         -1, 1, 0,
                         -1, 0, 1,
                         -1, 0,-1), ncol=3, byrow = TRUE)
  # grad.pt.1 <- matrix(c(-1, 0, 0,
  #                        0,-1, 0,
  #                        0, 0,-1,
  #                        0,-1, 1,
  #                        0,-1,-1,
  #                       -1,-1, 0,
  #                        1,-1, 0,
  #                        1,-1, 0,
  #                        1, 0,-1,
  #                        1, 0, 1,
  #                       -1,-1, 1,
  #                       -1,-1,-1,
  #                       1,-1,-1,
  #                       1,-1, 1), ncol=3, byrow = TRUE)
  # grad.pt.2 <- matrix(c( 1, 0, 0,
  #                        0, 1, 0,
  #                        0, 0, 1,
  #                        0, 1,-1,
  #                        0, 1, 1,
  #                        1, 1, 0,
  #                       -1, 1, 0,
  #                       -1, 1, 0,
  #                       -1, 0, 1,
  #                       -1, 0,-1,
  #                       1, 1,-1,
  #                       1, 1, 1,
  #                       -1, 1, 1,
  #                       -1, 1,-1), ncol=3, byrow = TRUE)
  valleys <- array(0,dim=dim(img))
  for (i in 1:nrow(which.vxls)) {
    for (j in 1:nrow(grad.pt.1)) {
      new.xyz.1 <- which.xyz[i, ] + grad.pt.1[j, ]
      new.xyz.2 <- which.xyz[i, ] + grad.pt.2[j, ]
      pts <- c(NA,img[which.xyz[i,1], which.xyz[i,2], which.xyz[i,3]],NA)
      if (mask[new.xyz.1[1], new.xyz.1[2], new.xyz.1[3]] == 1) {
        pts[1] <- img[new.xyz.1[1], new.xyz.1[2], new.xyz.1[3]]
      }
      if (mask[new.xyz.2[1], new.xyz.2[2], new.xyz.2[3]] == 1) {
        pts[3] <- img[new.xyz.2[1], new.xyz.2[2], new.xyz.2[3]]
      }
      if (pts[2] == min(pts, na.rm = TRUE)) {
        valleys[which.xyz[i,1], which.xyz[i,2], which.xyz[i,3]] <- 1
      }
    }
  }
  
  hills <- mask - valleys
  # fill in gaps
  which.hill <- which(hills == 1, arr.ind=TRUE)
  for (i in 1:nrow(which.hill)) {
    for (j in 1:nrow(grad.pt.1)) {
      new.xyz.1 <- which.hill[i, ] + grad.pt.1[j, ]
      new.xyz.2 <- which.hill[i, ] + grad.pt.2[j, ]
      if ((valleys[new.xyz.1[1], new.xyz.1[2], new.xyz.1[3]] == 1) &
          (valleys[new.xyz.2[1], new.xyz.2[2], new.xyz.2[3]] == 1)) {
        valleys[which.hill[i,1], which.hill[i,2], which.hill[i,3]] <- 1
      }
    }
  }
  hills <- mask - valleys
  # Cluster hills
  hill.cluster <- cluster.3d(in.array = hills, connectivity = 26L, min.size = 0)
  
  # find valleys with no hill
  valley.cluster <- cluster.3d(in.array = valleys, connectivity = 26L, min.size=0)
  for (i in 1:max(valley.cluster$table$cluster.ordered)) {
    valley.chk <- sum((valley.cluster$clusters == i) * hill.cluster$clusters)
  } #*** this aint working, of course there are no hills in valleys
  
  
  
    # debug
    init.nii("D:/data/duplex.20180219/m02/valleys.nii", dim(mask), pixdim, orient)
    init.nii("D:/data/duplex.20180219/m02/hills.nii", dim(mask), pixdim, orient)
    write.nii.volume("D:/data/duplex.20180219/m02/valleys.nii", 1, valleys)
    write.nii.volume("D:/data/duplex.20180219/m02/hills.nii", 1, hill.cluster$clusters)
    
    #
    
    hill.ls <- which(img.hills == 1, arr.ind=TRUE)
    valley.ls <- which(valleys == 1, arr.ind=TRUE)


    connected <- array(0, dim = img.dims[1:3])
    num.clusters <- 0
    idx <- numeric()
    for (x in 1:img.dims[1]) {
      for (y in 1:img.dims[2]) {
        for (z in 1:img.dims[3]) {
          if (img.hills[x, y, z] == 1) {
            num.clusters <- num.clusters + 1
            current.pt <- t(as.matrix(c(x, y, z)))
            idx <- as.integer(colSums(t(current.pt[, 1:3, drop = FALSE] - 1) * cdim) + 1L)
            connected[idx] <- num.clusters
            while (length(idx) != 0) {
              img.hills[idx] <- 0
              neighbors = as.vector(apply(as.matrix(idx), 1, "+", offsets))
              neighbors = unique(neighbors[which(neighbors > 0)])
              idx = neighbors[which(img.hills[neighbors] != 0)]
              connected[idx] <- num.clusters
            }
          }
        }
      }
    }

    # Remove Valleys that form small clusters
    connect.valley <- array(0, dim = img.dims[1:3])
    num.clusters <- 0
    idx <- numeric()
    for (x in 1:img.dims[1]) {
      for (y in 1:img.dims[2]) {
        for (z in 1:img.dims[3]) {
          if (valleys[x, y, z] == 1) {
            num.clusters <- num.clusters + 1
            current.pt <- t(as.matrix(c(x, y, z)))
            idx <- as.integer(colSums(t(current.pt[, 1:3, drop = FALSE] - 1) * cdim) + 1L)
            connect.valley[idx] <- num.clusters
            while (length(idx) != 0) {
              valleys[idx] <- 0
              neighbors = as.vector(apply(as.matrix(idx), 1, "+", offsets))
              neighbors = unique(neighbors[which(neighbors > 0)])
              idx = neighbors[which(valleys[neighbors] != 0)]
              connect.valley[idx] <- num.clusters
            }
          }
        }
      }
    }
    valley.size <- as.data.frame(table(connect.valley))
    valley.size <- valley.size[-1, ]
    valley.size <- valley.size[order(valley.size$Freq, decreasing = TRUE), ]
    for (i in 1:nrow(valley.size)) {
      if (valley.size$Freq[i] < min.size) {
        connect.valley[connect.valley == valley.size$connect.valley[i]] <- 0
      }
    }
    valleys <- (connect.valley > 0) * 1
    valley.ls <- which(valleys == 1, arr.ind=TRUE)

    # add valleys to cluster of nearest labelled voxel
    for (i in 1:nrow(valley.ls)) {
      nearest <- sqrt(rowSums(cbind((valley.ls[i,1] - hill.ls[ ,1])^2,
                                    (valley.ls[i,2] - hill.ls[ ,2])^2,
                                    (valley.ls[i,3] - hill.ls[ ,3])^2)))
      nearest <- which(nearest == min(nearest))[1]
      connected[matrix(valley.ls[i, ], ncol=3)] <- connected[matrix(hill.ls[nearest, ], ncol=3)]
    }

    # Clean up small clusters
    connected.size <- as.data.frame(table(connected))
    connected.size <- connected.size[-1, ]
    connected.size <- connected.size[order(connected.size$Freq, decreasing = TRUE), ]
    while (any(connected.size$Freq < min.size)) {
      relabel.ls <- which(connected == connected.size$connected[nrow(connected.size)], arr.ind=TRUE)
      other.ls <- which(connected != connected.size$connected[nrow(connected.size)] & connected > 0, arr.ind=TRUE)
      for (i in 1:nrow(relabel.ls)) {
        nearest <- sqrt(rowSums(cbind((relabel.ls[i,1] - other.ls[ ,1])^2,
                                      (relabel.ls[i,2] - other.ls[ ,2])^2,
                                      (relabel.ls[i,3] - other.ls[ ,3])^2)))
        nearest <- which(nearest == min(nearest))[1]
        connected[matrix(relabel.ls[i, ], ncol=3)] <- connected[matrix(other.ls[nearest, ], ncol=3)]
      }
      connected.size <- as.data.frame(table(connected))
      connected.size <- connected.size[-1, ]
      connected.size <- connected.size[order(connected.size$Freq, decreasing = TRUE), ]
    }

    # Renumber voxels according to descending size
    out <- connected
    for (i in 1:nrow(connected.size)) {
      out[connected == connected.size$connected[i]] <- i
    }

    # save output
    fname <- paste0(save.dir, "/", file.prefix, ".pos.cluster.nii")
    init.nii(fname, img.dims, pixdim, orient)
    write.nii.volume(fname, 1, out)
  }

  
}
