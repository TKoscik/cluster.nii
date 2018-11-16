cluster.valley <- function(data.3d, vol.3d,
                           mask, vol.mask,
                           min.size = 10,
                           tolerance = 0.95,
                           dir = c("native", "invert", "abs"),
                           save.dir,
                           file.prefix) {
  # # debug ----
  # rm(list=ls())
  # library(nifti.io)
  # data.3d <- "D:/data/duplex.20180219/m02/duplex.mask1.m02.coef.tvalue.nii"
  # vol.3d <- 2
  # mask <- "D:/data/duplex.20180219/m02/temp_mask.nii"
  # vol.mask <- 1
  # min.size = 10
  # dir = c("native", "invert", "abs")
  # source('D:/programs/nifti.cluster/R/cluster.3d.R')
  # source('D:/programs/nifti.cluster/R/get.neighbor.R')
  # #----

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
  neighbors <- get.neighbor(18L)
  grad.pairs <- matrix(c(3,17, 7,13, 9,11, 8,12, 6,14, 1,19, 15,5, 16,4, 18,2),
                       ncol=2, byrow = TRUE)
  valleys <- array(0,dim=dim(img))
  for (i in 1:nrow(which.xyz)) {
    for (j in 1:nrow(grad.pairs)) {
      new.xyz.1 <- which.xyz[i, ] + neighbors[grad.pairs[j,1], ]
      new.xyz.2 <- which.xyz[i, ] + neighbors[grad.pairs[j,2], ]
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
    for (j in 1:nrow(grad.pairs)) {
      new.xyz.1 <- which.hill[i, ] + neighbors[grad.pairs[j,1], ]
      new.xyz.2 <- which.hill[i, ] + neighbors[grad.pairs[j,2], ]
      if ((valleys[new.xyz.1[1], new.xyz.1[2], new.xyz.1[3]] == 1) &
          (valleys[new.xyz.2[1], new.xyz.2[2], new.xyz.2[3]] == 1)) {
        valleys[which.hill[i,1], which.hill[i,2], which.hill[i,3]] <- 1
      }
    }
  }
  hills <- mask - valleys
  # Cluster hills
  hills <- cluster.3d(in.array = hills, connectivity = 6L, min.size = 0)$clusters
  
  # assign valleys to adjacent hills
  neighbors <- get.neighbor(26L)
  which.valley <- data.frame(which(valleys == 1, arr.ind = T))
  which.valley$hill <- rep(0, nrow=nrow(which.valley))
  chk.valley <- TRUE
  while (chk.valley) {
    chk.valley <- FALSE
    unassigned <- which(which.valley$hill == 0)
    for (i in unassigned) {
      neighborhood <- matrix(c(which.valley[i,1] + neighbors[ ,1],
                       which.valley[i,2] + neighbors[ ,2],
                       which.valley[i,3] + neighbors[ ,3]), nrow=nrow(neighbors))
      nearby.hills <- hills[neighborhood]
      nearby.hills <- nearby.hills[nearby.hills != 0]
      if (length(nearby.hills) != 0) {
        unique.hills <- unique(nearby.hills)
        which.valley$hill[i] <- unique.hills[which.max(tabulate(match(nearby.hills, unique.hills)))]
        hills[which.valley[i,1], which.valley[i,2], which.valley[i,3]] <- which.valley$hill[i]
        chk.valley <- TRUE
      }
    }
  }
  
  # find remaining unlabelled valleys
  valleys <- ((valleys - (hills  > 0)*1) == 1)*1
  valleys <- cluster.3d(in.array = valleys, connectivity = 26L, min.size = min.size)$clusters
  
  # merge clusters
  valleys <- (valleys > 0)*max(hills) + valleys
  hills <- hills + valleys
  
  # sort clusters by size, descending order
  hill.table <- as.data.frame(table(hills))
  row.names(hill.table) <- NULL
  hill.table <- hill.table[-1, ]
  hill.table <- hill.table[order(hill.table$Freq, decreasing = TRUE), ]
  
  # assign small hills to nearby clusters
  tiny.hills <- as.numeric(hill.table$hills[hill.table$Freq < min.size])
  neighbors <- get.neighbor(26L)
  for (j in 1:length(tiny.hills)) {
    which.hill <- data.frame(which(hills == tiny.hills[j], arr.ind = T))
    which.hill$hill <- rep(0, nrow=nrow(which.valley))
    chk.hill <- TRUE
    while (chk.hill) {
      chk.hill <- FALSE
      unassigned <- which(which.hill$hill == 0)
      for (i in unassigned) {
        neighborhood <- matrix(c(which.hill[i,1] + neighbors[ ,1],
                                 which.hill[i,2] + neighbors[ ,2],
                                 which.hill[i,3] + neighbors[ ,3]), nrow=nrow(neighbors))
        nearby.hills <- hills[neighborhood]
        nearby.hills <- nearby.hills[nearby.hills != 0]
        if (length(nearby.hills) != 0) {
          unique.hills <- unique(nearby.hills)
          which.hill$hill[i] <- unique.hills[which.max(tabulate(match(nearby.hills, unique.hills)))]
          hills[which.hill[i,1], which.hill[i,2], which.hill[i,3]] <- which.hill$hill[i]
          chk.hill <- TRUE
        }
      }
    }
  }
  
  # order hills decreasing size
  hill.table <- as.data.frame(table(hills))
  row.names(hill.table) <- NULL
  hill.table <- hill.table[-1, ]
  hill.table <- hill.table[order(hill.table$Freq, decreasing = TRUE), ]
  hill.table <- hill.table[hill.table$Freq >= min.size, ]
  hills.ordered <- array(0, dim=dim(hills))
  for (i in 1:nrow(hill.table)) {
    hills.ordered[hills == hill.table$hills[i]] <- i
  }
  
# save output ------------------------------------------------------------------
    fname <- paste0(save.dir, "/", file.prefix, ".nii")
    init.nii(fname, dim(mask), pixdim, orient)
    write.nii.volume(fname, 1, hills.ordered)
}
