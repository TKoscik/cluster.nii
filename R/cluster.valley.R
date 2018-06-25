cluster.valley <- function(tmap,
                           tvol,
                           pmap,
                           pvol,
                           mask,
                           mvol,
                           alpha = 0.05,
                           min.size = 10,
                           save.dir,
                           file.prefix) {

  pixdim <- unlist(nii.hdr(tmap, "pixdim"))
  orient <- nii.orient(tmap)

   # Load data
  tmap <- read.nii.volume(tmap, tvol)
  pmap <- read.nii.volume(pmap, pvol)
  mask <- read.nii.volume(mask, mvol)
  img.dims <- dim(tmap)

  # remove NAs
  tmap[is.na(tmap)] <- 0
  pmap[is.na(pmap)] <- 1
  mask[is.na(mask)] <- 0
  mask <- (mask > 0) * 1

  # apply mask
  tmap <- tmap * mask
  pmap[mask == 0] <- 1

  # mask tmap with thresholded pmask
  pmask <- (pmap < alpha) * 1
  tmap <- tmap * pmask

  # create positive nad negative thresholded tmaps
  tpos <- tmap * (tmap > 0)
  tneg <- tmap * (tmap < 0) * -1

  # check if positive and negative values exist
  if (sum(tpos > 0) != 0) { do.pos <- TRUE } else { do.pos <- FALSE }
  if (sum(tneg > 0) != 0) { do.neg <- TRUE } else { do.neg <- FALSE }

  # initialize offsets for clustering
  n <- 3
  dimorder <- 1:3
  cdim <- cumprod(c(1, img.dims[dimorder][-n]))
  neighborhood <- t(matrix(c(1, 2, 2, 2, 1, 2, 2, 2, 1,
                             3, 2, 2, 2, 3, 2, 2, 2, 3), nrow = 3))
  center.pt <- t(matrix(c(2, 2, 2), nrow = 3))
  offsets <- as.integer(colSums(t(neighborhood[, 1:3, drop = FALSE] - 1) * cdim) + 1L) -
    as.integer(colSums(t(center.pt[, 1:3, drop = FALSE] - 1) * cdim) + 1L)
  neighborhood <- neighborhood - 2

  if (do.pos) {
    # Initial clustering to find small clusters and remove islands
    tpos.temp <- (tpos > 0)*1
    connected <- array(0, dim = img.dims[1:3])
    num.clusters <- 0
    idx <- numeric()
    for (x in 1:img.dims[1]) {
      for (y in 1:img.dims[2]) {
        for (z in 1:img.dims[3]) {
          if (tpos.temp[x, y, z] == 1) {
            num.clusters <- num.clusters + 1
            current.pt <- t(as.matrix(c(x, y, z)))
            idx <- as.integer(colSums(t(current.pt[, 1:3, drop = FALSE] - 1) * cdim) + 1L)
            connected[idx] <- num.clusters
            while (length(idx) != 0) {
              tpos.temp[idx] <- 0
              neighbors = as.vector(apply(as.matrix(idx), 1, "+", offsets))
              neighbors = unique(neighbors[which(neighbors > 0)])
              idx = neighbors[which(tpos.temp[neighbors] != 0)]
              connected[idx] <- num.clusters
            }
          }
        }
      }
    }
    init.cl.ls <- as.data.frame(table(connected))
    init.cl.ls <- init.cl.ls[-1, ]
    init.cl.ls <- init.cl.ls[order(init.cl.ls$Freq, decreasing = TRUE), ]
    # Remove small islands
    remove.ls <- as.numeric(init.cl.ls$connected[init.cl.ls$Freq < min.size])
    connected[connected %in% remove.ls] <- 0

    # reset map
    tpos <- tpos * (connected > 0)

    xdiff <- sign(tpos - tpos[c(2:img.dims[1], 1), , ])
    xdiff <- (xdiff - xdiff[c(2:img.dims[1], 1), , ])
    x1 <- (xdiff > 0) * 1
    x1 <- x1[c(img.dims[1], 1:(img.dims[1]-1)), , ]

    ydiff <- sign(tpos - tpos[ ,c(2:img.dims[2], 1), ])
    ydiff <- (ydiff - ydiff[ ,c(2:img.dims[2], 1), ])
    y1 <- (ydiff > 0) * 1
    y1 <- y1[ ,c(img.dims[2], 1:(img.dims[2]-1)), ]

    zdiff <- sign(tpos - tpos[ , ,c(2:img.dims[3], 1)])
    zdiff <- (zdiff - zdiff[ , ,c(2:img.dims[3], 1)])
    z1 <- (zdiff > 0) * 1
    z1 <- z1[ , ,c(img.dims[3], 1:(img.dims[3]-1))]

    valleys <- ((x1 + y1 + z1) > 0) * 1
    valley.ls <- which(valleys == 1, arr.ind=TRUE)
    tpos.hills <- (tpos > 0)*1
    tpos.hills[valleys == 1] <- 0
    hill.ls <- which(tpos.hills == 1, arr.ind=TRUE)

    # fix holes in valleys
    for (i in 1:nrow(hill.ls)) {
      neighbors <- valleys[rbind(hill.ls[i, ] - c(1, 0, 0),
                                 hill.ls[i, ] + c(1, 0, 0))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }

      neighbors <- valleys[rbind(hill.ls[i, ] - c(0, 1, 0),
                                 hill.ls[i, ] + c(0, 1, 0))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }

      neighbors <- valleys[rbind(hill.ls[i, ] - c(0, 0, 1),
                                 hill.ls[i, ] + c(0, 0, 1))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }

      neighbors <- valleys[rbind(hill.ls[i, ] - c(1, 1, 0),
                                 hill.ls[i, ] + c(1, 1, 0))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }

      neighbors <- valleys[rbind(hill.ls[i, ] - c(1, 0, 1),
                                 hill.ls[i, ] + c(1, 0, 1))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }

      neighbors <- valleys[rbind(hill.ls[i, ] - c(0, 1, 1),
                                 hill.ls[i, ] + c(0, 1, 1))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }
    }

    tpos.hills <- tpos > 0
    tpos.hills[valleys == 1] <- 0
    hill.ls <- which(tpos.hills == 1, arr.ind=TRUE)
    valley.ls <- which(valleys == 1, arr.ind=TRUE)


    connected <- array(0, dim = img.dims[1:3])
    num.clusters <- 0
    idx <- numeric()
    for (x in 1:img.dims[1]) {
      for (y in 1:img.dims[2]) {
        for (z in 1:img.dims[3]) {
          if (tpos.hills[x, y, z] == 1) {
            num.clusters <- num.clusters + 1
            current.pt <- t(as.matrix(c(x, y, z)))
            idx <- as.integer(colSums(t(current.pt[, 1:3, drop = FALSE] - 1) * cdim) + 1L)
            connected[idx] <- num.clusters
            while (length(idx) != 0) {
              tpos.hills[idx] <- 0
              neighbors = as.vector(apply(as.matrix(idx), 1, "+", offsets))
              neighbors = unique(neighbors[which(neighbors > 0)])
              idx = neighbors[which(tpos.hills[neighbors] != 0)]
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

  if (do.neg) {
    # Initial clustering to find small clusters and remove islands
    tneg.temp <- (tneg > 0)*1
    connected <- array(0, dim = img.dims[1:3])
    num.clusters <- 0
    idx <- numeric()
    for (x in 1:img.dims[1]) {
      for (y in 1:img.dims[2]) {
        for (z in 1:img.dims[3]) {
          if (tneg.temp[x, y, z] == 1) {
            num.clusters <- num.clusters + 1
            current.pt <- t(as.matrix(c(x, y, z)))
            idx <- as.integer(colSums(t(current.pt[, 1:3, drop = FALSE] - 1) * cdim) + 1L)
            connected[idx] <- num.clusters
            while (length(idx) != 0) {
              tneg.temp[idx] <- 0
              neighbors = as.vector(apply(as.matrix(idx), 1, "+", offsets))
              neighbors = unique(neighbors[which(neighbors > 0)])
              idx = neighbors[which(tneg.temp[neighbors] != 0)]
              connected[idx] <- num.clusters
            }
          }
        }
      }
    }
    init.cl.ls <- as.data.frame(table(connected))
    init.cl.ls <- init.cl.ls[-1, ]
    # Remove small islands
    remove.ls <- which(init.cl.ls$Freq < min.size)
    connected[connected %in% remove.ls] <- 0

    # reset map
    tneg <- tneg * (connected > 0)

    xdiff <- sign(tneg - tneg[c(2:img.dims[1], 1), , ])
    xdiff <- (xdiff - xdiff[c(2:img.dims[1], 1), , ])
    x1 <- (xdiff > 0) * 1
    x1 <- x1[c(img.dims[1], 1:(img.dims[1]-1)), , ]

    ydiff <- sign(tneg - tneg[ ,c(2:img.dims[2], 1), ])
    ydiff <- (ydiff - ydiff[ ,c(2:img.dims[2], 1), ])
    y1 <- (ydiff > 0) * 1
    y1 <- y1[ ,c(img.dims[2], 1:(img.dims[2]-1)), ]

    zdiff <- sign(tneg - tneg[ , ,c(2:img.dims[3], 1)])
    zdiff <- (zdiff - zdiff[ , ,c(2:img.dims[3], 1)])
    z1 <- (zdiff > 0) * 1
    z1 <- z1[ , ,c(img.dims[3], 1:(img.dims[3]-1))]

    valleys <- ((x1 + y1 + z1) > 0) * 1
    valley.ls <- which(valleys == 1, arr.ind=TRUE)
    tneg.hills <- (tneg > 0)*1
    tneg.hills[valleys == 1] <- 0
    hill.ls <- which(tneg.hills == 1, arr.ind=TRUE)

    # fix holes in valleys
    for (i in 1:nrow(hill.ls)) {
      neighbors <- valleys[rbind(hill.ls[i, ] - c(1, 0, 0),
                                 hill.ls[i, ] + c(1, 0, 0))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }

      neighbors <- valleys[rbind(hill.ls[i, ] - c(0, 1, 0),
                                 hill.ls[i, ] + c(0, 1, 0))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }

      neighbors <- valleys[rbind(hill.ls[i, ] - c(0, 0, 1),
                                 hill.ls[i, ] + c(0, 0, 1))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }

      neighbors <- valleys[rbind(hill.ls[i, ] - c(1, 1, 0),
                                 hill.ls[i, ] + c(1, 1, 0))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }

      neighbors <- valleys[rbind(hill.ls[i, ] - c(1, 0, 1),
                                 hill.ls[i, ] + c(1, 0, 1))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }

      neighbors <- valleys[rbind(hill.ls[i, ] - c(0, 1, 1),
                                 hill.ls[i, ] + c(0, 1, 1))]
      if (sum(neighbors) == 2) {
        valleys[matrix(hill.ls[i, ], ncol=3)] <- 1
      }
    }

    tneg.hills <- tneg > 0
    tneg.hills[valleys == 1] <- 0
    hill.ls <- which(tneg.hills == 1, arr.ind=TRUE)
    valley.ls <- which(valleys == 1, arr.ind=TRUE)

    connected <- array(0, dim = img.dims[1:3])
    num.clusters <- 0
    idx <- numeric()
    for (x in 1:img.dims[1]) {
      for (y in 1:img.dims[2]) {
        for (z in 1:img.dims[3]) {
          if (tneg.hills[x, y, z] == 1) {
            num.clusters <- num.clusters + 1
            current.pt <- t(as.matrix(c(x, y, z)))
            idx <- as.integer(colSums(t(current.pt[, 1:3, drop = FALSE] - 1) * cdim) + 1L)
            connected[idx] <- num.clusters
            while (length(idx) != 0) {
              tneg.hills[idx] <- 0
              neighbors = as.vector(apply(as.matrix(idx), 1, "+", offsets))
              neighbors = unique(neighbors[which(neighbors > 0)])
              idx = neighbors[which(tneg.hills[neighbors] != 0)]
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
    fname <- paste0(save.dir, "/", file.prefix, ".neg.cluster.nii")
    init.nii(fname, img.dims, pixdim, orient)
    write.nii.volume(fname, 1, out)
  }
}
