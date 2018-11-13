cluster.3d <- function(in.array, connectivity=26L, min.size=10) {
# Check inputs ---------------------------------------------------------------
  stopifnot(length(dim(in.array))==3)
  stopifnot(connectivity %in% c(6,18,26))
  
# Setup search neighbourhood -------------------------------------------------
  array.dims <- dim(in.array)
  
  cdim <- cumprod(c(1, array.dims[-3]))
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

# Run clustering algorithm ---------------------------------------------------
  clusters <- array(0, dim=array.dims)
  num.clusters <- 0
  idx <- numeric()
  for (x in 1:array.dims[1]) {
    for (y in 1:array.dims[2]) {
      for (z in 1:array.dims[3]) {
        if (in.array[x,y,z] == 1) {
          num.clusters <- num.clusters + 1
          current.pt <- t(as.matrix(c(x,y,z)))
          idx <- as.integer(colSums(t(current.pt[ ,1:3, drop=FALSE]-1)*cdim) + 1L)
          clusters[idx] <- num.clusters
          while (length(idx)!=0) {
            in.array[idx] <- 0
            neighbors = as.vector(apply(as.matrix(idx), 1, '+', offsets))
            neighbors = unique(neighbors[which(neighbors > 0)])
            idx = neighbors[which(in.array[neighbors]!=0)]
            clusters[idx] <- num.clusters
          }
        }
      }
    }
  }
  cluster.table <- as.data.frame(table(clusters))
  cluster.table <- cluster.table[cluster.table$Freq > min.size, ]
  if (nrow(cluster.table) == 1) {
    stop("No clusters found")
  } else {
# renumber clusters by descending size order ---------------------------------
    cluster.ordered <- array(0, dim=array.dims)
    cluster.order <- order(cluster.table$Freq[2:nrow(cluster.table)],
                           decreasing = TRUE)
    for (i in 1:length(cluster.order)) {
      cluster.ordered[which(clusters == cluster.order[i])] <- i
    }
    
# send output ----------------------------------------------------------------
    out <- list()
    out$clusters <- cluster.ordered
    out$table <- as.data.frame(table(cluster.ordered))
    return(out)
  }
}