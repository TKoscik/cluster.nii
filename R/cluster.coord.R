cluster.coord <- function(cluster.nii,
                          cluster.value = "all",
                          value.nii = NULL,
                          value.vol = NULL,
                          mm = TRUE) {

  ### debug
  # cluster.nii <- "/Shared/harshmanl/ckd_bids/derivatives/dwi/analyses/FA-boysOnly/dwi-FA_FA-boysOnly_coef-fdr_vol2_cl18_ppeak0.001_pplateau0.01_sz100.neg.nii"
  # cluster.value = "all"
  # value.nii <- "/Shared/harshmanl/ckd_bids/derivatives/dwi/analyses/FA-boysOnly/dwi-FA_FA-boysOnly_coef-tvalue.nii"
  # value.vol <- 2
  ###
  
  cluster <- read.nii.volume(cluster.nii, 1)
  cluster.table <- as.data.frame(table(cluster))[-1, ]
  colnames(cluster.table) <- c("Cluster Value", "Number of Voxels")
  
  # Setup output container ----
  cluster.table$`Geometric Center X` <- numeric(nrow(cluster.table)) * NA
  cluster.table$`Geometric Center Y` <- numeric(nrow(cluster.table)) * NA
  cluster.table$`Geometric Center Z` <- numeric(nrow(cluster.table)) * NA
  cluster.table$`Center of Gravity X` <- numeric(nrow(cluster.table)) * NA
  cluster.table$`Center of Gravity Y` <- numeric(nrow(cluster.table)) * NA
  cluster.table$`Center of Gravity Z` <- numeric(nrow(cluster.table)) * NA
  cluster.table$`Peak X` <- numeric(nrow(cluster.table)) * NA
  cluster.table$`Peak Y` <- numeric(nrow(cluster.table)) * NA
  cluster.table$`Peak Z` <- numeric(nrow(cluster.table)) * NA
  cluster.table$`Peak Value` <- numeric(nrow(cluster.table)) * NA
  
  # Get World Coordinates ----
  if (!is.null(value.nii)) {
    values <- read.nii.volume(value.nii, value.vol)
  }
  
  for (i in 1:nrow(cluster.table)) {
    cluster.coords <- which(cluster == cluster.table$`Cluster Value`[i], arr.ind = TRUE)
    
    cluster.table$`Geometric Center X`[i] <- mean(cluster.coords[ ,1])
    cluster.table$`Geometric Center Y`[i] <- mean(cluster.coords[ ,2])
    cluster.table$`Geometric Center Z`[i] <- mean(cluster.coords[ ,3])
    
    if (!is.null(value.nii)) {
      cluster.values <- values[cluster.coords]
      cluster.table$`Center of Gravity X`[i] <- weighted.mean(cluster.coords[ ,1], cluster.values)
      cluster.table$`Center of Gravity Y`[i] <- weighted.mean(cluster.coords[ ,2], cluster.values)
      cluster.table$`Center of Gravity Z`[i] <- weighted.mean(cluster.coords[ ,3], cluster.values)
      
      cluster.table$`Peak X`[i] <- cluster.coords[which(abs(cluster.values) == max(abs(cluster.values)))[1],1]
      cluster.table$`Peak Y`[i] <- cluster.coords[which(abs(cluster.values) == max(abs(cluster.values)))[1],2]
      cluster.table$`Peak Z`[i] <- cluster.coords[which(abs(cluster.values) == max(abs(cluster.values)))[1],3]
      cluster.table$`Peak Value`[i] <- cluster.values[which(abs(cluster.values) == max(abs(cluster.values)))[1]]
    }
  }

  # Convert to MM ----
  if (mm) {
    # Get transform ----
    tform <- nii.hdr(nii.file = cluster.nii,
                     field = c("srow_x", "srow_y", "srow_z"))
    tform <- rbind(tform$srow_x, tform$srow_y, tform$srow_z, c(0,0,0,1))
    
    cluster.table[ ,3:5] <- t(tform %*% rbind(cluster.table$`Geometric Center X`-1,
                                              cluster.table$`Geometric Center Y`-1,
                                              cluster.table$`Geometric Center Z`-1,
                                              rep(1, nrow(cluster.table))))[ ,-4]
    cluster.table[ ,6:8] <- t(tform %*% rbind(cluster.table$`Center of Gravity X`-1,
                                              cluster.table$`Center of Gravity Y`-1,
                                              cluster.table$`Center of Gravity Z`-1,
                                              rep(1, nrow(cluster.table))))[ ,-4]
    cluster.table[ ,9:11] <- t(tform %*% rbind(cluster.table$`Peak X`-1,
                                              cluster.table$`Peak Y`-1,
                                              cluster.table$`Peak Z`-1,
                                              rep(1, nrow(cluster.table))))[ ,-4]
  }

  # Output ----
  return(cluster.table)
}
