neg.pwr.thresh.decay <- function(roi.nii,
                                 roi.vol = 1,
                                 center,
                                 point1,
                                 point2 = NULL,
                                 exponent = -1,
                                 img.dims = NULL,
                                 pixdim = NULL,
                                 orient = NULL,
                                 save.dir,
                                 file.name = NULL) {
  
  # Debug
  # rm(list=ls())
  # gc()
  # roi.nii = "D:/brains/binMNI.nii"
  # roi.vol = 1
  # center = c(60,30,20) # left cerebellum
  # point1 = c(60, 5)
  # point2 = c(10, 20)
  # exponent = -1
  # img.dims = NULL
  # pixdim = NULL
  # orient = NULL
  # save.dir = "D:/brains"
  # file.name = NULL
  #
  
  roi.mask <- read.nii.volume(roi.nii, roi.vol)
  vxl.ls <- which(roi.mask != 0 & !is.na(roi.mask), arr.ind = TRUE)
  
  center <- matrix(center, ncol=3)
  vxl.dist <- sqrt((center[ ,1] - vxl.ls[ ,1])^2 +
                     (center[ ,2] - vxl.ls[ ,2])^2 +
                     (center[ ,3] - vxl.ls[ ,3])^2)
  
  if (is.null(point2)) {
    a <- 1
  } else {
    a <- (point2[2] - point1[2]) / (point2[1] ^ exponent - point1[1] ^ exponent)
  }
  b <- point1[2] - a * point1[1] ^ exponent
  
  threshold <- a * vxl.dist ^ exponent + b
  
  if (is.null(img.dims)) {
    img.dims <- nii.dims(roi.nii)[1:3]
  }
  if (is.null(pixdim)) {
    pixdim <- unlist(nii.hdr(roi.nii, "pixdim"))
  }
  if (is.null(orient)) {
    orient <- nii.orient(roi.nii)
  }
  
  out <- array(0, img.dims)
  out[vxl.ls] <- threshold
  out[is.infinite(out)] <- NA
  
  if (is.null(file.name)) {
    file.name <- paste0("neg.pwr.thresh.x", center[1],
                        ".y", center[2],
                        ".z", center[3],
                        ".e", exponent,
                        ".a", a,
                        ".b", b,
                        ".nii")
  }
  
  init.nii(paste0(save.dir, "/", file.name), img.dims, pixdim, orient)
  write.nii.volume(paste0(save.dir, "/", file.name), 1, out)
}