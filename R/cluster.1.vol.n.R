cluster.1.vol.n <- function(nii.file, save.dir=NULL, file.name=NULL) {
  if (missing(save.dir)) { save.dir <- getwd() }
  if (missing(file.name)) {
    ftemp <- unlist(strsplit(nii.file, split="[/]"))
    ftemp <- ftemp[length(ftemp)]
    ftemp <- paste0(ftemp, collapse="")
    ftemp <- unlist(strsplit(ftemp, split="[.]"))
    ftemp <- ftemp[-length(ftemp)]
    ftemp <- paste0(c(ftemp, "nvol.nii"), collapse=".")
    file.name <- ftemp
  }

  nii <- read.nii.volume(nii.file, 1)
  values <- unique(as.vector(nii))
  if (0 %in% values) { values <- values[-(values==0)] }
  img.dims <- c(nii.dims(nii.file)[1:3], length(values))
  if (prod(img.dims > .Machine$integer.max)) {stop("array size exceeds maximum integer index size in R")}

  init.nii(file.name = paste(save.dir, file.name, sep="/"),
           dims = img.dims,
           pixdim = unlist(nii.hdr(nii.file, "pixdim")),
           orient = nii.orient(nii.file))
  for (i in 1:img.dims[4]) {
    write.nii.volume(nii.file = paste(save.dir, file.name, sep="/"),
                     vol.num = i,
                     values = (nii==values[i]) * 1)
  }
}
