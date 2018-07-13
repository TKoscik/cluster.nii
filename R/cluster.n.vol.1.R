cluster.n.vol.1 <- function(nii.file, save.dir=NULL, file.name=NULL) {
  if (is.null(save.dir)) { save.dir <- getwd() }
  if (is.null(file.name)) {
    ftemp <- unlist(strsplit(nii.file, split="[/]"))
    ftemp <- ftemp[length(ftemp)]
    ftemp <- paste0(ftemp, collapse="")
    ftemp <- unlist(strsplit(ftemp, split="[.]"))
    ftemp <- ftemp[-length(ftemp)]
    ftemp <- paste0(c(ftemp, "1vol.nii"), collapse=".")
    file.name <- ftemp
  }

  img.dims <- nii.dims(nii.file)

  new.nii <- array(0, dim=img.dims[1:3])
  for (i in 1:img.dims[4]) {
    new.nii <- new.nii + read.nii.volume(nii.file, i) * i
  }

  init.nii(file.name = paste0(save.dir, "/", file.name),
           dims = img.dims[1:3],
           pixdim = unlist(nii.hdr(nii.file, "pixdim")),
           orient = nii.orient(nii.file))
  write.nii.volume(nii.file = paste0(save.dir, "/", file.name),
                   vol.num = 1,
                   values = new.nii)
}
