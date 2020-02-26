cluster.atlas <- function(cluster.nii, atlas.nii, atlas.key) {
  
  # debug
  # cluster.nii <- "/Shared/harshmanl/ckd_bids/derivatives/dwi/analyses/FA-boysOnly/dwi-FA_FA-boysOnly_coef-fdr_vol2_cl18_ppeak0.001_pplateau0.01_sz100.neg.nii"
  # atlas.nii <- c("/Shared/harshmanl/ckd_bids/derivatives/dwi/analyses/HCPICBM_1mm_label-DKT+SUIT.nii.gz",
  #                "/Shared/nopoulos/nimg_core/templates_human/HCPICBM/1mm/HCPICBM_1mm_label-JHU.nii.gz")
  # atlas.key <- c("/Shared/harshmanl/ckd_bids/derivatives/dwi/analyses/DKT+SUIT_key.csv",
  #                "/Shared/nopoulos/nimg_core/templates_human/HCPICBM/1mm/key-JHU.csv")
  # library(tools)
  # library(R.utils)
  # library(nifti.io)
  #
  
  n.atlas <- length(atlas.nii)
  atlas <- vector("list", n.atlas)
  key <- vector("list", n.atlas)
  for (i in 1:n.atlas) {
    if (file_ext(atlas.nii[i]) == "gz") {
      tname = gunzip(atlas.nii[i], remove=FALSE, temporary=TRUE, overwrite=TRUE)
      atlas[[i]] <- read.nii.volume(tname, 1)
    } else {
      atlas[[i]] <- read.nii.volume(atlas.nii[i], 1)
    }
    key[[i]] <- read.csv(atlas.key[i], stringsAsFactors = F)
    if (!(0 %in% key[[i]]$value)) {
      key[[i]][nrow(key[[i]])+1, ] = c(0, "Unlabelled", "unlabelled")
    }
  }
  
  cluster <- read.nii.volume(cluster.nii, 1)
  temp.table <- as.data.frame(table(cluster))[-1, ]
  n.cluster <- nrow(temp.table)
  cluster.table <- data.frame(Cluster.Value = temp.table$cluster,
                              Atlas = matrix("", nrow=n.cluster, ncol=2),
                              stringsAsFactors = FALSE)
  
  for (i in 1:n.cluster) {
    cluster.idx <- which(cluster == cluster.table$Cluster.Value[i], arr.ind=TRUE)
    for (j in 1:n.atlas) {
      atlas.vals <- as.data.frame(table(round(atlas[[j]][cluster.idx])))
      atlas.vals$roi <- as.character(key[[j]]$label[key[[j]]$value %in% atlas.vals$Var1])
      atlas.vals$pct <- atlas.vals$Freq / sum(atlas.vals$Freq) *100
      atlas.vals$string <- paste0(round(atlas.vals$pct, digits=2), "% ", atlas.vals$roi)
      cluster.table[i,j+1] <- paste(atlas.vals$string[order(atlas.vals$pct, decreasing = T)], collapse=", ")
    }
  }
  
  return(cluster.table)
}
