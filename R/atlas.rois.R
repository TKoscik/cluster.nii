atlas.rois <- function(cluster.nii,
                       cluster.vol="all",
                       atlas.nii,
                       atlas.csv,
                       return.table = TRUE,
                       save.table = FALSE,
                       save.dir = NULL,
                       prefix = NULL) {
  
  #debug 
  rm(list=ls())
  gc()
  cluster.nii  <- "/rdss/koscikt/projects/duplex_gamble/analyses/20171220_decision/decision_mediators.vol1.cl26.th.gt0.sz15.nii"
  cluster.vol="all"
  atlas.nii <- "/rdss/koscikt/brains/cluster_splitter_2mm.nii"
  atlas.csv <- "/rdss/koscikt/brains/cluster_splitter_key.csv"
  return.table = TRUE
  save.table = TRUE
  save.dir <- "/rdss/koscikt/projects/duplex_gamble/analyses/20171220_decision"
  prefix <- "test.table"
  #
  
  n.cluster <- nii.dims(cluster.nii)[4]
  if (cluster.vol=="all") { cluster.vol <- 1:n.cluster }

  atlas.nii <- read.nii.volume(atlas.nii, 1)
  atlas.csv <- read.csv(atlas.csv, header=TRUE, stringsAsFactors=FALSE)

  cluster.id <- numeric(0L)
  cluster.csv <- data.frame(cluster.vol = numeric(0), atlas.csv[0, ])
  for (i in cluster.vol) {
    temp.nii <- read.nii.volume(cluster.nii, i)
    atlas.nums <- atlas.nii[temp.nii==1]
    atlas.ls <- unique(atlas.nums)
    which.atlas <- atlas.ls[which.max(tabulate(match(atlas.nums, atlas.ls)))]
    which.row <- which(atlas.csv$atlas.label == which.atlas)
    cluster.csv <- rbind(cluster.csv, data.frame(cluster.vol=i, atlas.csv[which.row, ]))
  }

  if (save.table) {
    write.table(cluster.csv, file=paste0(save.dir, "/", prefix, ".csv"),
                row.names=FALSE, sep=",", quote=FALSE)
  }

  if (return.table) {
    return(cluster.csv)
  }
}
