# Data sets copied from the package drcData (Ritz, Streibig and Gerhard, GPL-2), which is not on CRAN.
# Run from the package root to recreate the files in data/.
sets <- c("chickweed", "daphnids", "earthworms", "germination", "lettuce", "liver.tumor", "ryegrass",
          "S.alba", "selenium", "TCDD", "vinclozolin")
for (d in sets){
  e <- new.env()
  data(list=d, package="drcData", envir=e)
  save(list=d, envir=e, file=file.path("data", paste0(d, ".rda")), compress="xz", version=2)
}
