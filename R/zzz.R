# Package lifecycle hooks. Loading the shared object is handled by
# useDynLib(gcmrec, .registration = TRUE) in NAMESPACE.

.onUnload <- function(libpath) {
  library.dynam.unload("gcmrec", libpath)
}
