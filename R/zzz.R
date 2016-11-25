

############ First.lib ###############

.onLoad <- function(lib, pkg){
   library.dynam("gcmrec", pkg, lib)
}

.onUnload <- function(libpath)
    library.dynam.unload("gcmrec", libpath)


############ End of .First.lib ###############



