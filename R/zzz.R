.noGenerics <- TRUE

.onUnload <- function(libpath)
    library.dynam.unload("dna", libpath)
