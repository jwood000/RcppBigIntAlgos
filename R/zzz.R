pkgEnv <- new.env(parent = emptyenv())
pkgEnv$nThreads <- NULL

## This will set the maximum number of cores
## and number of threads on a given machine
## when the package is loaded
.onLoad <- function(libname, pkgname) {
    tempThreads <- stdThreadMax()
    
    if (is.na(tempThreads)) {
        pkgEnv$nThreads <- 1L
    } else {
        pkgEnv$nThreads <- tempThreads
    }
    
    invisible()
}

