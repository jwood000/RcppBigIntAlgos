divisorsBig <- function(v, namedList = FALSE, nThreads = NULL) {
    getDivisorsC(v, namedList, nThreads, pkgEnv$nThreads)
}

quadraticSieve <- function(n, showStats = FALSE, nThreads = NULL) {
    QuadraticSieveContainer(n, showStats, nThreads, pkgEnv$nThreads)
}

stdThreadMax <- function() {
    nThreads <- cpp11GetNumThreads()
    
    if (nThreads < 1L || is.na(nThreads) || is.null(nThreads)) {
        nThreads = 1L
    }
    
    nThreads
}