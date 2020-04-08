divisorsBig <- function(v, namedList = FALSE) {
    GetDivisorsC(v, namedList, NULL, pkgEnv$nThreads)
}

quadraticSieve <- function(n, showStats = FALSE) {
    QuadraticSieveContainer(n, showStats, NULL, pkgEnv$nThreads)
}

stdThreadMax <- function() {
    nThreads <- cpp11GetNumThreads()
    
    if (nThreads < 1L || is.na(nThreads) || is.null(nThreads)) {
        nThreads = 1L
    }
    
    nThreads
}