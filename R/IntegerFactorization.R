divisorsBig <- function(v, namedList = FALSE, showStats = FALSE,
                        nThreads = NULL, skipExtPolRho = FALSE) {
    GetDivisorsC(v, namedList, nThreads, 
                 pkgEnv$nThreads, showStats, skipExtPolRho)
}

quadraticSieve <- function(n, showStats = FALSE, nThreads = NULL, skipExtPolRho = FALSE) {
    QuadraticSieveContainer(n, showStats, nThreads, pkgEnv$nThreads, skipExtPolRho)
}

stdThreadMax <- function() {
    nThreads <- cpp11GetNumThreads()
    
    if (nThreads < 1L || is.na(nThreads) || is.null(nThreads)) {
        nThreads = 1L
    }
    
    nThreads
}
