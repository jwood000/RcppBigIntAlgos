divisorsBig <- function(v, namedList = FALSE, showStats = FALSE,
                        skipPolRho = FALSE, skipECM = FALSE, nThreads = NULL) {
    
    GetDivisorsC(v, namedList, showStats, skipPolRho,
                 skipECM, nThreads, pkgEnv$nThreads)
}

quadraticSieve <- function(n, showStats = FALSE, nThreads = NULL) {
    QuadraticSieveContainer(n, showStats, nThreads, pkgEnv$nThreads)
}

primeFactorizeBig <- function(v, namedList = FALSE, showStats = FALSE,
                              skipPolRho = FALSE, skipECM = FALSE, nThreads = NULL) {
    
    PrimeFactorization(v, namedList, showStats, skipPolRho,
                       skipECM, nThreads, pkgEnv$nThreads)
}

stdThreadMax <- function() {
    nThreads <- cpp11GetNumThreads()
    
    if (nThreads < 1L || is.na(nThreads) || is.null(nThreads)) {
        nThreads = 1L
    }
    
    nThreads
}
