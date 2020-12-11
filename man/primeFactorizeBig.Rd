\name{primeFactorizeBig}
\alias{primeFactorizeBig}
\title{
Vectorized Prime Factorization with GMP
}
\description{
Quickly generates the prime factorization for many (possibly large) numbers, using trial division, \href{https://en.wikipedia.org/wiki/Pollard\%27s_rho_algorithm}{Pollard's rho algorithm}, \href{https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization}{Lenstra's Elliptic Curve method}, and finally the \href{https://en.wikipedia.org/wiki/Quadratic_sieve}{Quadratic Sieve}.
}
\usage{
primeFactorizeBig(v, namedList = FALSE, showStats = FALSE,
                  skipPolRho = FALSE, skipECM = FALSE, nThreads = NULL)
}
\arguments{
  \item{v}{Vector of integers, numerics, string values, or elements of class bigz.}
  \item{namedList}{Logical flag. If \code{TRUE} and the \code{length(v) > 1}, a named list is returned. The default is \code{FALSE}.}
  \item{showStats}{Logical flag passed to \code{\link{quadraticSieve}} (if needed). If \code{TRUE}, summary statistics will be displayed.}
  \item{skipPolRho}{Logical flag for skipping the extended pollard rho algorithm. The default is \code{FALSE}.}
  \item{skipECM}{Logical flag for skipping the extended elliptic curve algorithm. The default is \code{FALSE}.}
  \item{nThreads}{Number of threads to be used for \code{\link{quadraticSieve}} (if needed.) The default is \code{NULL}.}
}
\details{
Highly optimized algorithm to generate the prime factorization for many numbers.
}
\value{
\itemize{
\item{Returns an unnamed vector of class bigz if \code{length(v) == 1} regardless of the value of \code{namedList}.}
\item{If \code{length(v) > 1}, a named/unnamed list of vectors of class bigz will be returned.}
}
}
\references{
\href{https://en.wikipedia.org/wiki/Integer_factorization}{Integer Factorization}
}
\author{
Joseph Wood
}
\note{
Note, the function \code{primeFactorizeBig(n, skipECM = T, skipPolRho = T)} is the same as \code{quadraticSieve(n)}
}
\seealso{
\code{\link[RcppAlgos]{primeFactorize}}, \code{\link[numbers]{primeFactors}}, \code{\link[gmp:factor]{factorize}}
}
\examples{
## Get the prime factorization of a single number
primeFactorizeBig(100)

## Or get the prime factorization of many numbers
set.seed(29)
myVec <- sample(-1000000:1000000, 1000)
system.time(myFacs <- primeFactorizeBig(myVec))

## Return named list
myFacsWithNames <- primeFactorizeBig(myVec, namedList = TRUE)
}