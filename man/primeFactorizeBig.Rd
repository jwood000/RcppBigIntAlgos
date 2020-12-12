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
  \item{showStats}{Logical flag for showing summary statistics. The default is \code{FALSE}.}
  \item{skipPolRho}{Logical flag for skipping the extended pollard rho algorithm. The default is \code{FALSE}.}
  \item{skipECM}{Logical flag for skipping the extended elliptic curve algorithm. The default is \code{FALSE}.}
  \item{nThreads}{Number of threads to be used for the elliptic curve method and the quadratic sieve.s The default is \code{NULL}.}
}
\details{
This function should be preferred in most situations and is identical to \code{\link{quadraticSieve}} when both \code{skipECM} and  \code{skipPolRho} are set to \code{TRUE}.

It is optimized for factoring big and small numbers by dynamically using different algorithms based off of the input. It takes cares to not spend too much time in any of the methods and avoids wastefully switching to the quadratic sieve when the number is very large.

See \code{\link{quadraticSieve}} for information regarding \code{showStats}.
}
\value{
\itemize{
\item{Returns an unnamed vector of class bigz if \code{length(v) == 1} regardless of the value of \code{namedList}.}
\item{If \code{length(v) > 1}, a named/unnamed list of vectors of class bigz will be returned.}
}
}
\references{
\itemize{
\item{\href{https://www.iacr.org/archive/ches2006/10/10.pdf}{Gaj K. et al. (2006) Implementing the Elliptic Curve Method of Factoring in Reconfigurable Hardware. In: Goubin L., Matsui M. (eds) Cryptographic Hardware and Embedded Systems - CHES 2006. CHES 2006. Lecture Notes in Computer Science, vol 4249. Springer, Berlin, Heidelberg. https://doi.org/10.1007/11894063_10}}
\item{\href{https://en.wikipedia.org/wiki/Integer_factorization}{Integer Factorization}}
}
}
\author{
Joseph Wood
}
\note{
Note, the function \code{primeFactorizeBig(n, skipECM = T, skipPolRho = T)} is the same as \code{quadraticSieve(n)}
}
\seealso{
\code{\link[RcppAlgos]{primeFactorize}}, \code{\link[numbers]{primeFactors}}, \code{\link[gmp:factor]{factorize}}, \code{\link{quadraticSieve}}
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