\name{divisorsBig}
\alias{divisorsBig}
\title{
Vectorized Factorization (Complete) with GMP
}
\description{
Quickly generates the complete factorization for many (possibly large) numbers.
}
\usage{
divisorsBig(v, namedList = FALSE)
}
\arguments{
  \item{v}{Vector of integers, numerics, string values, or elements of class bigz.}
  \item{namedList}{Logical flag. If \code{TRUE} and the \code{length(v) > 1}, a named list is returned. The default is \code{FALSE}.}
}
\details{
Highly optimized algorithm to generate the complete factorization for many numbers. It is built specifically for the data type that is used in the gmp library (i.e. \href{https://gmplib.org/manual/Nomenclature-and-Types.html#Nomenclature-and-Types}{\code{mpz_t}}).

The main part of this algorithm is essentially the same algorithm that is implemented in \link[RcppAlgos]{divisorsRcpp} from the RcppAlgos package. A modified \href{https://en.wikipedia.org/wiki/Merge_sort}{merge sort} algorithm is implemented to better deal with the \code{mpz_t} data type. This algorithm avoids directly swapping elements of the main factor array of type \code{mpz_t} but instead generates a vector of indexing integers for ordering. 

See this stack overflow post for examples and benchmarks : \href{https://stackoverflow.com/a/49742904/4408538}{R Function for returning ALL factors}.
}
\value{
\itemize{
\item{Returns an unnamed vector of class bigz if \code{length(v) == 1} regardless of the value of \code{namedList}.}
\item{If \code{length(v) > 1}, a named/unnamed list of vectors of class bigz will be returned.}
}
}
\references{
\href{https://en.wikipedia.org/wiki/Divisor}{Divisor}
}
\author{
Joseph Wood
}

\seealso{
\code{\link[RcppAlgos]{divisorsRcpp}}, \code{\link[numbers]{divisors}}, \code{\link[gmp]{factorize}}
}
\examples{
## Get the complete factorization of a single number
divisorsBig(100)

## Or get the complete factorization of many numbers
set.seed(29)
myVec <- sample(-1000000:1000000, 1000)
system.time(myFacs <- divisorsBig(myVec))

## Return named list
myFacsWithNames <- divisorsBig(myVec, namedList = TRUE)
}