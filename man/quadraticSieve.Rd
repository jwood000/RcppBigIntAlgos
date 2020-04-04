\name{quadraticSieve}
\alias{quadraticSieve}
\title{
Prime Factorization with the Quadratic Sieve
}
\description{
Get the prime factorization of a number, \eqn{n}, using the \href{https://en.wikipedia.org/wiki/Quadratic_sieve}{Quadratic Sieve}.
}
\usage{
quadraticSieve(n, showStats = FALSE, nThreads = NULL)
}
\arguments{
  \item{n}{An integer, numeric, string value, or an element of class bigz.}
  \item{showStats}{Logical flag. If \code{TRUE}, summary statistics will be displayed.}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
}
\details{
First, \href{https://en.wikipedia.org/wiki/Trial_division}{trial division} is carried out to remove small prime numbers, then a modified version of \href{https://en.wikipedia.org/wiki/Pollard\%27s_rho_algorithm}{Pollard's rho algorithm} that is constrained is called to quickly remove further prime numbers. Next, we check to make sure that we are not passing a perfect power to the main quadratic sieve algorithm. After removing any perfect powers, we finally call the quadratic sieve with multiple polynomials in a recursive fashion until we have completely factored our number.
}
\value{
Vector of class bigz
}
\references{
\itemize{
\item{\href{http://library.msri.org/books/Book44/files/03carl.pdf}{Pomerance, C. (2008). Smooth numbers and the quadratic sieve. In Algorithmic
       Number Theory Lattices, Number Fields, Curves and Cryptography (pp. 69-81). 
       Cambridge: Cambridge University Press.}}
\item{\href{http://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866119-8/S0025-5718-1987-0866119-8.pdf}{Silverman, R. D. (1987). The Multiple Polynomial Quadratic Sieve. 
               Mathematics of Computation, 48(177), 329-339. doi:10.2307/2007894}}
\item{\href{http://micsymposium.org/mics_2011_proceedings/mics2011_submission_28.pdf}{Integer Factorization using the Quadratic Sieve}}
\item{\href{https://codegolf.stackexchange.com/questions/8629/fastest-semiprime-factorization/9088#9088}{From https://codegolf.stackexchange.com/ (Credit to user primo for answer)
       P., & Chowdhury, S. (2012, October 7). Fastest semiprime factorization. Retrieved October 06, 2017}}
}
}
\author{
Joseph Wood
}
\seealso{
\code{\link[gmp]{factorize}}
}
\examples{
mySemiPrime <- prod(nextprime(urand.bigz(2, 40, 17)))
quadraticSieve(mySemiPrime)
}
