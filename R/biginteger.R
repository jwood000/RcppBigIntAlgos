setGeneric("asNumeric", useAsDefault = function(x) {
    if(is.numeric(x)) x else if(is.atomic(x)) {
        storage.mode(x) <- "numeric"; x }
    else as(x, "numeric")
})

#----------------------------------------------------------
#
#  Author        : Immanuel Scholz (immanuel.scholz@gmx.de)
#		   Technische Universitaet Dresden
#
#  Brief         : Stub to call the dll functions
#
#  Licence       : GPL
#
#----------------------------------------------------------

print.bigz <- function(x, quote = FALSE, initLine = is.null(modulus(x)), ...)
{
  if((n <- length(x)) > 0) {
    if(initLine) {
      cat("Big Integer ('bigz') ")
      kind <- if(!is.null(nr <- attr(x, "nrow")))
        sprintf("%d x %d matrix", nr, n/nr)
      else if(n > 1) sprintf("object of length %d", n) else ""
      cat(kind,":\n", sep="")
    }
    print(as.character(x), quote = quote, ...)
  }
  else
    cat("bigz(0)\n")
  invisible(x)
}

as.bigz <- function(a, mod = NA)
{
    if(is.null(mod)) mod <- NA
    .Call("biginteger_as", a, mod, PACKAGE = "bigIntegerAlgos")
}
## the .as*() functions are exported for Rmpfr
.as.bigz <- function(a, mod = NA) {
    .Call("biginteger_as", a, mod, PACKAGE = "bigIntegerAlgos")
}
..as.bigz <- function(a, mod = NA) .Call("biginteger_as", a, mod, PACKAGE = "bigIntegerAlgos")

.as.char.bigz <-
    as.character.bigz <-
    function(x, b = 10L, ...) .Call("biginteger_as_character", x, b, PACKAGE = "bigIntegerAlgos")


##' format() Numbers such as to distinguish  bigz, integer, double, mpfr, etc
formatN <- function(x, ...) UseMethod("formatN")
formatN.integer <- function(x, ...) paste0(as.character(x, ...), "L")
formatN.bigz    <- function(x, ...) {
    r <- as.character(x, ...)
    if(any(noMod <- is.null(modulus(x))))
        r[noMod] <- paste0(r[noMod],"_Z")
    r
}
formatN.double	<- function(x, ...) {
    r <- vapply(x, format, "", ...)
    if(any(intLike <- !grepl("[^-0-9]",r)))
        r[intLike] <- paste0(r[intLike],".")
    r
}
##' Default Method: Use the standard format() --- e.g. for complex
formatN.default <- function(x, ...) format(x, ...)


as.double.bigz  <- function(x,...) .Call("biginteger_as_numeric", x, PACKAGE = "bigIntegerAlgos")
as.integer.bigz <- function(x,...) .Call("biginteger_as_integer", x, PACKAGE = "bigIntegerAlgos")

.bigz2num <- function(x) {
    r <- .Call("biginteger_as_numeric", x, PACKAGE = "bigIntegerAlgos")
    if(!is.null(d <- dim(x))) dim(r) <- d
    r
}
setMethod("asNumeric", "bigz", .bigz2num)

length.bigz <- function(x) .Call("biginteger_length", x, PACKAGE = "bigIntegerAlgos")

"length<-.bigz"<- function(x, value) .Call("biginteger_setlength", x, value, PACKAGE = "bigIntegerAlgos")

modulus      <- function(a) UseMethod("modulus")
modulus.bigz <- function(a) attr(a, "mod")

`modulus<-`      <- function(a, value) UseMethod("modulus<-")
`modulus<-.bigz` <- function(a, value) as.bigz(a, value)

getDivisors <- function(n) .Call("getDivisorsC", n, PACKAGE = "bigIntegerAlgos")
multPolyQuadSieve <- function(n) .Call("QuadraticSieveContainer", n, PACKAGE = "bigIntegerAlgos")
quadResidue <- function(n, p) .Call("QuadraticResidueContainer", n, p, PACKAGE = "bigIntegerAlgos")
