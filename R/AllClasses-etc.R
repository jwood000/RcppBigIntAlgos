### Really want  S4 methods for  all the binary operations.
### Otherwise we "never" make
###       <bigz> o  <bigq>
### or    <bigq> o  <bigz>    working --
##
## But unfortunately the above seems "impossible", see
##  see also  setMethod() in ./matrix-prods.R

## OTOH: This *still* helps to define single-dispatch methods for  asNumeric() :
##       {why does it work there ??}

setOldClass("bigz")#, prototype=as.bigz(integer()))
##                cannot use as.bigz() yet which is only defined in ./bigz.R
