#include "Cpp14MakeUnique.h"
#include <cstdlib>
#include <memory>
#include <gmp.h>

void Grow(std::unique_ptr<mpz_t[]> &myArray,
          std::size_t oldMax, std::size_t newMax) {
    
    auto tempFacBase = FromCpp14::make_unique<mpz_t[]>(newMax);
    
    for (std::size_t i = 0; i < oldMax; ++i) {
        mpz_init(tempFacBase[i]);
        mpz_set(tempFacBase[i], myArray[i]);
        mpz_clear(myArray[i]);
    }
    
    for (std::size_t i = oldMax; i < newMax; ++i)
        mpz_init(tempFacBase[i]);
    
    myArray.reset(tempFacBase.release());
}