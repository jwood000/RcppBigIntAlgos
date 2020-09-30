#include "SieveUtils.h"

std::vector<std::uint8_t> MyIntToBit(std::size_t x, std::size_t dig) {
    
    std::vector<std::uint8_t> binaryVec(dig);
    
    for (std::size_t i = 0; x > 0; ++i) {
        binaryVec[i] = x % 2;
        x >>= 1;
    }
    
    return binaryVec;
}


