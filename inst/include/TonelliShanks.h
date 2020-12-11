#ifndef TONELLISHANKS_GMP
#define TONELLISHANKS_GMP

#include "GmpxxCopy.h"

// Obtain quadratic residues of a mod p
// a: number to get square root of
// p: prime modulus
// quadRes [out]: the list of residues
// 
// note: a^((p-1)/2) mod p must be 1
// i.e. the Legendre symbol of a with respect to p must be 1

void TonelliShanksC(const mpz_class &myNum, const mpz_class &p, mpz_class &TS_1);

#endif
