#ifndef TONELLISHANKS_GMP_R
#define TONELLISHANKS_GMP_R 1
#include "bigvec.h"

/**
 * Obtain quadratic residues of a mod p
 * a: number to get square root of
 * p: prime modulus
 * quadRes [out]: the list of residues
 * 
 * note: a^((p-1)/2) mod p must be 1
 * i.e. the Legendre symbol of a with respect to p must be 1
 */
void TonelliShanksC (mpz_t a, mpz_t p, bigvec & quadRes);

#endif
