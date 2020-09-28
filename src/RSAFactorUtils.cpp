#include "RSAFactorUtils.h"
#include "StatsUtils.h"
#include <limits>

// Max number of iterations in the main loop
constexpr std::size_t POLLARD_RHO_REPS = 100000u;

std::size_t GetPower(mpz_class &nMpz) {
    
    mpz_class testRoot;
    
    std::size_t p = 2;
    std::size_t myPow = 1;
    
    for (std::size_t i = 0; i < primesDiffPR.size(); ) {
        mpz_root(testRoot.get_mpz_t(), nMpz.get_mpz_t(), p);
        mpz_pow_ui(testRoot.get_mpz_t(), testRoot.get_mpz_t(), p);
        
        if (cmp(testRoot, nMpz) == 0) {
            std::size_t powPrime = p;
            bool keepGoing = true;
            
            // Determine if root is a power of a prime
            while (keepGoing) {
                powPrime *= p;
                mpz_root(testRoot.get_mpz_t(), nMpz.get_mpz_t(), powPrime);
                mpz_pow_ui(testRoot.get_mpz_t(), testRoot.get_mpz_t(), powPrime);
                
                if (cmp(testRoot, nMpz) != 0)
                    keepGoing = false;
            }
            
            powPrime /= p;
            myPow *= powPrime;
            mpz_root(nMpz.get_mpz_t(), nMpz.get_mpz_t(), powPrime);
        }
        
        p += primesDiffPR[i++];
        if (!mpz_perfect_power_p(nMpz.get_mpz_t()))
            break;
    }
    
    // This means the powers involved are very large
    if (mpz_perfect_power_p(nMpz.get_mpz_t())) {
        mpz_class myNextP = static_cast<int>(p);
        
        for (;;) {
            mpz_root(testRoot.get_mpz_t(), nMpz.get_mpz_t(), p);
            mpz_pow_ui(testRoot.get_mpz_t(), testRoot.get_mpz_t(), p);
            
            if (cmp(testRoot, nMpz) == 0) {
                std::size_t powPrime = p;
                bool keepGoing = true;
                
                // Determine if root is a power of a prime
                while (keepGoing) {
                    powPrime *= p;
                    mpz_root(testRoot.get_mpz_t(), nMpz.get_mpz_t(), powPrime);
                    mpz_pow_ui(testRoot.get_mpz_t(), testRoot.get_mpz_t(), powPrime);
                    
                    if (cmp(testRoot, nMpz) != 0)
                        keepGoing = false;
                }
                
                powPrime /= p;
                myPow *= powPrime;
                mpz_root(nMpz.get_mpz_t(), nMpz.get_mpz_t(), powPrime);
            }
            
            mpz_nextprime(myNextP.get_mpz_t(), myNextP.get_mpz_t());
            p = myNextP.get_ui();
            
            if (!mpz_perfect_power_p(nMpz.get_mpz_t()))
                break;
        }
    }
    
    return myPow;
}

void PollardRhoWithConstraint(mpz_class &n, unsigned long int a, std::vector<mpz_class> &factors,
                              std::vector<std::size_t>& myLens, std::size_t myLimit,
                              std::size_t powMultiplier) {
    
    mpz_class x, z, y, p, t;
    y = x = z = 2;
    p = 1;
    
    std::size_t k = 1u;
    std::size_t q = 1u;
    std::size_t count = 0;
    
    while (cmp(n, 1) != 0) {
        for (;;) {
            do {
                x *= x;
                x %= n;
                x += a;
                
                t = z - x;
                mpz_mod(t.get_mpz_t(), t.get_mpz_t(), n.get_mpz_t());
                p *= t;
                p %= n;
                
                if (k % 32 == 1) {
                    t = gcd(p, n);
                    
                    if (cmp(t, 1) != 0)
                        goto factor_found;
                    
                    y = x;
                }
                
                ++count;
            } while (--k != 0 && count < myLimit);
            
            if (count == myLimit)
                goto myReturn;
            
            z = x;
            k = q;
            q <<= 1;
            
            for (std::size_t i = 0; i < k; ++i) {
                x *= x;
                x %= n;
                x += a;
            }
            
            y = x;
        }
        
    factor_found:
        do {
            y *= y;
            y %= n;
            y += a;
            t = gcd(z - y, n);
        } while (cmp(t, 1) == 0);
        
        n /= t;	/* divide by t, before t is overwritten */

        if (mpz_probab_prime_p(t.get_mpz_t(), MR_REPS) == 0) {
            PollardRhoWithConstraint(t, a + 1, factors,
                                     myLens, myLimit, powMultiplier);
        } else {
            factors.push_back(t);
            myLens.push_back(powMultiplier);
            
            while (mpz_divisible_p(n.get_mpz_t(), t.get_mpz_t())) {
                n /= t;
                myLens.back() += powMultiplier;
            }
        }
        
        if (mpz_probab_prime_p(n.get_mpz_t(), MR_REPS) != 0) {
            factors.push_back(n);
            n = 1;
            myLens.push_back(powMultiplier);
            break;
        }
        
        x %= n;
        z %= n;
        y %= n;
    }
    
myReturn:
    // If reach this point, we have hit the limit for the number of 
    // iterations that want to try for Pollard Rho... the rest of
    // the factorization will be up to the Quadratic Sieve. Also,
    // we must have a statement here, so we simply set x = n.
    x = n;
}

void GetBigPrimeFacs(mpz_class &n, std::vector<mpz_class> &factors,
                     std::vector<mpz_class> &result,
                     std::vector<std::size_t>& myLens, 
                     std::size_t nThreads, bool bShowStats,
                     std::size_t powMaster) {
    
    if (mpz_sizeinbase(n.get_mpz_t(), 10) < 24) {
        PollardRhoWithConstraint(n, 1, factors, myLens, 
                                 std::numeric_limits<std::size_t>::max(),
                                 powMaster);
    } else {
        QuadraticSieve(n, result, nThreads, bShowStats);
        
        for (std::size_t i = 0; i < 2; ++i) {
            const std::size_t myPow = ((mpz_perfect_power_p(result[i].get_mpz_t())) ?
                                        GetPower(result[i]) : 1) * powMaster;
            
            if (mpz_probab_prime_p(result[i].get_mpz_t(), MR_REPS) == 0) {
                std::vector<mpz_class> recurseTemp(2);
                
                if (bShowStats) {
                    Rcpp::Rcout << "\nSummary Statistics for Factoring:\n" << "    "
                                << result[i].get_str() << "\n" << std::endl;
                }
                
                GetBigPrimeFacs(result[i], factors, recurseTemp,
                                myLens, nThreads, bShowStats, myPow);
            } else {
                n /= result[i];
                factors.push_back(result[i]);
                myLens.push_back(myPow);

                while (mpz_divisible_p(n.get_mpz_t(), result[i].get_mpz_t()))
                    n /= result[i];
            }
        }
    }
}

void QuadSieveHelper(mpz_class &nMpz, std::vector<mpz_class> &factors,
                     std::vector<std::size_t> &lengths, std::size_t nThreads,
                     bool bShowStats, bool bSkipExtPR) {
    
    const auto t0 = std::chrono::steady_clock::now();
    
    // First we test for small factors.
    TrialDivision(nMpz, factors, lengths);
    
    if (bShowStats) {
        Rcpp::Rcout << "\nSummary Statistics for Factoring:\n" << "    "
                    << nMpz.get_str() << "\n" << std::endl;
    }
    
    if (cmp(nMpz, 1) != 0) {
        // We now test for larger primes using pollard's rho
        // algorithm, but constrain it to a limited number of checks
        PollardRhoWithConstraint(nMpz, 1, factors, lengths,
                                 POLLARD_RHO_REPS, 1);
        
        if (bShowStats) {
            Rcpp::Rcout << "|  Pollard Rho Time  |\n|--------------------|" << std::endl;
            OneColumnStats(std::chrono::steady_clock::now() - t0);
        }
        
        if (cmp(nMpz, 1) != 0) {
            std::vector<mpz_class> result(2);
            
            // Shield quadratic sieve from perfect powers
            const std::size_t myPow = (mpz_perfect_power_p(nMpz.get_mpz_t())) ? GetPower(nMpz) : 1;
            
            if (mpz_probab_prime_p(nMpz.get_mpz_t(), MR_REPS) != 0) {
                factors.push_back(nMpz);
                lengths.push_back(myPow);
            } else {
                const int digCount = mpz_sizeinbase(nMpz.get_mpz_t(), 10);
                
                // We add additional iterations for very large numbers. Sometimes,
                // these numbers have a disproportionately smaller prime factor and
                // can be factorized faster with Pollard's rho algo. The numbers
                // below were obtain empirically.
                const std::size_t adder = std::min((digCount > 70) ? (digCount - 70)
                                                       * 2500000 : 0, 100000000);
                
                if (adder > 0 && !bSkipExtPR) {
                    PollardRhoWithConstraint(nMpz, 1, factors, lengths,
                                             POLLARD_RHO_REPS + adder, 1);
                }
                
                if (bShowStats) {
                    OneColumnStats(std::chrono::steady_clock::now() - t0);
                    Rcpp::Rcout << "\n" << std::endl;
                }
                
                GetBigPrimeFacs(nMpz, factors, result, lengths,
                                nThreads, bShowStats, myPow);
            }
        }
    }
    
    if (bShowStats) {
        Rcpp::Rcout << "|     Total Time     |\n|--------------------|" << std::endl;
        OneColumnStats(std::chrono::steady_clock::now() - t0);
        Rcpp::Rcout << "\n" << std::endl;
    }
}
