#include "CppConvert.h"

#include "PrimeFactorUtils.h"
#include "StatsUtils.h"
#include "LenstraECM.h"
#include <numeric>

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
                mpz_pow_ui(
                    testRoot.get_mpz_t(), testRoot.get_mpz_t(), powPrime
                );

                if (cmp(testRoot, nMpz) != 0) {
                    keepGoing = false;
                }
            }

            powPrime /= p;
            myPow *= powPrime;
            mpz_root(nMpz.get_mpz_t(), nMpz.get_mpz_t(), powPrime);
        }

        p += primesDiffPR[i++];

        if (!mpz_perfect_power_p(nMpz.get_mpz_t())) {
            break;
        }
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
                    mpz_pow_ui(
                        testRoot.get_mpz_t(), testRoot.get_mpz_t(), powPrime
                    );

                    if (cmp(testRoot, nMpz) != 0) {
                        keepGoing = false;
                    }
                }

                powPrime /= p;
                myPow *= powPrime;
                mpz_root(nMpz.get_mpz_t(), nMpz.get_mpz_t(), powPrime);
            }

            mpz_nextprime(myNextP.get_mpz_t(), myNextP.get_mpz_t());
            p = myNextP.get_ui();

            if (!mpz_perfect_power_p(nMpz.get_mpz_t())) {
                break;
            }
        }
    }

    return myPow;
}

void PollardRhoWithConstraint(
    mpz_class &n, std::uint32_t a, std::vector<mpz_class> &factors,
    std::vector<std::size_t> &myLens, std::size_t myLimit,
    std::size_t powMultiplier
) {

    mpz_class x(2);
    mpz_class z(2);
    mpz_class y(2);
    mpz_class p(1);
    mpz_class t;
    mpz_class t2;

    std::size_t k = 1u;
    std::size_t q = 1u;
    std::size_t count = 0u;

    while (cmp(n, 1) != 0) {
        for (;;) {
            do {
                x = (x * x) % n + a;
                t = z - x;

                // Need to guarantee that t is positive so we must use mpz_mod
                mpz_mod(t.get_mpz_t(), t.get_mpz_t(), n.get_mpz_t());
                p *= t;
                p %= n;

                if (k % 32 == 1) {
                    t = gcd(p, n);

                    if (cmp(t, 1) != 0) {
                        goto factor_found;
                    }

                    y = x;
                }

                ++count;
            } while (--k != 0 && count < myLimit);

            if (count == myLimit) {
                goto myReturn;
            }

            z = x;
            k = q;
            q <<= 1;

            for (std::size_t i = 0; i < k; ++i) {
                x = (x * x) % n + a;
            }

            y = x;
        }

    factor_found:
        do {
            y = (y * y) % n + a;
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

void QuadraticSieveRecurse(mpz_class &n, std::vector<mpz_class> &factors,
                           std::vector<mpz_class> &results,
                           std::vector<std::size_t> &myLens,
                           std::size_t nThreads, bool bShowStats,
                           std::size_t powMaster) {

    if (mpz_sizeinbase(n.get_mpz_t(), 10) < 24) {
        PollardRhoWithConstraint(n, 1, factors, myLens,
                                 std::numeric_limits<std::size_t>::max(),
                                 powMaster);
    } else {
        QuadraticSieve(n, results, nThreads, bShowStats);

        for (std::size_t i = 0; i < 2; ++i) {
            const std::size_t myPow = (
                (mpz_perfect_power_p(results[i].get_mpz_t())) ?
                                        GetPower(results[i]) : 1
            ) * powMaster;

            if (mpz_probab_prime_p(results[i].get_mpz_t(), MR_REPS) == 0) {
                std::vector<mpz_class> recurseTemp(2);

                if (bShowStats) {
                    std::string res = "\nSummary Statistics for Factoring:\n"
                        "    " + results[i].get_str() + "\n";
                    Rprintf("%s\n", res.c_str());
                }

                QuadraticSieveRecurse(results[i], factors, recurseTemp,
                                      myLens, nThreads, bShowStats, myPow);
            } else {
                n /= results[i];
                factors.push_back(results[i]);
                myLens.push_back(myPow);

                while (mpz_divisible_p(n.get_mpz_t(),
                                       results[i].get_mpz_t())) {
                    n /= results[i];
                }
            }
        }
    }
}

void LenstraRecurse(mpz_class &n, std::vector<mpz_class> &factors,
                    std::vector<mpz_class> &results,
                    std::vector<mpz_class> &notFactored,
                    std::vector<std::size_t> &myLens,
                    const std::vector<unsigned long int> &primes,
                    std::size_t nThreads, bool bShowStats,
                    std::size_t powMaster, std::size_t totalCurves,
                    typeTimePoint checkPoint0) {

    if (mpz_sizeinbase(n.get_mpz_t(), 10) < 24) {
        PollardRhoWithConstraint(n, 1, factors, myLens,
                                 std::numeric_limits<std::size_t>::max(),
                                 powMaster);
    } else {
        std::size_t numCurves = 0;
        const std::size_t nDig = mpz_sizeinbase(n.get_mpz_t(), 10);
        const auto myKey = CurveLookup.upper_bound(nDig);
        const std::size_t maxLoopIter = myKey->second;
        const bool factor_found = LenstraECM(n, maxLoopIter, primes,
                                             results, numCurves, nThreads);

        if (bShowStats) {
            totalCurves += numCurves;
            TwoColumnStats(
                std::chrono::steady_clock::now() - checkPoint0,
                totalCurves, 0, false
            );
        }

        if (factor_found) {
            for (std::size_t i = 0; i < 2; ++i) {
                const std::size_t myPow = (
                    (mpz_perfect_power_p(results[i].get_mpz_t())) ?
                                               GetPower(results[i]) : 1
                ) * powMaster;

                if (mpz_probab_prime_p(results[i].get_mpz_t(), MR_REPS) == 0) {
                    std::vector<mpz_class> recurseTemp(2);
                    LenstraRecurse(
                        results[i], factors, recurseTemp, notFactored,
                        myLens, primes, nThreads, bShowStats, myPow,
                        totalCurves, checkPoint0
                    );
                } else {
                    n /= results[i];
                    factors.push_back(results[i]);
                    myLens.push_back(myPow);

                    while (mpz_divisible_p(
                            n.get_mpz_t(), results[i].get_mpz_t()
                           )) {
                        n /= results[i];
                    }
                }
            }
        } else {
            notFactored.push_back(n);
        }
    }
}

void FactorECM(mpz_class &n, std::vector<mpz_class> &factors,
               std::vector<mpz_class> &notFactored,
               std::vector<std::size_t> &myLens,
               std::size_t nThreads, bool bShowStats,
               std::size_t powMaster) {

    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t nDig = mpz_sizeinbase(n.get_mpz_t(), 10);

    const auto myKey = CurveLookup.upper_bound(nDig);
    const std::size_t maxLoopIter = myKey->second;
    const unsigned long int maxCurves = GetMaxCurves(maxLoopIter);
    const std::vector<unsigned long int> primes = GenerateNPrimes(maxCurves);

    std::size_t numCurves = 0;
    std::vector<mpz_class> results(2);

    if (bShowStats) {
        Rprintf("|  Lenstra ECM Time  |  Number of Curves  |\n"
                "|--------------------|--------------------|\n");
        TwoColumnStats(std::chrono::steady_clock::now() - t0, 0, 0, false);
    }

    LenstraRecurse(n, factors, results, notFactored, myLens, primes,
                   nThreads, bShowStats, powMaster, numCurves, t0);
}

void QuadSieveHelper(mpz_class &nMpz, std::vector<mpz_class> &factors,
                     std::vector<std::size_t> &lengths, std::size_t nThreads,
                     bool bShowStats, bool bSkipPR, bool bSkipECM) {

    const auto t0 = std::chrono::steady_clock::now();

    // First we test for small factors.
    TrialDivision(nMpz, factors, lengths);

    if (bShowStats) {
        std::string res = "\nSummary Statistics for Factoring:\n"
            "    " + nMpz.get_str() + "\n";
        Rprintf("%s\n", res.c_str());
    }

    if (nMpz > 1) {
        // We now test for larger primes using pollard's rho
        // algorithm, but constrain it to a limited number of checks
        PollardRhoWithConstraint(nMpz, 1, factors, lengths,
                                 POLLARD_RHO_REPS, 1);

        if (bShowStats && !bSkipPR) {
            Rprintf("|  Pollard Rho Time  |\n|--------------------|\n");
            OneColumnStats(std::chrono::steady_clock::now() - t0);
        }

        if (nMpz > 1) {
            // Shield quadratic sieve from perfect powers
            const std::size_t myPow = (
                mpz_perfect_power_p(nMpz.get_mpz_t())
            ) ? GetPower(nMpz) : 1;

            if (mpz_probab_prime_p(nMpz.get_mpz_t(), MR_REPS) != 0) {
                factors.push_back(nMpz);
                lengths.push_back(myPow);

                if (bShowStats) {
                    Rprintf("\n\n");
                }
            } else {
                if (!bSkipPR) {
                    const int digCount = mpz_sizeinbase(nMpz.get_mpz_t(), 10);

            // We add additional iterations for very large numbers. Sometimes,
            // these numbers have a disproportionately smaller prime factor and
            // can be factorized faster with Pollard's rho algo. The numbers
            // below were obtain empirically.
                    const std::size_t adder = std::min(
                        (digCount > 70) ? (digCount - 70) * 80000 : 0, 2000000
                    );

                    if (adder > 0) {
                        PollardRhoWithConstraint(nMpz, 1, factors, lengths,
                                                 POLLARD_RHO_REPS + adder, 1);
                    }

                    if (bShowStats) {
                        OneColumnStats(std::chrono::steady_clock::now() - t0);
                        Rprintf("\n\n");
                    }
                }

                if (bSkipECM) {
                    std::vector<mpz_class> results(2);
                    QuadraticSieveRecurse(
                        nMpz, factors, results, lengths,
                        nThreads, bShowStats, myPow
                    );
                } else {
                    std::vector<mpz_class> notFactored;
                    FactorECM(nMpz, factors, notFactored,
                              lengths, nThreads, bShowStats, myPow);

                    if (bShowStats) {
                        Rprintf("\n\n");
                    }

                    for (auto n: notFactored) {
                        std::vector<mpz_class> results(2);
                        QuadraticSieveRecurse(n, factors, results, lengths,
                                              nThreads, bShowStats, myPow);
                    }
                }
            }
        } else {
            if (bShowStats) {
                Rprintf("\n\n");
            }
        }
    }

    if (bShowStats) {
        Rprintf("|     Total Time     |\n|--------------------|\n");
        OneColumnStats(std::chrono::steady_clock::now() - t0);
        Rprintf("\n\n");
    }
}

SEXP PrimeFactorizeHuge(mpz_class &nMpz, std::size_t nThreads,
                        bool bShowStats, bool bSkipPR, bool bSkipECM) {

    if (sgn(nMpz) == 0) {
        cpp11::stop("Cannot factorize 0");
    }

    const std::size_t IsNegative = (sgn(nMpz) < 0) ? 1 : 0;

    if (IsNegative) {
        nMpz = abs(nMpz);
    }

    if (cmp(nMpz, 1) == 0) {
        if (IsNegative) {
            mpz_class mpzNegOne = -1;
            cpp11::sexp myNegOne = Rf_allocVector(RAWSXP, intSize * 4);

            char* r = (char*) (RAW(myNegOne));
            ((int*) (r))[0] = 1;

            CppConvert::rawExport(
                &r[intSize], mpzNegOne.get_mpz_t(), intSize * 3
            );

            myNegOne.attr("class") = "bigz";
            return myNegOne;
        } else {
            cpp11::sexp noPrimeFacs = Rf_allocVector(RAWSXP, intSize);
            Rbyte* rawPt = RAW(noPrimeFacs);

            for (std::size_t i = 0; i < intSize; ++i) {
                rawPt[i] = 0;
            }

            noPrimeFacs.attr("class") = "bigz";
            return noPrimeFacs;
        }
    }

    std::vector<std::size_t> lengths;
    std::vector<mpz_class> factors;

    QuadSieveHelper(nMpz, factors, lengths, nThreads,
                    bShowStats, bSkipPR, bSkipECM);

    // Sort the prime factors as well as order the
    // lengths vector by the order of the factors array
    CppConvert::QuickSort(factors, 0, factors.size() - 1, lengths);
    const std::size_t totalNum = std::accumulate(
        lengths.cbegin(), lengths.cend(), 0u
    ) + IsNegative;

    std::size_t size = intSize;
    std::vector<std::size_t> mySizes(totalNum);

    mpz_class negOne = -1;

    if (IsNegative) {
        const std::size_t tempSize = intSize *
            (2 + (mpz_sizeinbase(negOne.get_mpz_t(), 2) + numb - 1) / numb);
        size += tempSize;
        mySizes[0] = tempSize;
    }

    // adding each bigint's needed size
    for (std::size_t i = 0, count = IsNegative; i < factors.size(); ++i) {
        for (std::size_t j = 0; j < lengths[i]; ++j, ++count) {
            const std::size_t tempSize = intSize *
                (2 + (mpz_sizeinbase(factors[i].get_mpz_t(),
                                     2) + numb - 1) / numb);
            size += tempSize;
            mySizes[count] = tempSize;
        }
    }

    cpp11::sexp ans = Rf_allocVector(RAWSXP, size);
    char* r = (char*) (RAW(ans));
    ((int*) (r))[0] = totalNum; // first int is vector-size-header

    // current position in pos[] (starting after vector-size-header)
    std::size_t pos = intSize;

    if (IsNegative) {
        pos += CppConvert::rawExport(
            &r[pos], negOne.get_mpz_t(), mySizes.front()
        );
    }

    for (std::size_t i = 0, count = IsNegative; i < factors.size(); ++i) {
        for (std::size_t j = 0; j < lengths[i]; ++j, ++count) {
            pos += CppConvert::rawExport(
                &r[pos], factors[i].get_mpz_t(), mySizes[count]
            );
        }
    }

    ans.attr("class") = "bigz";
    return ans;
}
