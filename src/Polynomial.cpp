#include "Polynomial.h"
#include "SolutionSearch.h"
#include "StatsUtils.h"
#include <fstream>
#include <set>

#include <Rcpp.h>

Polynomial::Polynomial(std::size_t _mpzContainerSize, std::size_t _facSize,
                       std::size_t _SaPThresh, std::size_t _polyLimit,
                       bool _bShowStats, mpz_class myNum) : 
            mpzFacSize(_facSize), SaPThresh(_SaPThresh),
            polyLimit(_polyLimit), bShowStats(_bShowStats) {
    
    powsOfSmooths.reserve(_SaPThresh);
    powsOfPartials.reserve(_SaPThresh);
    myStart.assign(_facSize * 2, 0u);

    smoothInterval.reserve(_mpzContainerSize);
    largeCoFactors.reserve(_mpzContainerSize);
    partialInterval.reserve(_mpzContainerSize);

    nPolys = 0;
    nPartial = 0;
    nSmooth = 0;
    coFactorInd = 0;
    
    if (bShowStats) {
        RcppThread::Rcout << "\nSummary Statistics for Factoring:\n" << "    "
                          << myNum.get_str() << "\n\n"
                          << "|        Time        | Complete | Polynomials |   Smooths"
                          << "  |  Partials  |\n";
        RcppThread::Rcout << "|--------------------|----------|-------------|---------"
                          << "---|------------|" << std::endl;
    }
}

void Polynomial::SievePolys(const std::vector<std::size_t> &SieveDist,
                            const std::vector<std::size_t> &facBase, const std::vector<int> &LnFB,
                            const std::vector<mpz_class> &NextPrime,
                            mpz_class LowBound, mpz_class myNum, int theCut, int DoubleLenB,
                            int vecMaxSize, std::size_t strt) {

    for (; nPolys <= polyLimit; ++nPolys) {
        ++mpzFacSize;

        SinglePoly(SieveDist, facBase, LnFB, powsOfSmooths, powsOfPartials,
                   coFactorIndexVec, myStart, partFactorsMap, partIntvlMap,
                   keepingTrack, smoothInterval, largeCoFactors,
                   partialInterval, NextPrime[nPolys], LowBound, myNum,
                   nPartial, nSmooth, coFactorInd, theCut, DoubleLenB,
                   mpzFacSize, vecMaxSize, strt);
    }
}

void Polynomial::FactorFinish(const std::vector<std::size_t> &SieveDist,
                              const std::vector<std::size_t> &facBase, const std::vector<int> &LnFB,
                              std::vector<mpz_class> &mpzFacBase, mpz_class NextPrime,
                              mpz_class LowBound, mpz_class myNum,
                              int theCut, int DoubleLenB, int vecMaxSize, std::size_t strt,
                              std::chrono::time_point<std::chrono::steady_clock> checkPoint0) {
    
    if (SaPThresh - facBase.size() > 45) {
        Rcpp::print(Rcpp::wrap("dawg"));
    }
    
    std::size_t currLim = nSmooth + nPartial;
    std::size_t lastLim = currLim;
    const std::size_t facSize = facBase.size();
    
    auto checkPoint1 = std::chrono::steady_clock::now();
    auto checkPoint2 = checkPoint1;
    auto checkPoint3 = checkPoint1;
    auto showStatsTime = (checkPoint1 - std::chrono::steady_clock::now());

    while (currLim <= SaPThresh) {
        for (bool LegendreTest = true; LegendreTest; ) {
            mpz_nextprime(NextPrime.get_mpz_t(), NextPrime.get_mpz_t());

            if (mpz_legendre(myNum.get_mpz_t(), NextPrime.get_mpz_t()) == 1)
                LegendreTest = false;
        }

        mpzFacBase.push_back(NextPrime);
        ++mpzFacSize;

        SinglePoly(SieveDist, facBase, LnFB, powsOfSmooths, powsOfPartials,
                   coFactorIndexVec, myStart, partFactorsMap, partIntvlMap,
                   keepingTrack, smoothInterval, largeCoFactors,
                   partialInterval, NextPrime, LowBound, myNum,
                   nPartial, nSmooth, coFactorInd, theCut, DoubleLenB,
                   mpzFacSize, vecMaxSize, strt);

        ++nPolys;
        currLim = nSmooth + nPartial;

        const auto checkPoint4 = std::chrono::steady_clock::now();

        if ((checkPoint4 - checkPoint2) > checkInterTime) {
            // Check for user interrupt and udpate timepoint
            RcppThread::checkUserInterrupt();
            checkPoint2 = std::chrono::steady_clock::now();
        }

        if (bShowStats && (checkPoint4 - checkPoint3) > showStatsTime) {
            MakeStats(currLim, SaPThresh, nPolys, nSmooth,
                      nPartial, lastLim, checkPoint4 - checkPoint0);

            checkPoint3 = std::chrono::steady_clock::now();
            UpdateStatTime(currLim, facSize, checkPoint4 - checkPoint0, showStatsTime);
        }
    }
    
    SaPThresh += 10;
}

void Polynomial::GetSolution(const std::vector<mpz_class> &mpzFacBase,
                             const std::vector<std::size_t> &facBase, mpz_t *const factors,
                             mpz_t mpzNum, std::size_t nThreads,
                             std::chrono::time_point<std::chrono::steady_clock> checkPoint0) {

    // Not every prime in mpzFacBase is utilized. For example, with our
    // attempt at factoring rsa99, there were over 4 million elements
    // in mpzFacBase, however there were only ~40000 smooths + partials.
    // This caused an inordinate memory allocation haulting the
    // factorization. The set, setIndex, will be used to drastically
    // reduce the size of the factor base and related structures making
    // the linear algebra portion much more efficient.
    std::set<int> setIndex;
    auto indIt = setIndex.begin();
    const std::size_t facSize = facBase.size();

    for (std::size_t i = 0; i < nSmooth; ++i) {
        if (powsOfSmooths[i].front() > facSize) {
            setIndex.emplace_hint(indIt, powsOfSmooths[i].front());
            indIt = setIndex.end();
        }
    }

    // Here we get the 4 largest powers as we know these relate to
    // primes not contained in facBase.
    for (std::size_t i = 0; i < nPartial; ++i) {
        std::partial_sort(powsOfPartials[i].begin(), powsOfPartials[i].begin() + 4,
                          powsOfPartials[i].end(), std::greater<int>());

        if (powsOfPartials[i].front() > facSize)
            setIndex.insert(powsOfPartials[i].front());

        if (powsOfPartials[i][2] > facSize)
            setIndex.insert(powsOfPartials[i][2]);
    }

    const std::size_t nRows = nSmooth + nPartial;
    const std::size_t nonTrivSize = setIndex.size() + facSize;

    const std::size_t nCols = nonTrivSize + coFactorInd + 1;
    auto nonTrivialFacs = FromCpp14::make_unique<mpz_t[]>(nCols);

    for (std::size_t i = 0; i < facSize; ++i)
        mpz_init_set_ui(nonTrivialFacs[i], facBase[i]);

    std::unordered_map<int, std::size_t> mapIndex;
    mapIndex.reserve(setIndex.size());

    std::size_t setInd = facSize;

    for (auto s: setIndex) {
        mpz_init_set(nonTrivialFacs[setInd++], mpzFacBase[s - 1].get_mpz_t());
        mapIndex[s] = setInd;
    }

    for (std::size_t i = 0, j = setInd; i < coFactorInd; ++i, ++j)
        mpz_init_set(nonTrivialFacs[j], largeCoFactors[i].get_mpz_t());

    auto newTestInt = FromCpp14::make_unique<mpz_t[]>(nRows);
    std::vector<std::uint8_t> mat(nRows * nCols, 0u);

    for (std::size_t r = 0, row = 0; r < nSmooth; ++r, row += nCols) {

        // Remap the powers corresponding to powers not contained in facBase
        for (std::size_t j = 0; j < 2; ++j)
            powsOfSmooths[r][j] = mapIndex[powsOfSmooths[r][j]];

        for (const auto p: powsOfSmooths[r])
            ++mat[row + p];

        mpz_init_set(newTestInt[r], smoothInterval[r].get_mpz_t());
    }

    for (std::size_t i = 0, r = nSmooth,
         row = nCols * nSmooth; i < nPartial; ++i, ++r, row += nCols) {

        // Remap the powers corresponding to powers not contained in facBase
        for (std::size_t j = 0; j < 4; ++j)
            powsOfPartials[i][j] = mapIndex[powsOfPartials[i][j]];

        for (const auto p: powsOfPartials[i])
            ++mat[row + p];

        mat[row + nonTrivSize + coFactorIndexVec[i] + 1] = 2u;
        mpz_init_set(newTestInt[r], partialInterval[i].get_mpz_t());
    }

    SolutionSearch(mat, nRows, nCols, mpzNum, nonTrivialFacs.get(),
                   newTestInt.get(), factors, nThreads);

    for (std::size_t i = 0; i < nRows; ++i)
        mpz_clear(newTestInt[i]);

    for (std::size_t i = 0; i < nCols; ++i)
        mpz_clear(nonTrivialFacs[i]);

    nonTrivialFacs.reset();
    newTestInt.reset();

    if (bShowStats && mpz_cmp_ui(factors[0], 0)) {
        const auto checkPoint3 = std::chrono::steady_clock::now();
        std::size_t lastLim = nSmooth + nPartial;

        MakeStats(SaPThresh, SaPThresh, nPolys, nSmooth,
                  nPartial, lastLim, checkPoint3 - checkPoint0);

        RcppThread::Rcout << "\n" << std::endl;
    }
}
