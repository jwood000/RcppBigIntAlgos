#include <vector>
#include <gmpxx.h>
#include <thread>
#include <string>
#include <cmath>

constexpr std::size_t maxIterPerThrd = 70;
constexpr std::size_t minIterPerThrd = 10;

unsigned long int GetMaxCurves(std::size_t maxLoopIter) {

    unsigned long int B1 = 8;
    unsigned long int B2 = 13;

    for (std::size_t i = 0; i < maxLoopIter; ++i) {
        unsigned long int B2_temp = B2;
        B2 += B1;
        B1 = B2_temp;
    }

    return B2;
}

unsigned long int UpperBoundEst(double tar) {

    double myMax = tar;
    double myMin = tar;
    double x = tar;

    while ((x / std::log(x)) < tar) {
        myMin = x;
        x *= x;
        myMax = x;
    }

    double size = myMax - myMin;
    double size2 = std::trunc(size / 2);
    double temp = myMin + size2;
    double dist = (tar - (temp / std::log(temp)));
    double lower = 0;
    double upper = 0;

    if (dist > 0) {
        lower = myMin + size2 - 1;
        upper = myMin + size;
    } else {
        lower = myMin + 1;
        upper = myMin + size2 + 1;
    }

    while (size2 > 1 && !(dist == 0)) {
        size2 = std::trunc((upper - lower) / 2);
        temp = lower + size2;
        dist = (tar - (temp / std::log(temp)));

        if (dist > 0) {
            lower = temp - 1;
        } else {
            upper = temp + 1;
        }
    }

    return temp;
}

std::vector<unsigned long int> GenerateNPrimes(std::size_t limit) {

    std::size_t upperBound = UpperBoundEst(limit);
    std::vector<char> primes(upperBound + 1, 1);
    std::vector<unsigned long int> myps(limit);

    const std::size_t fsqr = std::sqrt(static_cast<double>(upperBound));

    for (std::size_t j = 4u; j <= upperBound; j += 2) {
        primes[j] = 0;
    }

    for (std::size_t lastP = 3u; lastP <= fsqr;) {
        for (std::size_t j = lastP * lastP; j <= upperBound; j += 2 * lastP) {
            primes[j] = 0;
        }

        lastP += 2;

        while (!primes[lastP]) {
            lastP += 2;
        }
    }

    myps.front() = 2u;

    for (unsigned long int j = 3u, i = 1; i < limit; j += 2) {
        if (primes[j]) {
            myps[i++] = j;
        }
    }

    return myps;
}

void ecm_add(const mpz_class &x1, const mpz_class &x2, const mpz_class &y1,
             const mpz_class &y2, const mpz_class &z1, const mpz_class &z2,
             const mpz_class &n, mpz_class &t1, mpz_class &t2,
             mpz_class &t3, mpz_class &t4) {

    t1 = (x1 - x2) * (y1 + y2);
    t2 = (x1 + x2) * (y1 - y2);
    t3 = t1 + t2;
    t4 = t1 - t2;

    t1 = (z2 * t3 * t3) % n;

    if (t1 < 0) {
        t1 += n;
    }

    t2 = (z1 * t4 * t4) % n;

    if (t2 < 0) {
        t2 += n;
    }
}

void ecm_double(const mpz_class &x1, const mpz_class &x2,
                const mpz_class &c1, const mpz_class &c2,
                const mpz_class &n, mpz_class &t1,
                mpz_class &t2, mpz_class &t3, mpz_class &t4) {

    t1 = x1 + x2;
    t1 *= t1;

    t2 = x1 - x2;
    t2 *= t2;

    t3 = t1 - t2;
    t4 = 4 * c2 * t2;

    t1 = (t1 * t4) % n;

    if (t1 < 0) {
        t1 += n;
    }

    t2 = (t3 * (t4 + c1 * t3)) % n;

    if (t2 < 0) {
        t2 += n;
    }
}

void ecm_multiply(const mpz_class &k, const mpz_class &z1,
                  const mpz_class &z2, const mpz_class &c1,
                  const mpz_class &c2, mpz_class &x1,
                  mpz_class &x2, const mpz_class &n, mpz_class &t1,
                  mpz_class &t2, mpz_class &t3, mpz_class &t4,
                  mpz_class &q1, mpz_class &q2, mpz_class &b,
                  mpz_class &test) {
    q1 = x1;
    q2 = x2;

    ecm_double(x1, x2, c1, c2,
               n, t1, t2, t3, t4);

    x1 = t1;
    x2 = t2;

    b = k / 2;
    test = -b;
    test &= b;

    while (cmp(b, test) > 0) {
        b ^= test;
        test = -b;
        test &= b;
    }

    for (; cmp(b, 0) > 0; b /= 2) {
        test = b & k;

        if (cmp(test, 0)) {
            ecm_add(x1, x2, q1, q2, z1,
                    z2, n, t1, t2, t3, t4);
            q1 = t1;
            q2 = t2;

            ecm_double(x1, x2, c1, c2,
                       n, t1, t2, t3, t4);
            x1 = t1;
            x2 = t2;
        } else {
            ecm_add(x1, x2, q1, q2, z1,
                    z2, n, t1, t2, t3, t4);
            x1 = t1;
            x2 = t2;

            ecm_double(q1, q2, c1, c2,
                       n, t1, t2, t3, t4);
            q1 = t1;
            q2 = t2;
        }
    }
}

void InnerLoop(std::size_t strt, std::size_t myEnd, const mpz_class &k,
               mpz_class &g, const mpz_class &n, std::vector<char> &res,
               int idx) {

    mpz_class u, v, x, z, b, t;
    mpz_class x1, x2, c1, c2, q1, q2;
    mpz_class t1, t2, t3, t4, t5;

    for (unsigned long int s = strt; s < myEnd; ++s) {
        u = s * s - 5;
        v = 4 * s;

        x = u * u * u;
        z = v * v * v;
        t = v - u;

        mpz_powm_ui(t.get_mpz_t(), t.get_mpz_t(),
                    3, n.get_mpz_t());

        c1 = (t * (3 * u + v)) % n;
        c2 = (4 * x * v) % n;

        x1 = x;
        x2 = z;

        ecm_multiply(k, x1, x2, c1, c2, x, z,
                     n, t1, t2, t3, t4, q1, q2, b, t5);
        g = gcd(q2, n);

        if (g > 1 && g < n) {
            res[idx] = 1;
            return void();
        }
    }

    res[idx] = 0;
}

void SetThreadsIters(int sectionLength,
                     std::size_t maxThreads,
                     std::size_t &nThrdsThisIter,
                     std::size_t &chunk) {

    if (sectionLength > static_cast<int>(maxThreads * maxIterPerThrd)) {
        nThrdsThisIter = maxThreads;
        chunk = maxIterPerThrd;
    } else if (sectionLength >
               static_cast<int>(maxThreads * maxThreads * minIterPerThrd)) {
        nThrdsThisIter = maxThreads;
        chunk = (sectionLength + nThrdsThisIter - 1) / nThrdsThisIter;
    } else {
        for (std::size_t thrd = 1; thrd <= maxThreads; ++thrd) {
            const std::size_t test = minIterPerThrd * thrd;

            if ((sectionLength / test) < thrd) {
                nThrdsThisIter = thrd;
                chunk = (sectionLength + thrd - 1) / thrd;
                break;
            }
        }
    }
}

bool LenstraECM(const mpz_class &n, std::size_t maxLoopIter,
                const std::vector<unsigned long int> &primes,
                std::vector<mpz_class> &factors,
                std::size_t &numCurves, std::size_t nThreads) {

    unsigned long int B1 = 8;
    unsigned long int B2 = 13;
    bool factor_found = false;
    mpz_class k, g;

    for (std::size_t i = 0, j = 0;
         i < maxLoopIter && !factor_found; ++i, j = 0) {
        double p = primes[j++];
        k = 1;

        while (p < B1) {
            const double mypow = std::trunc(std::log(B1) / std::log(p));
            const double temp = std::trunc(std::pow(p, mypow));
            k *= temp;
            p = primes[j++];
        }

        std::size_t strt = B1;
        int sectionLength = B2 - B1;
        numCurves = B1;

        std::size_t nThrdsThisIter = nThreads;
        std::size_t chunk = sectionLength / nThreads;
        SetThreadsIters(sectionLength, nThreads, nThrdsThisIter, chunk);

        while (sectionLength > 0 && nThrdsThisIter > 1) {
            std::vector<std::thread> myThreads;
            std::vector<mpz_class> vecFactors(nThrdsThisIter);
            std::vector<char> vecSuccess(nThrdsThisIter);
            std::size_t myEnd = strt + chunk;

            for (std::size_t thrd = 0; thrd < (nThrdsThisIter - 1);
                 ++thrd, strt = myEnd, myEnd += chunk) {

                myThreads.emplace_back(
                    std::cref(InnerLoop), strt, myEnd, std::cref(k),
                    std::ref(vecFactors[thrd]), std::cref(n),
                    std::ref(vecSuccess), thrd
                );
            }

            const std::size_t lastEnd = (myEnd > B2) ? B2 : myEnd;

            myThreads.emplace_back(
                std::cref(InnerLoop), strt, myEnd, std::cref(k),
                std::ref(vecFactors.back()), std::cref(n),
                std::ref(vecSuccess), nThrdsThisIter - 1
            );

            for (auto &thr: myThreads) {
                thr.join();
            }

            numCurves += ((nThrdsThisIter - 1) * chunk + lastEnd - strt);
            sectionLength -= ((nThrdsThisIter - 1) * chunk + lastEnd - strt);
            strt = myEnd;

            const bool bSuccess = std::any_of(
                vecSuccess.begin(), vecSuccess.end(), [](char vS) {return vS;}
            );

            if (bSuccess) {
                for (std::size_t thrd = 0; thrd < nThrdsThisIter; ++thrd) {
                    if (vecSuccess[thrd]) {
                        factors[0] = vecFactors[thrd];
                        factors[1] = n / vecFactors[thrd];
                        factor_found = true;
                    }
                }
            }

            SetThreadsIters(sectionLength, nThreads, nThrdsThisIter, chunk);
        }

        if (sectionLength > 0) {
            std::vector<char> res(1);
            InnerLoop(strt, B2, k, g, n, res, 0);

            if (res.front()) {
                factors[0] = g;
                factors[1] = n / g;
                factor_found = true;
            }
        }

        double B_temp = B2;
        B2 += B1;
        B1 = B_temp;
    }

    return factor_found;
}
