#include "SieveUtils.h"

int64_t int64_invert(int64_t n, int64_t p) {
    
    const int64_t const_p = p;
    int64_t x = 0;
    n %= p;
    
    for (int64_t u = 1; n; ) {
        int64_t temp = x - ((p / n) * u);
        x = u;
        u = temp;
        temp = p % n;
        p = n;
        n = temp;
    }
    
    return (x + const_p) % const_p;
}

int64_t GetPowIndex(const mpz_class &myNum, mpz_class &temp,
                    unsigned long int myR, int64_t primePow) {
    
    double myInv = int64_invert(myR + myR, primePow);
    temp = myNum - (static_cast<double>(myR) * 
                    static_cast<double>(myR));
    temp *= myInv;
    temp += myR;
    temp %= static_cast<unsigned long int>(primePow);
    return (temp.get_si());
}

std::vector<std::size_t> SetHenselLift(const std::vector<int> &facBase,
                                       const mpz_class &myNum,
                                       int minPrime, int TwiceLenB) {
    
    const std::size_t facSize = facBase.size();
    std::vector<std::size_t> Hensel;
    mpz_class p, TS_1, TS_2, temp;
    
    for (std::size_t i = 1; i < facSize; ++i) {
        double primePow = facBase[i];
        
        p = facBase[i];
        TonelliShanksC(myNum, p, TS_1);
        TS_2 = p - TS_1;
        unsigned long int r1 = TS_1.get_ui();
        unsigned long int r2 = TS_2.get_ui();
        
        while (primePow < minPrime) {
            r1 = GetPowIndex(myNum, temp, r1, primePow);
            r2 = GetPowIndex(myNum, temp, r2, primePow);
            primePow *= static_cast<double>(facBase[i]);
        }
        
        for (; primePow < TwiceLenB;) {
            r1 = GetPowIndex(myNum, temp, r1, primePow);
            r2 = GetPowIndex(myNum, temp, r2, primePow);
            Hensel.push_back(r1);
            Hensel.push_back(r2);
            primePow *= static_cast<double>(facBase[i]);
        }
    }
    
    return Hensel;
}

int int_invert(unsigned int n, unsigned int p) {
    
    int x = 0;
    
    for (int u = 1; n;) {
        int temp = x - static_cast<int>((p / n) * u);
        x = u;
        u = temp;
        temp = p % n;
        p = n;
        n = temp;
    }
    
    return x;
}

// Getting quadratic residues. See tonellishanks.cc for more details
std::vector<std::size_t> GetSieveDist(const std::vector<int> &facBase,
                                      const mpz_class &myNum) {
    
    const std::size_t facSize = facBase.size();
    std::vector<std::size_t> SieveDist(facSize);
    mpz_class p, TS_1;
    
    for (std::size_t i = 1; i < facSize; ++i) {
        p = facBase[i];
        TonelliShanksC(myNum, p, TS_1);
        SieveDist[i] = TS_1.get_ui();
    }
    
    return SieveDist;
}

std::vector<int> GetPrimesQuadRes(const mpz_class &myN, double LimB, double fudge1,
                                  double sqrLogLog, std::size_t myTarget) {
    
    const std::size_t uN = LimB;
    std::vector<char> primes(uN + 1, 1);
    std::vector<int> myps;
    
    myps.reserve(LimB * 2.0 / std::log(LimB));
    const std::size_t fsqr = std::floor(std::sqrt(LimB));
    
    for (std::size_t j = 4; j <= uN; j += 2)
        primes[j] = 0;
    
    for (std::size_t lastP = 3; lastP <= fsqr;) {
        for (std::size_t j = lastP * lastP; j <= uN; j += 2 * lastP)
            primes[j] = 0;
        
        std::size_t k = lastP + 2;
        std::size_t ind = 2;
        
        while (!primes[k]) {
            k += 2;
            ind += 2;
        }
        
        lastP += ind;
    }
    
    myps.push_back(2);
    mpz_class currP, nextP;
    
    for (int j = 3; j <= static_cast<int>(uN); j += 2) {
        if (primes[j]) {
            currP = j;
            
            if (mpz_legendre(myN.get_mpz_t(), currP.get_mpz_t()) == 1)
                myps.push_back(j);
        }
    }
    
    while (myps.size() < myTarget) {
        fudge1 += 0.005;
        LimB = std::exp((0.5 + fudge1) * sqrLogLog);
        
        currP = myps.back();
        mpz_nextprime(nextP.get_mpz_t(),
                      currP.get_mpz_t());
        
        while (cmp(nextP, LimB) < 0) {
            currP = nextP;
            mpz_nextprime(nextP.get_mpz_t(),
                          currP.get_mpz_t());
            
            if (mpz_legendre(myN.get_mpz_t(), currP.get_mpz_t()) == 1)
                myps.push_back(currP.get_si());
        }
    }
    
    // Ensure that the vecMaxSize is utilized most
    // efficiently based off the size of facBase
    if (myps.back() > (8 * L1Cache)) {
        double myDec = std::fmod(static_cast<double>(myps.back()) / static_cast<double>(L1Cache), 1.0);
        
        if (myDec < 0.2) {
            const int smallerTarget = (myps.back() / L1Cache) * L1Cache;
    
            while (myps.back() > smallerTarget) {
                myps.pop_back();
            }
        }
    }
    
    return myps;
}
