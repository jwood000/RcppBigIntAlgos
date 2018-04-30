/* Factoring with Multiple Polynomial Quadratic Sieve.
 
In addition to the references in the man file, the links below are very helpful:
    - 1:
       URL: http://www.cs.virginia.edu/crab/QFS_Simple.pdf
       author: Eric Landquist
       date: December 14, 2001
       title: The Quadratic Sieve Factoring Algorithm
    - 2:
       URL: https://blogs.msdn.microsoft.com/devdev/2006/06/19/factoring-large-numbers-with-quadratic-sieve/
       author: MSDN Archive
       date: June 19, 2006
       title: Factoring large numbers with quadratic sieve
    - 3:
       URL: http://www.math.colostate.edu/~hulpke/lectures/m400c/quadsievex.pdf
*/

#include <vector>
#include <algorithm>
#include <cmath>
#include <inttypes.h>
#include "tonellishanks.h"
#include "quadraticsieve.h"

typedef std::vector<int64_t> v1d;
typedef std::vector<v1d> v2d;
typedef std::vector<v2d> v3d;

const unsigned long int hundredMillion = (unsigned long int) 100000000;

static inline v1d outersect (v1d x, v1d y) {
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    unsigned long int lenX = x.size(), lenY = y.size();
    v1d v(lenX + lenY);
    v1d::iterator it, it2;
    it = std::set_difference(x.begin(), x.end(), y.begin(), y.end(), v.begin());
    it2 = std::set_difference(y.begin(), y.end(), x.begin(), x.end(), it);
    v.resize(it2 - v.begin());
    return v;
}

static void reduceMatrix (unsigned long int n1,
                          unsigned long int n2, v2d &nullMat, v1d &myCols) {
    int64_t temp;
    unsigned long int myMin = 0, myMax = 0;
    v1d::iterator it;
    v1d myOnes;

    for (std::size_t i = 0; i < n1; i++) {
        myOnes.reserve(n2);
        for (std::size_t j = myMax; j < n2; j++)
            if (nullMat[j][i] == 1)
                myOnes.push_back((int64_t) j);

        if (myOnes.size() > 0) {
            myMin = myOnes[0];
            if (myMin != myMax) {
                for (std::size_t j = 0; j < n1; j++) {
                    temp = nullMat[myMin][j];
                    nullMat[myMin][j] = nullMat[myMax][j];
                    nullMat[myMax][j] = temp;
                }
            }
            myOnes.erase(myOnes.begin());
            for (it = myOnes.begin(); it < myOnes.end(); it++)
                for (std::size_t j = 0; j < n1; j++)
                    nullMat[*it][j] = (nullMat[*it][j] + nullMat[myMax][j]) % 2;

            myMax++;
        }
        myOnes.clear();
    }

    unsigned long int newLen;
    
    if (myMax < n2 && myMax != 0) {
        for (std::size_t j = myMax; j < n2; j++)
            nullMat.erase(nullMat.begin() + myMax);

        newLen = nullMat.size();
    } else {
        newLen = 0;
    }

    bool allZero;
    
    if (newLen > 0) {
        unsigned long int k = 0;
        while (k < newLen) {
            allZero = true;
            for (std::size_t j = 0; j < n1; j++) {
                if (nullMat[k][j] != 0) {
                    allZero = false;
                    break;
                }
            }
            if (allZero) {
                nullMat.erase(nullMat.begin() + k);
                newLen--;
                continue;
            }
            if (nullMat[k][k] != 1) {
                for (std::size_t j = 0; j < n1; j++) {
                    if (nullMat[k][j] == 1) {
                        myMin = j;
                        break;
                    }
                }
                for (std::size_t j = 0; j < newLen; j++) {
                    temp = nullMat[j][k];
                    nullMat[j][k] = nullMat[j][myMin];
                    nullMat[j][myMin] = temp;
                }
                temp = myCols[k];
                myCols[k] = myCols[myMin];
                myCols[myMin] = temp;
            }
            k++;
        }
    }
}

static inline v1d myIntToBit (unsigned long int x, 
                              unsigned long int dig) {
    unsigned long int i = 0;
    v1d binaryVec(dig);
    
    while (x > 0) {
        binaryVec[i] = (int64_t) (x % 2);
        x /= 2;
        i++;
    }
    
    return binaryVec;
}

static void solutionSearch (v2d mat, mpz_t n, v1d FB,
                            mpz_t * test, mpz_t factors[]) {
    
    unsigned long int nr1 = mat.size(), numCol = mat[0].size(), r = 0;
    v2d nullMat;

    for (std::size_t j = 0; j < numCol; j++) {
        std::size_t i = 0;
        while ((i < nr1) && (mat[i][j] == 0)) {i++;}
        if (i < nr1) {
            nullMat.push_back(v1d(nr1, 0));
            for (std::size_t k = 0; k < nr1; k++)
                nullMat[r][k] = mat[k][j] % 2;
            r++;
        }
    }
    
    unsigned long int nr2 = nullMat.size();
    numCol = nr1;
    v1d myCols(numCol, 0);
    for (std::size_t i = 0; i < numCol; i++)
        myCols[i] = i;
    
    reduceMatrix (numCol, nr2, nullMat, myCols);
    unsigned long int newNrow = nullMat.size();

    unsigned long int tLen;
    v1d::iterator it, it2;
    v2d myList(numCol, v1d());
    v1d freeVariables, temp;
    freeVariables.reserve(numCol);

    if (numCol > newNrow) {
        for (std::size_t i = newNrow; i < numCol; i++)
            freeVariables.push_back(myCols[i]);
        
        for (it = freeVariables.begin(); it < freeVariables.end(); it++)
            myList[*it].push_back(*it);
    }

    bool allBiggerNewNrow;

    if (newNrow > 0) {
        for (std::size_t i = newNrow; i > 0; i--) {
            temp.reserve(numCol);
            for (std::size_t j = 0; j < numCol; j++)
                if (nullMat[i-1][j] == 1)
                    temp.push_back((int64_t) j);

            if (temp.size() == 1) {
                for (std::size_t j = 0; j < newNrow; j++)
                    nullMat[j][temp[0]] = 0;
                
                myList[myCols[i-1]].clear();
                myList[myCols[i-1]].push_back(0);
            } else {
                temp.clear();
                temp.reserve(numCol);
                allBiggerNewNrow = true;
                for (std::size_t j = i; j < numCol; j++) {
                    if (nullMat[i-1][j] == 1) {
                        temp.push_back((int64_t) j);
                        if (allBiggerNewNrow)
                            if (j < newNrow)
                                allBiggerNewNrow = false;
                    }
                }

                if (allBiggerNewNrow) {
                    myList[myCols[i-1]].clear();
                    for (it = temp.begin(); it < temp.end(); it++) {
                        tLen = myList[myCols[*it]].size();
                        if (tLen > 0)
                            for (std::size_t j = 0; j < tLen; j++)
                                myList[myCols[i-1]].push_back(myList[myCols[*it]][j]);
                    }
                } else {
                    for (it = temp.begin(); it < temp.end(); it++)
                        myList[myCols[i-1]] = outersect(myList[myCols[i-1]], myList[myCols[*it]]);
                }
            }
        }
    }
    
    unsigned long int myLim, myCheck, lenFree = freeVariables.size();
    v1d posAns, ansVec, yExponents(mat[0].size(), 0), posVec(numCol, 0);
    ansVec.reserve(numCol);
    bool bExit = false;
    mpz_t mpzTemp1, mpzTemp2, mpzMin, xMpz, yMpz;

    if (lenFree > 0) {
        mpz_init (mpzTemp1);
        mpz_init (mpzTemp2);
        mpz_init (mpzMin);
        mpz_init (xMpz);
        mpz_init (yMpz);
        mpz_ui_pow_ui (mpzTemp1, 2, lenFree);
        mpz_sub_ui (mpzTemp1, mpzTemp1, 1);

        if (mpz_cmp_ui(mpzTemp1, hundredMillion) > 0)
            myLim = hundredMillion;
        else
            myLim = mpz_get_ui(mpzTemp1);

        for (std::size_t i = 1; i <= myLim; i++) {
            if (bExit) {break;}
            posAns = myIntToBit(i, lenFree);
            
            for (std::size_t j = 0; j < myList.size(); j++) {
                for (it = myList[j].begin(); it < myList[j].end(); it++) {
                    for (std::size_t k = 0; k < lenFree; k++) {
                        if (freeVariables[k] == *it) {
                            posVec[j] += posAns[k] % 2;
                            break;
                        }
                    }
                }
                if (posVec[j]==1)
                    ansVec.push_back((int64_t) j);
            }

            if (ansVec.size() > 0) {
                myCheck = 0;
                for (std::size_t j = 0; j < mat[0].size(); j++) {
                    for (it = ansVec.begin(); it < ansVec.end(); it++)
                        yExponents[j] += mat[*it][j];

                    myCheck += (yExponents[j] % 2);
                    yExponents[j] /= 2;
                }

                if (myCheck == 0) {
                    yExponents.erase(yExponents.begin());
                    mpz_set_ui(xMpz, 1);
                    mpz_set_ui(yMpz, 1);

                    for (it = ansVec.begin(); it < ansVec.end(); it++)
                         mpz_mul(xMpz, xMpz, test[*it]);

                    for (std::size_t j = 0; j < yExponents.size(); j++) {
                        mpz_ui_pow_ui(mpzTemp1, FB[j], yExponents[j]);
                        mpz_mul(yMpz, yMpz, mpzTemp1);
                    }

                    mpz_mod(xMpz, xMpz, n);
                    mpz_mod(yMpz, yMpz, n);
                    mpz_sub(mpzTemp1, xMpz, yMpz);
                    mpz_gcd(mpzTemp1, mpzTemp1, n);
                    mpz_add(mpzTemp2, xMpz, yMpz);
                    mpz_gcd(mpzTemp2, mpzTemp2, n);

                    if (mpz_cmp(mpzTemp1, mpzTemp2) < 0)
                        mpz_set(mpzMin, mpzTemp1);
                    else
                        mpz_set(mpzMin, mpzTemp2);

                    if (mpz_cmp_ui(mpzMin, 1) > 0) {
                        if (mpz_cmp(mpzTemp1, mpzTemp2) < 0) {
                            mpz_set(factors[0], mpzTemp1);
                            mpz_set(factors[1], mpzTemp2);
                        } else {
                            mpz_set(factors[1], mpzTemp1);
                            mpz_set(factors[0], mpzTemp2);
                        }
                        bExit = true;
                    }
                }
            }
            posVec = v1d(numCol, 0);
            yExponents = v1d(mat[0].size(), 0);
            ansVec.clear();
            ansVec.reserve(numCol);
        }
    }

    mpz_clear(mpzTemp1); mpz_clear(mpzTemp2);
    mpz_clear(mpzMin); mpz_clear(xMpz); mpz_clear(yMpz);
}

static v1d getPrimesQuadRes (mpz_t myN, double n) {
    std::vector<char> primes(n+1, 1);
    v1d myps;
    
    myps.reserve(std::floor((double)(2*n/std::log((double)n))));

    int k, ind, lastP = 3;
    int fsqr = (int) std::floor((double)std::sqrt((double)n));

    for (std::size_t j = 4; j <= n; j += 2)
        primes[j] = 0;

    while (lastP <= fsqr) {
        for (int j = lastP*lastP; j <= n; j += 2*lastP)
            primes[j] = 0;
        
        k = lastP + 2;
        ind = 2;
        while (!primes[k]) {
            k += 2;
            ind += 2;
        }
        lastP += ind;
    }

    mpz_t test, jmpz, temp;
    mpz_init(test);
    mpz_init(jmpz);
    mpz_init(temp);

    myps.push_back(2);

    for (int64_t j = 3; j <= n; j += 2) {
        if (primes[j]) {
            mpz_set_si(jmpz, j);
            mpz_set_si(temp, j);
            mpz_sub_ui(temp, jmpz, 1);
            mpz_div_2exp(temp, temp, 1);
            mpz_powm(test,myN,temp,jmpz);
            if (mpz_cmp_ui(test, 1) == 0)
                myps.push_back(j);
        }
    }

    mpz_clear(jmpz); mpz_clear(temp); mpz_clear(test);
    return myps;
}

static v3d sieveLists (int facLim, v1d FBase, int64_t vecLen, mpz_t sqrD[], v2d sieveD) {
    
    v3d outList(facLim, v2d(2, v1d()));
    uint64_t uiFB;
    unsigned long int uFacLim = (unsigned long int) facLim;
    int64_t siFB, vStrt1, vStrt2;
    mpz_t modTest;
    mpz_init(modTest);
    vStrt1 = vStrt2 = 0;
    
    for (std::size_t i = 1; i < uFacLim; i++) {
        uiFB = (uint64_t) FBase[i];
        siFB = FBase[i];
        for (int64_t j = 0; j < vecLen; j++) {
            mpz_mod_ui(modTest, sqrD[j], uiFB);
            if (mpz_cmp_ui(modTest, 0) == 0) {
                vStrt1 = j;
                break;
            }
        }
        
        for (int64_t j = vStrt1; j < vecLen; j += siFB)
            outList[i][0].push_back(j);

        for (int64_t j = vStrt1 + 1; j < vecLen; j++) {
            mpz_mod_ui(modTest, sqrD[j], uiFB);
            if (mpz_cmp_ui(modTest, 0) == 0) {
                vStrt2 = j;
                break;
            }
        }
        
        for (int64_t j = vStrt2; j < vecLen; j += siFB)
            outList[i][1].push_back(j);
    }
    
    mpz_clear(modTest);
    return outList;
}

void quadraticSieve (mpz_t myNum, double fudge1,
                     double fudge2, int64_t LenB, mpz_t factors[]) {
    
    unsigned long int digCount = mpz_sizeinbase(myNum, 10);
    unsigned long int bits = mpz_sizeinbase(myNum, 2);

    unsigned long int myTarget, facSize, facSize2;
    double sqrLogLog, LimB, lognum = bits / log2(std::exp(1.0));
    sqrLogLog = std::sqrt(lognum * std::log(lognum));
    mpz_t currP, nextP, resTest, CP1;
    mpz_init(currP);
    mpz_init(nextP);
    mpz_init(CP1);
    mpz_init(resTest);
    v1d facBase;
    double dblDigCount = (double) digCount;
    double dblMyTarget;

    // These values were obtained from "The Multiple Polynomial
    // Quadratic Sieve" by Robert D. Silverman
    // DigSize <- c(24,30,36,42,48,54,60,66)
    // FBSize <- c(100,200,400,900,1200,2000,3000,4500)
    // MSize <- c(5,20,35,100,125,250,350,500)
    //
    // rawCoef <- round(unname(lm(FBSize ~ poly(DigSize, 4, raw = TRUE))$coefficients), 7)
    // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
    // rawCoef
    //    intercept          x^1          x^2          x^3          x^4
    // 3637.0670996 -391.8275012   15.1541456   -0.2475566    0.0016806

    if (fudge1 == 0) {
        fudge1 = -0.4;
        LimB = std::exp((double) ((0.5 + fudge1) * sqrLogLog));
        
        dblMyTarget = -391.8275012 * dblDigCount + 15.1541456 * std::pow(dblDigCount, 2.0);
        dblMyTarget += -0.2475566*std::pow(dblDigCount, 3.0) + 0.0016806 * std::pow(dblDigCount, 4.0) + 3637.0671;
        dblMyTarget = std::ceil(dblMyTarget);
        myTarget = (unsigned long int) dblMyTarget;

        while (LimB < myTarget) {
            LimB = std::exp((double) ((0.5 + fudge1) * sqrLogLog));
            fudge1 += 0.001;
        }

        facBase = getPrimesQuadRes(myNum, LimB);
        facSize = facBase.size();

        while (facSize < myTarget) {
            fudge1 += 0.005;
            LimB = std::exp((double) ((0.5 + fudge1) * sqrLogLog));
            mpz_set_ui(currP, facBase.back());
            mpz_nextprime(nextP, currP);
            while (mpz_cmp_ui(nextP, LimB) < 0) {
                mpz_set(currP, nextP);
                mpz_nextprime(nextP, currP);
                mpz_sub_ui(CP1, currP, 1);
                mpz_div_2exp(CP1, CP1, 1);
                mpz_powm(resTest,myNum,CP1,currP);
                if (mpz_cmp_ui(resTest, 1) == 0)
                    facBase.push_back(mpz_get_si(currP));
            }
            facSize = facBase.size();
        }
    } else {
        LimB = std::exp((double)((.5+fudge1)*sqrLogLog));
        facBase = getPrimesQuadRes(myNum, LimB);
    }
    // rawCoef <- round(unname(lm(MSize ~ poly(DigSize, 4, raw = TRUE))$coefficients), 7)
    // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
    // rawCoef
    //    intercept           x^1           x^2           x^3           x^4
    // -213.1466450    23.3947361    -0.9494686     0.0165574    -0.0000767

    if (LenB == 0) {
        double dblLenB;
        dblLenB = 23.394736 * dblDigCount - 0.9494686 * std::pow(dblDigCount, 2.0);
        dblLenB += 0.0165574 * std::pow(dblDigCount, 3.0) - 0.0000767 * std::pow(dblDigCount, 4.0) - 213.1466450;
        dblLenB = std::ceil(dblLenB);
        LenB = (int64_t) dblLenB;
        LenB *= 1000;
    }
    
    facSize2 = facSize = facBase.size();

    // facBase2 will be used if multiple polynomials are needed.
    // With every iteration, a prime will be added to the factor
    // base for additional factorization. The original facBase
    // should not be tampered with, hence the need for facBase2.
    v1d facBase2 = facBase;

    mpz_t sqrtInt;
    mpz_init(sqrtInt);
    mpz_sqrt(sqrtInt, myNum);
    v1d myInterval;

    int64_t Lower = -1 * LenB, Upper = LenB, LenB2 = 2 * LenB + 1;
    uint64_t uLenB2 = (uint64_t) LenB2;
    myInterval.reserve(LenB2);
    for (int64_t i = Lower; i <= Upper; i++)
        myInterval.push_back(i);

    v2d SieveDist(facSize, v1d(2));
    SieveDist[0][0] = 1;
    SieveDist[0][1] = 1;

    // Getting quadratic residues. See tonellishanks.cc for more
    // details. The array "TS" was used here to make the code
    // more concise and since everything will be stored in
    // SieveDist, TS can easily be cleared from memory when done.
    mpz_t TS[13];
    for (std::size_t i = 0; i < 13; i++)
        mpz_init(TS[i]);

    int pow2, iter1, iter2;
    mpz_set_ui(TS[12], 2);

    for (std::size_t i = 1; i < facSize; i++) {
        mpz_set_ui(TS[5], facBase[i]);
        mpz_set_ui(TS[0], facBase[i]);
        mpz_sub_ui(TS[0], TS[0], 1);
        mpz_set(TS[1], TS[0]);
        mpz_set_ui(TS[7], 2);
        pow2 = (int) mpz_scan1 (TS[1], 0);
        mpz_div_2exp (TS[1], TS[1], pow2);

        if (pow2 == 1) {
            mpz_add_ui (TS[4], TS[5], 1);
            mpz_div_2exp (TS[4], TS[4], 2);
            mpz_powm (TS[2], myNum, TS[4], TS[5]);
            mpz_neg (TS[4], TS[2]);
            mpz_mod (TS[3], TS[4], TS[5]);
        } else {
            mpz_div_2exp (TS[4], TS[0], 1);
            mpz_powm (TS[6], TS[7], TS[4], TS[5]);
            while (mpz_cmp_ui(TS[6], 1) == 0) {
                mpz_add_ui(TS[7], TS[7], 1);
                mpz_powm (TS[6], TS[7], TS[4], TS[5]);
            }

            mpz_add_ui(TS[4], TS[1], 1);
            mpz_div_2exp(TS[4], TS[4], 1);
            mpz_powm(TS[10], myNum, TS[4], TS[5]);
            mpz_powm(TS[8], myNum, TS[1], TS[5]);
            mpz_powm(TS[9], TS[7], TS[1], TS[5]);

            iter1 = pow2;
            iter2 = 1;
            mpz_mod(TS[11], TS[8], TS[5]);

            while ((mpz_cmp_ui(TS[11], 1) != 0) && (iter2 != 0)) {
                iter2 = 0;
                mpz_mod(TS[11], TS[8], TS[5]);
                while (mpz_cmp_ui(TS[11], 1) != 0) {
                    iter2++;
                    mpz_pow_ui(TS[4], TS[12], iter2);
                    mpz_powm(TS[11], TS[8], TS[4], TS[5]);
                }
                if (iter2 != 0) {
                    mpz_pow_ui(TS[4], TS[12], iter1 - iter2 - 1);
                    mpz_powm(TS[4], TS[9], TS[4], TS[5]);
                    mpz_mul(TS[4], TS[4], TS[10]);
                    mpz_mod(TS[10], TS[4], TS[5]);

                    mpz_pow_ui(TS[4], TS[12], iter1 - iter2);
                    mpz_powm(TS[9], TS[9], TS[4], TS[5]);

                    mpz_mul(TS[4], TS[8], TS[9]);
                    mpz_mod(TS[8], TS[4], TS[5]);
                    iter1 = iter2;
                }
                mpz_set_ui(TS[11], 0);
            }
            mpz_set(TS[2], TS[10]);
            mpz_sub(TS[4], TS[5], TS[10]);
            mpz_mod(TS[3], TS[4], TS[5]);
        }
        SieveDist[i][0] = (int64_t) mpz_get_si(TS[2]);
        SieveDist[i][1] = (int64_t) mpz_get_si(TS[3]);
    }

    // Finished generating residues.. now free memory
    for (std::size_t i = 0; i < 13; i++)
        mpz_clear(TS[i]);

    mpz_t *largeInterval;
    mpz_t *sqrDiff;

    largeInterval = (mpz_t *) malloc(LenB2 * sizeof(mpz_t));
    sqrDiff = (mpz_t *) malloc(LenB2 * sizeof(mpz_t));

    for (std::size_t i = 0; i < uLenB2; i++) {
        mpz_init(largeInterval[i]);
        mpz_init(sqrDiff[i]);
    }

    mpz_sub_ui(largeInterval[0], sqrtInt, Upper);
    for (std::size_t i = 1; i < uLenB2; i++)
        mpz_add_ui(largeInterval[i], largeInterval[i-1], 1);

    mpz_t temp;
    mpz_init(temp);

    for (std::size_t i = 0; i < uLenB2; i++) {
        mpz_pow_ui(temp, largeInterval[i], 2);
        mpz_sub(sqrDiff[i], temp, myNum);
    }

    v3d FBDivSieve;
    FBDivSieve = sieveLists(facSize, facBase, (int64_t) LenB2, sqrDiff, SieveDist);

    mpz_mul_ui(temp, myNum, 2);
    mpz_sqrt(temp, temp);
    mpz_mul_ui(temp, temp, Upper);

    if (fudge2 == 0) {
        if (digCount < 30)
            fudge2 = 0.7;
        else if (digCount < 35)
            fudge2 = 0.9;
        else if (digCount < 40)
            fudge2 = 1.05;
        else if (digCount < 45)
            fudge2 = 1.1;
        else if (digCount < 50)
            fudge2 = 1.2;
        else
            fudge2 = 1.5;
    }

    double theCut = fudge2 * mpz_sizeinbase(temp, 10);
    mpz_t *mpzFacBase;
    mpzFacBase = (mpz_t *) malloc(facSize * sizeof(mpz_t));
    std::vector<double> LnFB, myLogs(LenB2, 0);
    LnFB.reserve(facSize);
    v1d::iterator it;

    std::size_t ind = 0;
    for (it = facBase.begin(); it < facBase.end(); it++, ind++) {
        mpz_init_set_ui(mpzFacBase[ind], *it);
        LnFB.push_back(std::log((double) *it));
    }

    unsigned long int tempSize, evenStrt;
    int64_t minPrime;
    mpz_mul_ui(temp, sqrtInt, Upper);
    minPrime = (int64_t) mpz_sizeinbase(temp, 10) * 2;
    v2d indexDiv(LenB2, v1d());

    if (mpz_even_p(sqrDiff[0]) != 0)
        evenStrt = 0;
    else
        evenStrt = 1;

    for (std::size_t j = evenStrt; j < uLenB2; j += 2)
        indexDiv[j].push_back(0);

    for (std::size_t i = 1; i < facSize; i++) {
        for (std::size_t k = 0; k <= 1; k++) {
            tempSize = FBDivSieve[i][k].size();
            for (std::size_t j = 0; j < tempSize; j++)
                indexDiv[FBDivSieve[i][k][j]].push_back((int64_t) i);
        }
        if (facBase[i] > minPrime) {
            tempSize = FBDivSieve[i][0].size();
            for (std::size_t j = 0; j < tempSize; j++)
                myLogs[FBDivSieve[i][0][j]] += LnFB[i];

            tempSize = FBDivSieve[i][1].size();
            for (std::size_t j = 0; j < tempSize; j++)
                myLogs[FBDivSieve[i][1][j]] += LnFB[i];
        }
    }

    v1d largeLogs;

    for (int64_t i = 0; i < LenB2; i++)
        if (myLogs[i] > theCut)
            largeLogs.push_back(i);

    unsigned long int largeLogsSize = largeLogs.size();
    v2d myMat(largeLogsSize, v1d(facSize+1, 0));

    for (std::size_t i = 0; i < largeLogsSize; i++) {
        if (mpz_sgn(sqrDiff[largeLogs[i]]) < 0) {
            myMat[i][0] = 1;
            mpz_abs(sqrDiff[largeLogs[i]], sqrDiff[largeLogs[i]]);
        }
    }

    mpz_t rem, quot;
    mpz_init(rem);
    mpz_init(quot);
    v1d sFacs;
    v1d myLargeLogs;
    bool divides = true;

    if (largeLogsSize > 0) {
        for (std::size_t j = 0; j < largeLogsSize; j++) {
            tempSize = indexDiv[largeLogs[j]].size();
            for (std::size_t i = 0; i < tempSize; i++) {
                while (divides) {
                    mpz_fdiv_qr_ui(quot, rem, sqrDiff[largeLogs[j]],
                                            facBase[indexDiv[largeLogs[j]][i]]);
                    divides = (mpz_cmp_ui(rem, 0) == 0);
                    if (divides) {
                        mpz_set(sqrDiff[largeLogs[j]], quot);
                        myMat[j][indexDiv[largeLogs[j]][i] + 1]++;
                    }
                }
                divides = true;
            }
            if (mpz_cmp_ui(sqrDiff[largeLogs[j]], 1) == 0) {
                // Found a smooth number
                sFacs.push_back((int64_t) j);
                myLargeLogs.push_back(largeLogs[j]);
            }
        }
    }

    unsigned long int lenM = sFacs.size();
    v3d listMatrix;
    v2d listLargeLogs;
    v2d tempMat(lenM, v1d(facSize + 1, 0));

    for (std::size_t i = 0; i < lenM; i++)
        for (std::size_t j = 0; j <= facSize; j++)
            tempMat[i][j] = myMat[sFacs[i]][j];

    listMatrix.push_back(tempMat);
    listLargeLogs.push_back(myLargeLogs);

    mpz_t A, B, B2, C, Atemp, Btemp, Atemp2;
    mpz_init(A); mpz_init(B); mpz_init(C); mpz_init(Atemp2);
    mpz_init(Atemp); mpz_init(Btemp); mpz_init(B2);

    mpz_mul_2exp(Atemp, myNum, 1);
    mpz_sqrt(Atemp, Atemp);
    mpz_div_ui(Atemp, Atemp, Upper);
    mpz_sqrt(Atemp, Atemp);

    int64_t maxFBase = *std::max_element(facBase.begin(), facBase.end());
    if (mpz_cmp_ui(Atemp, maxFBase) < 0)
        mpz_set_ui(Atemp, maxFBase);

    bool LegendreTest;
    mpz_t quadRes[2];
    mpz_init(quadRes[0]); mpz_init(quadRes[1]);
    v2d polySieveD;
    unsigned long int numPolys = 0, numSmooth = lenM;
    v1d myAtemps, myQuadRes;
    v1d myIntervalSqrd(LenB2);

    if (facSize2 > lenM)
        for (std::size_t i = 0; i < uLenB2; i++)
            myIntervalSqrd[i] = myInterval[i] * myInterval[i];
    

    // Find enough smooth numbers to guarantee a non-trivial solution
    int extraFacs = 0;
    while (mpz_cmp_ui(factors[0], 0) == 0) {
        while (numSmooth <= facSize + extraFacs) {
            LegendreTest = true;
            while (LegendreTest) {
                mpz_nextprime(Atemp, Atemp);
                mpz_sub_ui(temp, Atemp, 1);
                mpz_div_2exp(temp, temp, 1);
                mpz_powm(temp, myNum, temp, Atemp);
                if (mpz_cmp_ui(temp, 1) == 0)
                    LegendreTest = false;
            }

            myAtemps.push_back((int64_t) mpz_get_d(Atemp));
            facBase2.push_back((int64_t) mpz_get_d(Atemp));
            facSize2++;

            mpz_pow_ui(A, Atemp, 2);
            TonelliShanksC(myNum, Atemp, quadRes);

            if (mpz_cmp(quadRes[0], quadRes[1]) > 0)
                mpz_set(Btemp, quadRes[0]);
            else
                mpz_set(Btemp, quadRes[1]);

            myQuadRes.push_back((int64_t) mpz_get_d(Btemp));
            mpz_mul_2exp(temp, Btemp, 1);
            mpz_invert(temp, temp, Atemp);
            mpz_pow_ui(B2, Btemp, 2);
            mpz_sub(B2, myNum, B2);
            mpz_mul(B2, B2, temp);
            mpz_add(B2, B2, Btemp);
            mpz_mod(B2, B2, A);

            mpz_pow_ui(C, B2, 2);
            mpz_sub(C, C, myNum);
            mpz_divexact(C, C, A);

            numPolys++;
            polySieveD = v2d(facSize, v1d(2));

            for (std::size_t i = 0; i < facSize; i++) {
                mpz_invert(Atemp2, A, mpzFacBase[i]);
                for (std::size_t j = 0; j <= 1; j++) {
                    mpz_ui_sub(temp, SieveDist[i][j], B2);
                    mpz_mul(temp, temp, Atemp2);
                    mpz_mod_ui(temp, temp, facBase[i]);
                    polySieveD[i][j] = (int64_t) mpz_get_si(temp);
                }
            }

            for (std::size_t i = 0; i < uLenB2; i++) {
                mpz_mul_si(temp, B2, myInterval[i]);
                mpz_mul_2exp(temp, temp, 1);
                mpz_add(temp, temp, C);
                mpz_set_d(Atemp2, (double) myIntervalSqrd[i]);
                mpz_mul(Atemp2, Atemp2, A);
                mpz_add(sqrDiff[i], Atemp2, temp);
            }

            FBDivSieve = sieveLists(facSize, facBase, (int64_t) LenB2, sqrDiff, polySieveD);

            if (mpz_even_p(sqrDiff[0]) != 0)
                evenStrt = 0;
            else
                evenStrt = 1;

            indexDiv = v2d(LenB2, v1d());
            for (std::size_t j = evenStrt; j < uLenB2; j += 2)
                indexDiv[j].push_back(0);

            std::fill(myLogs.begin(), myLogs.end(), 0.0);

            for (std::size_t i = 1; i < facSize; i++) {
                for (std::size_t k = 0; k <= 1; k++) {
                    tempSize = FBDivSieve[i][k].size();
                    for (std::size_t j = 0; j < tempSize; j++)
                        indexDiv[FBDivSieve[i][k][j]].push_back((int64_t) i);
                }
                if (facBase[i] > minPrime) {
                    tempSize = FBDivSieve[i][0].size();
                    for (std::size_t j = 0; j < tempSize; j++)
                        myLogs[FBDivSieve[i][0][j]] += LnFB[i];

                    tempSize = FBDivSieve[i][1].size();
                    for (std::size_t j = 0; j < tempSize; j++)
                        myLogs[FBDivSieve[i][1][j]] += LnFB[i];
                }
            }

            largeLogs.clear();

            for (int64_t i = 0; i < LenB2; i++) 
                if (myLogs[i] > theCut)
                    largeLogs.push_back(i);

            largeLogsSize = largeLogs.size();
            myMat = v2d(largeLogsSize, v1d(facSize + 1, 0));

            for (std::size_t i = 0; i < largeLogsSize; i++) {
                if (mpz_sgn(sqrDiff[largeLogs[i]]) < 1) {
                    myMat[i][0] = 1;
                    mpz_abs(sqrDiff[largeLogs[i]], sqrDiff[largeLogs[i]]);
                }
            }

            sFacs.clear();
            myLargeLogs.clear();

            if (largeLogsSize > 0) {
                divides = true;

                for (std::size_t j = 0; j < largeLogsSize; j++) {
                    tempSize = indexDiv[largeLogs[j]].size();
                    for (std::size_t i = 0; i < tempSize; i++) {
                        while (divides) {
                            mpz_fdiv_qr_ui(quot, rem, sqrDiff[largeLogs[j]],
                                           facBase[indexDiv[largeLogs[j]][i]]);
                            divides = (mpz_cmp_ui(rem, 0) == 0);
                            if (divides) {
                                mpz_set(sqrDiff[largeLogs[j]], quot);
                                myMat[j][indexDiv[largeLogs[j]][i] + 1]++;
                            }
                        }
                        divides = true;
                    }
                    if (mpz_cmp_ui(sqrDiff[largeLogs[j]], 1) == 0) {
                        // Found a smooth number
                        sFacs.push_back((int64_t) j);
                        myLargeLogs.push_back(largeLogs[j]);
                    }
                }
            }

            lenM = sFacs.size();
            numSmooth += lenM;
            tempMat = v2d(lenM, v1d(facSize2 + 1, 0));

            for (std::size_t i = 0; i < lenM; i++) {
                for (std::size_t j = 0; j <= facSize; j++)
                    tempMat[i][j] = myMat[sFacs[i]][j];
                tempMat[i][facSize2] = 2;
            }

            listMatrix.push_back(tempMat);
            listLargeLogs.push_back(myLargeLogs);
        }

        mpz_t *newTestInt;
        newTestInt = (mpz_t *) malloc(numSmooth * sizeof(mpz_t));
        v2d newMat = v2d(numSmooth, v1d(facSize2 + 1, 0));
        unsigned long int row = 0, polyOne;

        for (std::size_t i = 0; i < listLargeLogs[0].size(); i++) {
            mpz_init_set(newTestInt[row], largeInterval[listLargeLogs[0][i]]);
            for (std::size_t j = 0; j <= facSize; j++)
                newMat[row][j] = listMatrix[0][i][j];
            row++;
        }

        unsigned long int fSize = facSize + 1;

        for (std::size_t k = 0; k < numPolys; k++) {
            polyOne = k + 1;
            fSize++;

            mpz_set_si(Atemp, myAtemps[k]);
            mpz_pow_ui(A, Atemp, 2);
            mpz_set_si(Btemp, myQuadRes[k]);

            mpz_mul_2exp(temp, Btemp, 1);
            mpz_invert(temp, temp, Atemp);
            mpz_pow_ui(B2, Btemp, 2);
            mpz_sub(B2, myNum, B2);
            mpz_mul(B2, B2, temp);
            mpz_add(B2, B2, Btemp);
            mpz_mod(B2, B2, A);

            for (std::size_t i = 0; i < listLargeLogs[polyOne].size(); i++) {
                mpz_mul_si(temp, A, myInterval[listLargeLogs[polyOne][i]]);
                mpz_init(newTestInt[row]);
                mpz_add(newTestInt[row], temp, B2);
                for (std::size_t j = 0; j < fSize; j++)
                    newMat[row][j] = listMatrix[polyOne][i][j];
                row++;
            }
        }

        solutionSearch (newMat, myNum, facBase2, newTestInt, factors);
        
        for (std::size_t i = 0; i < numSmooth; i++)
            mpz_clear(newTestInt[i]);

        extraFacs += 5;
    }

    for (std::size_t i = 0; i < uLenB2; i++) {
        mpz_clear(largeInterval[i]);
        mpz_clear(sqrDiff[i]);
    }

    mpz_clear(rem); mpz_clear(quot); mpz_clear(temp); mpz_clear(sqrtInt);
    mpz_clear(A); mpz_clear(B); mpz_clear(C); mpz_clear(Atemp);
    mpz_clear(quadRes[0]); mpz_clear(quadRes[1]);
}
