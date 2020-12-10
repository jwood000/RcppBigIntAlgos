#ifndef STATS_UTILS_H
#define STATS_UTILS_H

#include <RcppThread.h>
#include <chrono>

// Used for checking whether user has interrupted computation
constexpr auto checkInterTime = std::chrono::milliseconds(1000);
constexpr auto milliCutOff = std::chrono::milliseconds(1000);
constexpr auto secondCutOff = std::chrono::seconds(60);
constexpr auto minuteCutOff = std::chrono::minutes(60);

constexpr auto fifteenSeconds = std::chrono::seconds(15);
constexpr std::size_t maxHours = 999 * 60;
constexpr std::size_t maxMinutes = (maxHours * 60) + 59;

constexpr std::size_t colOneWidth = 20u;
using typeTimePoint = std::chrono::time_point<std::chrono::steady_clock>;

static void MakeStrLen(std::string & myStr, std::size_t myLen) {
    
    while (myStr.size() < myLen) {
        myStr = " " + myStr + " ";
    }
    
    myStr.resize(myLen);
}

template <typename typeTime>
std::string GetTime(typeTime timeDiff) {
    std::string myTime;
    
    std::size_t nMilliSec = std::chrono::duration_cast<std::chrono::milliseconds>(timeDiff).count();
    std::size_t nSeconds = std::chrono::duration_cast<std::chrono::seconds>(timeDiff).count();
    std::size_t nMinutes = std::chrono::duration_cast<std::chrono::minutes>(timeDiff).count();
    
    // Include hours
    if ((timeDiff) > minuteCutOff) {
        const std::size_t nHours = (nMinutes > maxMinutes) ? maxHours : nMinutes / 60;
        myTime = std::to_string(nHours) + "h ";
        nMinutes -= nHours * 60;
        nSeconds -= nHours * 60 * 60;
        nMilliSec -= nHours * 60 * 60 * 1000;
    }
    
    // Include minutes
    if ((timeDiff)  > secondCutOff) {
        myTime += std::to_string(nMinutes) + "m ";
        nSeconds -= nMinutes * 60;
        nMilliSec -= nMinutes * 60 * 1000;
    }
    
    // Include seconds
    if ((timeDiff)  > milliCutOff) {
        myTime += std::to_string(nSeconds) + "s ";
        nMilliSec -= nSeconds * 1000;
    }
    
    myTime += std::to_string(nMilliSec) + "ms";
    return myTime;
}

template <typename typeTime>
void MakeStats(std::size_t loopLimit, std::size_t nPolys,
               std::size_t nSmooth, std::size_t nPartial, typeTime timeDiff) {
    
    std::string strPerc = std::to_string(100 * (nSmooth + nPartial) / loopLimit) + "%";
    std::string strPolys = std::to_string(nPolys);
    std::string strSmooth = std::to_string(nSmooth);
    std::string strPartial = std::to_string(nPartial);
    
    MakeStrLen(strPerc, 10);
    MakeStrLen(strPolys, 13);
    MakeStrLen(strSmooth, 12);
    MakeStrLen(strPartial, 12);
    
    std::string myTime = GetTime(timeDiff);
    MakeStrLen(myTime, colOneWidth);
    
    RcppThread::Rcout << "\r|" << myTime << "|" << strPerc << "|" << strPolys << "|"
                      << strSmooth << "|" << strPartial << "|";
}

template <typename typeTime>
void UpdateStatTime(std::size_t n, std::size_t facSize,
                    typeTime timeDiff, typeTime &showStatsTime) {
    if (n > 0) {
        const std::size_t percentComplete = ((100 * n) / facSize) + 1;
        const auto onePercentTime = timeDiff / percentComplete;
        
        if (onePercentTime > fifteenSeconds) {
            showStatsTime = fifteenSeconds;
        } else {
            if (onePercentTime < std::chrono::seconds(1)) {
                showStatsTime = 5 * onePercentTime;
            } else {
                showStatsTime = onePercentTime;
            }
        }
    } else {
        showStatsTime = std::chrono::milliseconds(500);
    }
}

template <typename typeTime>
void OneColumnStats(typeTime timeDiff) {
    std::string myTime = GetTime(timeDiff);
    MakeStrLen(myTime, colOneWidth);
    RcppThread::Rcout << "\r|" << myTime << "|";
}

template <typename typeTime>
void TwoColumnStats(typeTime timeDiff, std::size_t valOne,
                    std::size_t valTwo, bool matrix = true) {
    
    std::string myTime = GetTime(timeDiff);
    MakeStrLen(myTime, colOneWidth);
    std::string myDim = (matrix) ? std::to_string(valOne) + " x "
                                   + std::to_string(valTwo) : std::to_string(valOne); 
    MakeStrLen(myDim, colOneWidth);
    RcppThread::Rcout << "\r|" << myTime << "|" << myDim << "|";
}

#endif
