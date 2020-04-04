#ifndef STATS_UTILS_H
#define STATS_UTILS_H

#include <RcppThread.h>
#include <chrono>

// Used for checking whether user has interrupted computation
constexpr auto checkInterTime = std::chrono::milliseconds(1000);
constexpr auto milliCutOff = std::chrono::milliseconds(1000);
constexpr auto secondCutOff = std::chrono::seconds(60);
constexpr auto minuteCutOff = std::chrono::minutes(60);

constexpr std::size_t maxHours = 999 * 60;
constexpr std::size_t maxMinutes = (maxHours * 60) + 59;

void makeStrLen(std::string & myStr, std::size_t myLen) {
    
    while (myStr.size() < myLen) {
        myStr = " " + myStr + " ";
    }
    
    myStr.resize(myLen);
}

template <typename typeTime>
void makeStats(std::size_t currLim, std::size_t loopLimit, std::size_t nPolys,
               std::size_t nSmooth, std::size_t nPartial, std::size_t &lastLim,
               typeTime timeDiff) {
    
    std::string strPerc = std::to_string(100 * currLim / loopLimit) + "%";
    std::string strPolys = std::to_string(nPolys);
    std::string strSmooth = std::to_string(nSmooth);
    std::string strPartial = std::to_string(nPartial);
    
    makeStrLen(strPerc, 10);
    makeStrLen(strPolys, 13);
    makeStrLen(strSmooth, 12);
    makeStrLen(strPartial, 12);
    
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
    makeStrLen(myTime, 20);
    
    RcppThread::Rcout << "\r|" << myTime << "|" << strPerc << "|" << strPolys << "|"
                      << strSmooth << "|" << strPartial << "|";
    lastLim = currLim;
}

template <typename typeTime>
void updateStatTime(std::size_t n, std::size_t facSize,
                    typeTime timeDiff, typeTime &showStatsTime) {
    if (n > 0) {
        const std::size_t percentComplete = ((100 * n) / facSize) + 1;
        const auto onePercentTime = timeDiff / percentComplete;
        
        if (onePercentTime > std::chrono::seconds(15)) {
            showStatsTime = std::chrono::seconds(15);
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

#endif
