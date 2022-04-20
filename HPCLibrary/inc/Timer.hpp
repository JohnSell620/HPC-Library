/*
* Timer.hpp
* Description: Timer class for benchmarking runs.
* @author Johnny Sellers
* @version 0.1 04/28/2017
*/
#if !defined(TIMER_HPP)
#define TIMER_HPP

#include <chrono>

class Timer {
public:
    Timer() : startTime(), stopTime() {}

    std::chrono::system_clock::time_point start() {
        return (startTime = std::chrono::system_clock::now());
    }
    
    std::chrono::system_clock::time_point stop() {
        return (stopTime  = std::chrono::system_clock::now());
    }

    double elapsed() {
        return std::chrono::duration_cast<std::chrono::milliseconds>(stopTime-startTime).count();
    }

private:
    std::chrono::time_point<std::chrono::system_clock> startTime, stopTime;
};

#endif // TIMER_HPP
