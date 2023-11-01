//
// Created by Kai Zhao on 10/30/20.
//

#ifndef SZ_TIMER_HPP
#define SZ_TIMER_HPP

#include <string>
#include <iostream>
#include <chrono>

namespace SZ3 {
    class Timer {
    public:
        Timer() = default;

        Timer(bool initstart) {
            if (initstart) {
                start();
            }
        }

        void start() {
            begin = std::chrono::steady_clock::now();
        }

        double stop() {
            end = std::chrono::steady_clock::now();
            return std::chrono::duration<double>(end - begin).count();
        }

        double stop(const std::string &msg) {
            double seconds = stop();
#if SZ3_DEBUG_TIMINGS
            std::cout << msg << " time = " << seconds << "s" << std::endl;
#endif
            return seconds;
        }

    private:
        std::chrono::time_point<std::chrono::steady_clock> begin, end;
    };
};


#endif //SZ_TIMER_HPP
