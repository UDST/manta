#ifndef BENCHMARKER_H_
#define BENCHMARKER_H_

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

using Timestamp = std::chrono::time_point<std::chrono::high_resolution_clock>;
using Duration = std::chrono::nanoseconds;


class Benchmarker {
    public:
        Benchmarker(const std::string desc, int depth);

        void startMeasuring();
        void stopMeasuring();
        void endBenchmark();
        void stopAndEndBenchmark();

    private:
        bool on;
        Timestamp lastTimestamp;
        Duration elapsed;
        std::string description;
        std::string margin;

        static std::ofstream outStream;
        static int amountOpened;
};


#endif  // BENCHMARKER_H_

