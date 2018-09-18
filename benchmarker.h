#ifndef BENCHMARKER__H
#define BENCHMARKER__H

#include <string>
#include <chrono>
#include <iostream>
#include <fstream>

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


extern Benchmarker mainBench;
extern Benchmarker intersectionBench;
extern Benchmarker peopleBench;


#endif  // BENCHMARKER__H

