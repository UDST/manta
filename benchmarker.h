#ifndef BENCHMARKER__H
#define BENCHMARKER__H

#include <string>
#include <chrono>
#include <iostream>
#include <fstream>

using Timestamp = std::chrono::time_point<std::chrono::steady_clock>;


class Benchmarker {
    public:
        Benchmarker(const std::string desc, int depth);

        void begin();
        void end();

    private:
        bool on;
        Timestamp beginning;
        std::string description;
        std::string margin;

        static std::ofstream outStream;
        static int amountOpened;
};


#endif  // BENCHMARKER__H

