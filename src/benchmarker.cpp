#include "benchmarker.h"


std::ofstream Benchmarker::outStream;
int Benchmarker::amountOpened = 0;

Benchmarker::Benchmarker(const std::string desc, int depth) :
    on(false),
    elapsed(Duration::zero()),
    description(desc),
    margin(depth * 4, ' ')
{
    if (amountOpened == 0) outStream.open("timestamps.info");
    ++amountOpened;
}

void Benchmarker::startMeasuring()
{
    if (on) return;
    on = true;

    lastTimestamp = std::chrono::high_resolution_clock::now(); 

    return;
}

void Benchmarker::stopMeasuring()
{
    if (!on) return;

    elapsed += std::chrono::duration_cast<Duration>(
            std::chrono::high_resolution_clock::now() - lastTimestamp);

    on = false;
}

void Benchmarker::endBenchmark()
{
    outStream << margin
        << " << Ended: " << description
        << " (elapsed time: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << ")"
        << std::endl;

    --amountOpened;
    if (amountOpened == 0) outStream.close();
}

void Benchmarker::stopAndEndBenchmark()
{
    stopMeasuring();
    endBenchmark();
}

Benchmarker mainBench("main function", 0);
Benchmarker intersectionBench("Intersections total time", 0);
Benchmarker peopleBench("People total time", 0);
