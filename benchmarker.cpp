#include "benchmarker.h"


std::ofstream Benchmarker::outStream;
int Benchmarker::amountOpened = 0;

Benchmarker::Benchmarker(const std::string desc, int depth) :
    on(false),
    description(desc),
    margin(depth * 4, ' ')
{
    if (amountOpened == 0) outStream.open("timestamps.info");

    ++amountOpened;
}

void Benchmarker::begin()
{
    if (on) return;
    on = true;

    beginning = std::chrono::steady_clock::now();
    outStream << margin << " >> Started: " << description << std::endl;

    return;
}

void Benchmarker::end()
{
    if (!on) return;
    Timestamp end = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - beginning);
    outStream << margin
        << " << Ended: " << description
        << " (elapsed time: " << elapsedTime.count() <<")"
        << std::endl;

    --amountOpened;
    if (amountOpened == 0) outStream.close();
    on = false;
}

