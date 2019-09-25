#include <iostream>

#include "./dataExporter.h"

DataExporter::DataExporter(void) : current_phase_(Phase::Initialization)
{
    for (int phase_int = Phase::Initialization; phase_int <= Phase::Export; phase_int++) {
        Phase phase = static_cast<Phase>(phase_int);
        phases_total_times_[phase] = 0;
    }

    timer_.start();

    current_phase_ = Phase::Initialization;
}


void DataExporter::SwitchMeasuringFromTo(const Phase & current_phase, const Phase & new_phase)
{
    assert(current_phase_ == current_phase);

    phases_total_times_.at(current_phase) += timer_.restart();
    current_phase_ = new_phase;
}

void DataExporter::ExportTimes(void)
{
    assert(current_phase_ == Phase::Export);

    phases_total_times_.at(current_phase_) += timer_.restart();
    for (int phase_int = Phase::Initialization; phase_int <= Phase::Export; phase_int++) {
        Phase phase = static_cast<Phase>(phase_int);
        std::cout
            << "[Time] " << TagForPhase(phase) << " took "
            << phases_total_times_.at(phase) << "ms."
            << std::endl;
    }

}

std::string DataExporter::TagForPhase(const Phase & phase) const
{
    if (phase == Phase::Initialization)
        return "initialization";
    if (phase == Phase::Routing)
        return "routing";
    if (phase == Phase::Simulation)
        return "simulation";

    return "export";
}
