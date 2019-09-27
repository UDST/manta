#include <iostream>
#include <fstream>

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

void DataExporter::ExportPersonsSummary(
    const std::vector<LC::B18TrafficPerson> & traffic_persons,
    const std::vector<float> & persons_travelled_distances) const
{
    assert(current_phase_ == Phase::Export);
    assert(traffic_persons.size() == persons_travelled_distances.size());

    std::ofstream peoples_file("people.csv");

    peoples_file
        << "p"
        << ",init_intersection"
        << ",end_intersection"
        << ",time_departure"
        << ",num_steps"
        << ",co"
        << ",gas"
        << ",distance"
        << ",a"
        << ",b"
        << ",T"
        << ",avg_v(mph)"
        << std::endl;

    for (size_t p = 0; p < traffic_persons.size(); ++p) {
        peoples_file
            << p
            << "," << traffic_persons.at(p).init_intersection
            << "," << traffic_persons.at(p).end_intersection
            << "," << traffic_persons.at(p).time_departure
            << "," << traffic_persons.at(p).num_steps
            << "," << traffic_persons.at(p).co
            << "," << traffic_persons.at(p).gas
            << "," << persons_travelled_distances.at(p)
            << "," << traffic_persons.at(p).a
            << "," << traffic_persons.at(p).b
            << "," << traffic_persons.at(p).T
            << "," << (traffic_persons[p].cumulative_velocity / traffic_persons[p].num_steps) * 3600 / 1609.34
            << std::endl;
    }

    peoples_file.close();
}

void DataExporter::ExportRoutes(const std::vector<std::vector<uint>> & persons_routes) const
{
    std::ofstream persons_routes_file("routes.csv");

    persons_routes_file << "p:route" << std::endl;

    for (size_t p = 0; p < persons_routes.size(); ++p) {
        persons_routes_file << p << ":[";

        const auto & route = persons_routes.at(p);
        if (!route.empty()) {
            size_t route_index = 0;
            persons_routes_file << route.at(route_index);
            route_index++;
            while (route_index < route.size()) {
                persons_routes_file << "," << route.at(route_index);
                route_index++;
            }
        }

        persons_routes_file << "]" << std::endl;
    }

    persons_routes_file.close();
}

std::string DataExporter::TagForPhase(const Phase & phase) const
{
    if (phase == Phase::Initialization)
        return "Initialization";
    if (phase == Phase::Routing)
        return "Routing";
    if (phase == Phase::Simulation)
        return "Simulation";
    if (phase == Phase::Export)
        return "Export";

    throw std::runtime_error("DataExporter::TagForPhase -> Unsupported tag.");
}
