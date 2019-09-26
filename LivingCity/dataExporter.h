#ifndef LIVINGCITY_DATAEXPORTER_H_
#define LIVINGCITY_DATAEXPORTER_H_

#include <assert.h>

#include <map>
#include <string>
#include <vector>
#include <QTime>


#include "./traffic/b18TrafficPerson.h"

enum Phase {Initialization, Routing, Simulation, Export};

class DataExporter
{

public:
    DataExporter(void);

    void SwitchMeasuringFromTo(const Phase & current_phase, const Phase & new_phase);

    void ExportTimes(void);

    void ExportPersonsSummary(
        const std::vector<LC::B18TrafficPerson> & traffic_persons,
        const std::vector<float> & persons_travelled_distances) const;

private:
    Phase current_phase_;
    QTime timer_;

    // Stores the amount of milliseconds of each phase
    std::map<Phase, int> phases_total_times_;

    std::string TagForPhase(const Phase & phase) const;
};


#endif  // LIVINGCITY_DATAEXPORTER_H_
