#ifndef SIMULATOR_CONFIGURATION_H_
#define SIMULATOR_CONFIGURATION_H_

#include <QString>
#include <qcoreapplication.h>
#include <QSettings>

#include "./routing.h"

namespace LC
{

class SimulatorConfiguration
{

public:
    SimulatorConfiguration(const char * settings_path);

    bool TryToReadPreviousPaths(void) const;
    float DeltaTime(void) const;
    float SimulationStartingHour(void) const;
    float SimulationEndingHour(void) const;
    Routing SimulationRouting(void) const;
    int AmountOfPasses(void) const;
    int LimitOfPeople(void) const;
    QString NetworkPath(void) const;
    bool UseCPU(void) const;
    bool AddRandomPeople(void) const;

    bool ShouldExportTimes(void) const;
    bool ShouldExportPeopleSummary(void) const;
    
private:
    QSettings settings_;
};

}
#endif  // SIMULATOR_CONFIGURATION_H_
