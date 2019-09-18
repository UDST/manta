#include "./simulatorConfiguration.h"

namespace LC
{

SimulatorConfiguration::SimulatorConfiguration(const char * settings_path) :
    settings_(settings_path, QSettings::IniFormat) {}

bool SimulatorConfiguration::TryToReadPreviousPaths(void) const
{
    return settings_.value("USE_PREV_PATHS", false).toBool();
}

float SimulatorConfiguration::DeltaTime(void) const
{
    return settings_.value("TIME_STEP", .5).toFloat();
}


float SimulatorConfiguration::SimulationStartingHour(void) const
{
    return settings_.value("START_HR", 5).toFloat();
}


float SimulatorConfiguration::SimulationEndingHour(void) const
{
    return settings_.value("END_HR", 12).toFloat();
}

Routing SimulatorConfiguration::SimulationRouting(void) const
{
    const bool useJohnsonRouting = settings_.value("USE_JOHNSON_ROUTING", false).toBool();
    const bool useSPRouting = settings_.value("USE_SP_ROUTING", false).toBool();
    const Routing routing =
        useSPRouting ? Routing::SP
        : useJohnsonRouting ? Routing::Johnson
        : Routing::Dijkstra;
    return routing;
}

int SimulatorConfiguration::AmountOfPasses(void) const
{
    return settings_.value("NUM_PASSES", 1).toInt();
}

int SimulatorConfiguration::LimitOfPeople(void) const
{
    return settings_.value("LIMIT_NUM_PEOPLE", -1).toInt();
}
    
QString SimulatorConfiguration::NetworkPath(void) const
{
    return settings_.value("NETWORK_PATH").toString();
}

bool SimulatorConfiguration::UseCPU(void) const
{
    // The CPU implementation has not been updated with the latest intersection changes.
    return false;
}

bool SimulatorConfiguration::AddRandomPeople(void) const
{
    return settings_.value("ADD_RANDOM_PEOPLE").toBool();
}


}
