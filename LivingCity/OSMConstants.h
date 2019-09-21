#ifndef LIVING_CITY_OSM_CONSTANTS_H_
#define LIVING_CITY_OSM_CONSTANTS_H_

#include <string>

enum class OSMConstant
{
    Empty,
    MotorwayJunction,
    TrafficSignals,
    StopJunction,
    TurningCircle
};

OSMConstant mapStringToOSMConstant(const std::string stringId);


#endif  // LIVING_CITY_OSM_CONSTANTS_H_
