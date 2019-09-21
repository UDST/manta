#include "OSMConstants.h"

OSMConstant mapStringToOSMConstant(const std::string stringId)
{
    if (stringId == "motorway_junction")
        return OSMConstant::MotorwayJunction;

    if (stringId == "traffic_signals")
        return OSMConstant::TrafficSignals;

    if (stringId == "stop")
        return OSMConstant::StopJunction;

    if (stringId == "turning_circle")
        return OSMConstant::TurningCircle;

    return OSMConstant::Empty;
}
