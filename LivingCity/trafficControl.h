#ifndef LIVING_CITY_TRAFFIC_CONTROL_H_
#define LIVING_CITY_TRAFFIC_CONTROL_H_

#include "./OSMConstants.h"


enum class TrafficControl
{
    TrafficLight,
    Unsupervised
};

TrafficControl mapOSMConstantToTrafficControl(const OSMConstant & osmConstant);


#endif  // LIVING_CITY_TRAFFIC_CONTROL_H_
