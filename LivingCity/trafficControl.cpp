#include "./trafficControl.h"

#include "./OSMConstants.h"


TrafficControl mapOSMConstantToTrafficControl(const OSMConstant & osmConstant)
{
    return osmConstant == OSMConstant::TrafficSignals
        ? TrafficControl::TrafficLight
        : TrafficControl::Unsupervised;
}
