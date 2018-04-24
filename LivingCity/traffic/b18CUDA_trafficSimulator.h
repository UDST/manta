/************************************************************************************************
 *
 *		CUDA hearder
 *
 *		@author igarciad
 *
 ************************************************************************************************/

#ifndef B18_TRAFFIC_SIMULATION_H
#define B18_TRAFFIC_SIMULATION_H

#include "b18TrafficPerson.h"
#include "b18EdgeData.h"
#include <vector>


extern void b18InitCUDA(std::vector<LC::B18TrafficPerson> &trafficPersonVec,
                        std::vector<uint> &indexPathVec,
                        std::vector<LC::B18EdgeData> &edgesData, std::vector<uchar> &laneMap,
                        std::vector<uchar> &trafficLights,
                        std::vector<LC::B18IntersectionData> &intersections,
                        float startTimeH, float endTimeH,
                        std::vector<float>& accSpeedPerLinePerTimeInterval,
                        std::vector<float>& numVehPerLinePerTimeInterval);
extern void b18GetDataCUDA(std::vector<LC::B18TrafficPerson> &trafficPersonVec,
                           std::vector<uchar> &trafficLights);
extern void b18GetSampleTraffic(std::vector<float>& accSpeedPerLinePerTimeInterval, 
                                std::vector<float>& numVehPerLinePerTimeInterval);
extern void b18FinishCUDA(void);
extern void b18SimulateTrafficCUDA(float currentTime, uint numPeople,
                                   uint numIntersections);

#endif // B18_TRAFFIC_SIMULATION_H

