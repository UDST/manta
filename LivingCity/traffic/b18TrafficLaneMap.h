/************************************************************************************************
*
*		LC Project - B18 Traffic lane map
*
*
*		@desc Class that contains the structure of the lane maps
*		@author igaciad
*
************************************************************************************************/

#ifndef LC_B18_TRAFFIC_LANEMAP_H
#define LC_B18_TRAFFIC_LANEMAP_H


#include "../misctools/misctools.h"
#include "RoadGraph/roadGraph.h"
#include "b18EdgeData.h"

namespace LC {


class B18TrafficLaneMap {
 public:
   B18TrafficLaneMap();
   ~B18TrafficLaneMap();

  void createLaneMap(RoadGraph &inRoadGraph, std::vector<uchar> &laneMap,
                     std::vector<B18EdgeData> &edgesData, std::vector<B18IntersectionData> &intersections,
                     std::vector<uchar> &trafficLights,
                     std::map<uint, RoadGraph::roadGraphEdgeDesc_BI> &laneMapNumToEdgeDesc,
                     std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum);

  void resetIntersections(std::vector<B18IntersectionData> &intersections,
                          std::vector<uchar> &trafficLights);

  
};
}

#endif  // LC_B18_TRAFFIC_LANEMAP_H
