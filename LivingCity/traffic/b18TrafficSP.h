/************************************************************************************************
*		@desc Class that finds the path for each person using Johnsons
*		@author igaciad
************************************************************************************************/
#ifndef LC_B18_TRAFFIC_SP_H
#define LC_B18_TRAFFIC_SP_H

#include "b18TrafficPerson.h"
#include "../RoadGraph/roadGraph.h"
#include "sp/graph.h"
#include "sp/csv.h"
#include "omp.h"

#include "sp/config.h"
#include "sp/mpi_wrapper.h"

namespace LC {
class B18TrafficSP {

 public:

  static void generateRoutesSP(
      LC::RoadGraph::roadBGLGraph_BI &roadGraph,
      std::vector<B18TrafficPerson> &trafficPersonVec,
      std::vector<uint>& indexPathVec,
      std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum,
      int weigthMode = 0,
      float sample = 1.0f);
};
}

#endif  // LC_B18_TRAFFIC_SP_H
