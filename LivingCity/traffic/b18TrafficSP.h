/************************************************************************************************
*		@desc Class that finds the path for each person using Johnsons
*		@author igaciad
************************************************************************************************/
#ifndef LC_B18_TRAFFIC_SP_H
#define LC_B18_TRAFFIC_SP_H

#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <vector>

#include "b18TrafficPerson.h"
#include "../RoadGraph/roadGraph.h"
#include "sp/graph.h"
#include "sp/external/csv.h"
#include "omp.h"

#include "sp/config.h"
#include "sp/mpi_wrapper.h"

namespace LC {
class B18TrafficSP {

 public:
	 
  static std::vector<abm::graph::vertex_t> compute_routes(int mpi_rank,
                                                   int mpi_size,
                                                   std::shared_ptr<abm::Graph> graph_,
                                                   std::vector<std::array<abm::graph::vertex_t, 2>> all_od_pairs_);

  static std::vector<std::array<abm::graph::vertex_t, 2>> make_od_pairs(std::vector<B18TrafficPerson> trafficPersonVec,
                                                                               int nagents);

  static void generateRoutesSP(
      LC::RoadGraph::roadBGLGraph_BI &roadGraph,
      std::vector<B18TrafficPerson> &trafficPersonVec,
      std::vector<uint>& indexPathVec,
      std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum,
      int weigthMode = 0,
      float sample = 1.0f);

};
} //namespace LC

#endif  // LC_B18_TRAFFIC_SP_H
