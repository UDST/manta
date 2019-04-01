/************************************************************************************************
*		@desc Class that finds the path for each person using Johnsons
*		@author igaciad
************************************************************************************************/
#ifndef LC_B18_TRAFFIC_SP_H
#define LC_B18_TRAFFIC_SP_H

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "src/benchmarker.h"
#include "src/linux_host_memory_logger.h"
#include "b18TrafficPerson.h"
#include "../RoadGraph/roadGraph.h"
#include "sp/graph.h"
//#include "sp/external/csv.h"
#include "omp.h"

#include "sp/config.h"
#include "sp/mpi_wrapper.h"

namespace LC {
class B18TrafficSP {

 public:
	 
	 /*
  static std::vector<abm::graph::vertex_t> compute_routes(int mpi_rank,
                                                   	  int mpi_size,
                                                   	  std::shared_ptr<abm::Graph> graph_,
                                                   	  std::vector<std::array<abm::graph::vertex_t, 2>> all_od_pairs_);
*/
  static std::vector<abm::graph::vertex_t> compute_routes(int mpi_rank,
                                                          int mpi_size,
                                                          const std::shared_ptr<abm::Graph>& graph_,
                                                          const std::vector<std::array<abm::graph::vertex_t, 2>>& all_od_pairs_);
  
  static std::vector<std::array<abm::graph::vertex_t, 2>> make_od_pairs(std::vector<B18TrafficPerson> trafficPersonVec,
                                                                        int nagents);

 static std::vector<std::array<abm::graph::vertex_t, 2>> read_od_pairs(const std::string& filename, int nagents);
 // static abm::graph::Graph read_graph_osm(const std::string& filename);

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
