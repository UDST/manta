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
	 
  static std::vector<abm::graph::vertex_t> compute_routes(int mpi_rank,
                                                          int mpi_size,
                                                          const std::shared_ptr<abm::Graph>& graph_,
                                                          const std::vector<std::array<abm::graph::vertex_t, 2>>& od_pairs);

  static std::vector<std::array<abm::graph::vertex_t, 2>> make_od_pairs(std::vector<B18TrafficPerson> trafficPersonVec,
                                                                        int nagents);

  static std::vector<std::array<abm::graph::vertex_t, 2>> read_od_pairs(const std::string& filename, int nagents);
  
  static std::vector<float> read_dep_times(const std::string& filename);

  static std::vector<abm::graph::edge_id_t> RoutingWrapper(
  const std::vector<std::array<abm::graph::vertex_t, 2>> all_od_pairs_,
  const std::shared_ptr<abm::Graph>& street_graph,
  const std::vector<float>& dep_times,
  const float start_time_mins,
  const float end_time_mins);

  static void initialize_person_to_init_edge(
    std::vector<abm::graph::edge_id_t>& all_paths,
    const std::shared_ptr<abm::Graph>& street_graph);

  static void filterODByTimeRange(
    const std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs,
    const std::vector<float> dep_times,
    const float start_time_mins,
    const float end_time_mins,
    std::vector<abm::graph::vertex_t>& filtered_od_pairs_sources_,
    std::vector<abm::graph::vertex_t>& filtered_od_pairs_targets_,
    std::vector<float>& filtered_dep_times_);

  static void convertVector(std::vector<abm::graph::edge_id_t> paths_SP, std::vector<uint>& indexPathVec, std::vector<uint> &edgeIdToLaneMapNum, const std::shared_ptr<abm::Graph>& graph_);

  explicit B18TrafficSP(const std::shared_ptr<abm::Graph>& graph) : graph_{graph} {};
 private:
  //all od pairs
  std::shared_ptr<std::vector<std::array<abm::graph::vertex_t, 2>>> all_od_pairs_;
  
  //filtered od pairs
  std::shared_ptr<std::vector<std::array<abm::graph::vertex_t, 2>>> filtered_od_pairs_;

  //graph (street network)
  std::shared_ptr<abm::Graph> graph_;

  //all paths
  std::vector<abm::graph::vertex_t> all_paths_;
};
} //namespace LC

#endif  // LC_B18_TRAFFIC_SP_H
