#include "b18TrafficSP.h"

#include <boost/graph/exterior_property.hpp>
#include "src/linux_host_memory_logger.h"
#include "roadGraphB2018Loader.h"

#define ROUTE_DEBUG 0
//#define DEBUG_JOHNSON 0

namespace LC {


////////////////
/////////////////////////////
using namespace boost;
using namespace std::chrono;

inline bool fileExists(const std::string& fileName) {
  std::ifstream f(fileName.c_str());
  return f.good();
}

typedef exterior_vertex_property<RoadGraph::roadBGLGraph_BI, float>
DistanceProperty;
typedef DistanceProperty::matrix_type DistanceMatrix;
typedef DistanceProperty::matrix_map_type DistanceMatrixMap;

bool myfunction (std::array<abm::graph::vertex_t, 2> i, std::array<abm::graph::vertex_t, 2> j) {
  return (i==j);
}

// Convert OD pairs to SP graph format
std::vector<std::array<abm::graph::vertex_t, 2>> B18TrafficSP::make_od_pairs(std::vector<B18TrafficPerson> trafficPersonVec, 
									     int nagents) {
  bool status = true;
  std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs;
  try {
    abm::graph::vertex_t v1, v2;
    abm::graph::weight_t weight;
    printf("trafficPersonSize = %d\n", trafficPersonVec.size());
    for (int person = 0; person < trafficPersonVec.size(); person++) {
    //for (int person = 0; person < 1; person++) {
      v1 = trafficPersonVec[person].init_intersection;
      v2 = trafficPersonVec[person].end_intersection;
      std::array<abm::graph::vertex_t, 2> od = {v1, v2};
      od_pairs.emplace_back(od);
    }
    if (nagents != std::numeric_limits<int>::max())
      od_pairs.resize(nagents);
    
    sort( od_pairs.begin(), od_pairs.end() );
    //od_pairs.erase( unique( od_pairs.begin(), od_pairs.end() ), od_pairs.end() );

    /*
    //make sure all OD pairs are unique
    std::sort(od_pairs.begin(), od_pairs.end());
    auto last = std::unique(od_pairs.begin(), od_pairs.end());
      // using default comparison:
  std::array<abm::graph::vertex_t, 2>::iterator it;
  it = std::unique (od_pairs.begin(), od_pairs.end());   // 10 20 30 20 10 ?  ?  ?  ?
                                                         //                ^

  od_pairs.resize( std::distance(od_pairs.begin(),it) ); // 10 20 30 20 10
  */
  // using predicate comparison:
  //std::unique (od_pairs.begin(), od_pairs.end(), myfunction);   // (no changes)
    //for (int i = 0; i < od_pairs.size(); i++) {
	    //printf("vertex %d = %llu %llu\n", i, od_pairs[i].first, od_pairs[i].second);
	    //printf("vertex %d = %llu %llu\n", i, std::get<0>(od_pairs[i]), std::get<1>(od_pairs[i]));
    //}
    //od_pairs.erase(last, od_pairs.end());
    /*
    for (int i = 0; i < od_pairs.size(); i++) {
	    printf("vertex %d = %llu %llu\n", i, std::get<0>(od_pairs[i]), std::get<1>(od_pairs[i]));
    }
    */

  } catch (std::exception& exception) {
    std::cout << "Looping through trafficPersonVec doesn't work " << exception.what() << "\n";
    status = false;
  }
  return od_pairs;
}

// Read OD pairs file format
std::vector<std::array<abm::graph::vertex_t, 2>> B18TrafficSP::read_od_pairs(const std::string& filename, int nagents) {
  bool status = true;
  std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs;
  try {
    csvio::CSVReader<2> in(filename);
    in.read_header(csvio::ignore_extra_column, "origin", "destination");
    abm::graph::vertex_t v1, v2;
    abm::graph::weight_t weight;
    while (in.read_row(v1, v2)) {
      //std::array<abm::graph::vertex_t, 2> od = {v1, v2};
      std::array<abm::graph::vertex_t, 2> od = {v1, v2};
      //printf("v1 = %u, v2 = %u\n", v1, v2);
      od_pairs.emplace_back(od);
      RoadGraphB2018::demandB2018.push_back(DemandB2018(1, v1, v2)); //there is only one person for each OD pair
    }
    RoadGraphB2018::totalNumPeople = RoadGraphB2018::demandB2018.size();
    if (nagents != std::numeric_limits<int>::max())
      od_pairs.resize(nagents);
  } catch (std::exception& exception) {
    std::cout << "Read OD file: " << exception.what() << "\n";
    status = false;
  }
  return od_pairs;
}

// Read OD pairs file format
std::vector<float> B18TrafficSP::read_dep_times(const std::string& filename) {
  bool status = true;
  std::vector<float> dep_time_vec;
  try {
    csvio::CSVReader<1> in(filename);
    in.read_header(csvio::ignore_extra_column, "dep_time");
    float dep_time;
    while (in.read_row(dep_time)) {
      //printf("dep time %f\n", dep_time);
      dep_time_vec.emplace_back(dep_time);
    }
  } catch (std::exception& exception) {
    std::cout << "Read OD file: " << exception.what() << "\n";
    status = false;
  }
  //printf("dep time size %d\n", dep_time_vec.size());
  return dep_time_vec;
}


void B18TrafficSP::filterODByHour(std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs, std::vector<float> dep_times, float start_time, float end_time, std::vector<std::array<abm::graph::vertex_t, 2>> &filtered_od_pairs_, std::vector<float> &filtered_dep_times_) {
    filtered_od_pairs_.clear();
    filtered_dep_times_.clear();
    int filt_size = 0;
    for (int x = 0; x < od_pairs.size(); x++) {
        //printf("dep time = %f\n", dep_times[x]);
        //printf("od_pair v1 = %u\n, od_pair v2 = %u\n", od_pairs[x][0], od_pairs[x][1]);
        if ((dep_times[x] >= start_time) && (dep_times[x] < end_time)) {
            //printf("dep time = %f\n", dep_times[x]);
            std::array<abm::graph::vertex_t, 2> od = {od_pairs[x][0], od_pairs[x][1]};
            filtered_od_pairs_.emplace_back(od);
            //printf("od_pairs v1 = %d, od_pairs v2 = %d\n", filtered_od_pairs_[filt_size][0], filtered_od_pairs_[filt_size][1]);
            filtered_dep_times_.emplace_back(dep_times[x]);
            filt_size++;
        }
    }
}

//void B18TrafficSP::convertVector(std::vector<abm::graph::vertex_t> paths_SP, std::vector<uint>& indexPathVec, std::map<std::shared_ptr<abm::Graph::Edge>, uint> &edgeDescToLaneMapNumSP, const std::shared_ptr<abm::Graph>& graph_) {
void B18TrafficSP::convertVector(std::vector<abm::graph::vertex_t> paths_SP, std::vector<uint>& indexPathVec, tsl::robin_map<std::shared_ptr<abm::Graph::Edge>, uint> &edgeDescToLaneMapNumSP, const std::shared_ptr<abm::Graph>& graph_) {

    //indexPathVec.clear();
    for (auto& x: paths_SP) {
        if (x != -1) {
            indexPathVec.emplace_back(edgeDescToLaneMapNumSP[graph_->edges_[graph_->edge_ids_to_vertices[x]]]);
        } else {
            indexPathVec.emplace_back(-1);
        }
    }
	//std::cout << "indexPathVec size = " << indexPathVec.size() << "\n";
}

void B18TrafficSP::createTimeMatrix(std::vector<std::vector<abm::graph::vertex_t>> pathsMatrix, std::vector<std::vector<float>>& timeMatrix, const std::shared_ptr<abm::Graph>& graph_) {
    //let's assume there is a vector of vectors that are output ([[all_paths1, all_paths2, all_paths3],])
    float routeTime = 0;
    for (int num_route = 0; num_route < pathsMatrix.size(); num_route++) {
        for (int edge = 0; edge < pathsMatrix[num_route].size(); edge++) {
            if (pathsMatrix[num_route][edge] != -1) {
                routeTime += graph_->edge_costs_[pathsMatrix[num_route][edge]];
            } else {
                timeMatrix[num_route].emplace_back(routeTime);
                routeTime = 0;
            }
        }
    }
}

int exponentiate(float x) {return exp(x);}

void B18TrafficSP::updateRouteShareMatrix(std::vector<std::vector<float>>& timeMatrix, std::vector<std::vector<float>>& routeShareMatrix) {

    std::vector<float> beta_vec = {.5, .3, .7}; //made up values for b1,b2,b3
    for (int edge = 0; edge < timeMatrix[edge].size(); edge++) {
        std::vector<float> util_vec;
        std::vector<float> prob_vec;
        for (int num_route = 0; num_route < timeMatrix.size(); num_route++) {
                float util = beta_vec[num_route]*timeMatrix[num_route][edge];
                util_vec.emplace_back(util);
        }
        
        std::vector<float> denominator_vec;
        denominator_vec.resize(util_vec.size()); // unfortunately this is necessary
        std::transform(util_vec.begin(), util_vec.end(), denominator_vec.begin(), exponentiate);
        float sum_util = accumulate(denominator_vec.begin(),denominator_vec.end(),0);
        for (int num_route = 0; num_route < timeMatrix.size(); num_route++) {
                float prob = exp(util_vec[num_route]) / sum_util;
                prob_vec.emplace_back(prob);
        }
        routeShareMatrix.emplace_back(prob_vec);
    }
}

void B18TrafficSP::createTimeMap(std::map<std::tuple<abm::graph::vertex_t, abm::graph::vertex_t>, std::vector<std::vector<abm::graph::vertex_t>>> &subgraphPathsMap, std::map<std::tuple<abm::graph::vertex_t, abm::graph::vertex_t>, std::vector<float>> &subgraphTimeMap, const std::shared_ptr<abm::Graph>& graph_) {
    for (auto const& x: subgraphPathsMap) {
        for (int num_route = 0; num_route < x.second.size(); num_route++) {
            subgraphTimeMap[x.first][num_route] += graph_->edge_costs_[graph_->edge_ids_[x.first]];
        }
    }
}

void B18TrafficSP::updateRouteShareMap(std::map<std::tuple<abm::graph::vertex_t, abm::graph::vertex_t>, std::vector<float>> &subgraphTimeMap, std::map<std::tuple<abm::graph::vertex_t, abm::graph::vertex_t>, std::vector<float>> &subgraphRouteShareMap) {
    std::vector<float> beta_vec = {.5, .3, .7}; //made up values for b1,b2,b3
    for (auto const& x: subgraphTimeMap) {
        std::vector<float> util_vec;
        std::vector<float> prob_vec;
        for (int num_route = 0; num_route < x.second.size(); num_route++) {
                float util = beta_vec[num_route]*x.second[num_route];
                util_vec.emplace_back(util);
        }
        
        std::vector<float> denominator_vec;
        denominator_vec.resize(util_vec.size()); // unfortunately this is necessary
        std::transform(util_vec.begin(), util_vec.end(), denominator_vec.begin(), exponentiate);
        float sum_util = accumulate(denominator_vec.begin(),denominator_vec.end(),0);
        for (int num_route = 0; num_route < denominator_vec.size(); num_route++) {
                float prob = exp(util_vec[num_route]) / sum_util;
                prob_vec.emplace_back(prob);
        }
        subgraphRouteShareMap[x.first] = prob_vec;
    }
}

int randNumGen(std::vector<float> odRouteShare) {
    int i, RandNumIndex;
    srand((unsigned)time(NULL));
    int N = 100;
 
    RandNumIndex = (int) N * rand() / (RAND_MAX + 1); 

    for (int index = 0; index < odRouteShare.size(); index++) {
            if (RandNumIndex <= odRouteShare[index]) {
                RandNumIndex = index;
                break;
            }
    }
    return RandNumIndex;
}
     
        
void B18TrafficSP::filterODByHourAndSubgraph(std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs,
                                             std::vector<float> dep_times,
                                             float start_time,
                                             float end_time,
                                             std::vector<std::array<abm::graph::vertex_t, 2>> &filtered_od_pairs_,
                                             std::vector<float> &filtered_dep_times_,
                                             std::map<std::tuple<abm::graph::vertex_t, abm::graph::vertex_t>, std::tuple<abm::graph::vertex_t, abm::graph::vertex_t>> &localToSubMap,
                                             std::map<std::tuple<abm::graph::vertex_t, abm::graph::vertex_t>, std::vector<std::vector<abm::graph::vertex_t>>> &subgraphPathsMap,
                                             std::map<std::tuple<abm::graph::vertex_t, abm::graph::vertex_t>, std::vector<float>> &subgraphRouteShareMap,
                                             std::vector<abm::graph::vertex_t> &reconstructed_all_paths,
                                             std::vector<abm::graph::vertex_t> &LtoS,
                                             std::vector<abm::graph::vertex_t> &StoL) {
        
        B18TrafficSP::filterODByHour(od_pairs, dep_times, start_time, end_time, filtered_od_pairs_, filtered_dep_times_);

        for (int edge_verts = 0; edge_verts < filtered_od_pairs_.size(); edge_verts++) {
            std::tuple<abm::graph::vertex_t, abm::graph::vertex_t> od = std::make_pair(filtered_od_pairs_[edge_verts][0], filtered_od_pairs_[edge_verts][1]);
            std::tuple<abm::graph::vertex_t, abm::graph::vertex_t> subgraph_od = localToSubMap[od];

            std::vector<float> subgraph_od_route_share = subgraphRouteShareMap[subgraph_od];

            int route_index = randNumGen(subgraph_od_route_share);

            std::vector<abm::graph::vertex_t> od_chosen_path = subgraphPathsMap[subgraph_od][route_index];


            //keep reconstructing allPaths vector
            reconstructAllPaths(reconstructed_all_paths, od, LtoS, StoL, od_chosen_path);
            
        }
}

void B18TrafficSP::reconstructAllPaths(std::vector<abm::graph::vertex_t> &reconstructed_all_paths,
                                       std::tuple<abm::graph::vertex_t, abm::graph::vertex_t> od,
                                       std::vector<abm::graph::vertex_t> &LtoS,
                                       std::vector<abm::graph::vertex_t> &StoL,
                                       std::vector<abm::graph::vertex_t> od_chosen_path) {
    
    reconstructed_all_paths.emplace_back(std::get<0>(od));
    
    for (int x = 0; x < LtoS.size(); x++) {
        reconstructed_all_paths.emplace_back(LtoS[x]);
    }
        
    for (int x = 0; x < od_chosen_path.size(); x++) {
        reconstructed_all_paths.emplace_back(od_chosen_path[x]);
    }
    
    for (int x = 0; x < StoL.size(); x++) {
        reconstructed_all_paths.emplace_back(StoL[x]);
    }
  
    reconstructed_all_paths.emplace_back(std::get<1>(od));
    
    reconstructed_all_paths.emplace_back(-1);

}


std::vector<abm::graph::vertex_t> B18TrafficSP::compute_routes(int mpi_rank,
                                                          int mpi_size,
                                                          const std::shared_ptr<abm::Graph>& graph_,
                                                          const std::vector<std::array<abm::graph::vertex_t, 2>>& od_pairs) {
  //! All paths
  std::vector<abm::graph::vertex_t> all_paths_;
  //! All paths indices
  std::vector<std::array<abm::graph::vertex_t, 3>> all_paths_idx_;

  //std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs;

/*
#ifdef USE_MPI
  // Create MPI pair type
  MPI_Datatype pair_t;
  MPI_Type_vector(2, 1, 1, MPI_LONG_LONG_INT, &pair_t);
  MPI_Type_commit(&pair_t);

  // Calculate chunk size to split router
  int chunk_size = od_pairs.size() / mpi_size;
  MPI_Bcast(&chunk_size, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
  od_pairs.resize(chunk_size);
  // Send route chunks to different compute nodes
  MPI_Scatter(od_pairs.data(), chunk_size, pair_t, od_pairs.data(),
              od_pairs.size(), pair_t, 0, MPI_COMM_WORLD);

  // Calculate the remaining chunk of od_pairs and add to rank 0
  int chunk_remainder = od_pairs.size() % mpi_size;
  if (mpi_rank == 0) {
    od_pairs.insert(od_pairs.begin(), od_pairs.end() - chunk_remainder,
                    od_pairs.end());
  }
#else
*/
  //od_pairs = od_pairs;
//#endif
  // Paths (vector of edges)
  std::vector<abm::graph::vertex_t> paths;
  paths.reserve(graph_->nedges());
  // Indices of start of path and length for each agent
  std::vector<std::array<abm::graph::vertex_t, 3>> paths_idx;
  paths_idx.reserve(od_pairs.size());

#pragma omp parallel for schedule(dynamic)
  for (abm::graph::vertex_t i = 0; i < od_pairs.size(); ++i) {
    //std::cout << "index " << i << "\n";
    //const auto sp = graph_->dijkstra_edges(graph_->vertex_map_[od_pairs[i][0]], graph_->vertex_map_[od_pairs[i][1]]);
    const auto sp = graph_->dijkstra_edges(od_pairs[i][0], od_pairs[i][1]);
    //printf("od pair 1 o = %d, od pair 1 d = %d\n", od_pairs[i][0], od_pairs[i][1]);
#pragma omp critical
    {
      paths_idx.emplace_back(std::array<abm::graph::vertex_t, 3>(
          {i, static_cast<abm::graph::vertex_t>(paths.size()),
           static_cast<abm::graph::vertex_t>(sp.size())}));
            paths.insert(std::end(paths), std::begin(sp), std::end(sp));
    }
  }
  // Get all paths and indices
  std::cout << "paths size " << paths.size() << "\n";
  all_paths_ = abm::gather_vectors_ids(paths);
  // for (int i = 0; i < all_paths_.size(); i++) {
  //   printf("all_paths = %d\n", all_paths_[i]);
  //}
  all_paths_idx_ = abm::gather_vector_arrays(paths_idx);
/*
#ifdef USE_MPI
  MPI_Type_free(&pair_t);
#endif
*/
  return all_paths_;
}

}  // Closing namespace LC

