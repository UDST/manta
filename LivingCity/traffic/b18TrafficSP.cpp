#include "b18TrafficSP.h"

#include <boost/graph/exterior_property.hpp>
#include "src/linux_host_memory_logger.h"
#include "roadGraphB2018Loader.h"
#include "accessibility.h"

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
std::vector<std::array<abm::graph::vertex_t, 2>> B18TrafficSP::read_od_pairs(
  const std::string& filename, int nagents,
  const float startSimulationH,
  const float endSimulationH) {
  bool status = true;
  std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs;
  csvio::CSVReader<3> in(filename);
  in.read_header(csvio::ignore_extra_column, "dep_time", "origin", "destination");
  abm::graph::vertex_t v1, v2;
  abm::graph::weight_t weight;
  float dep_time;
  int count_outside_filter = 0;
  while (in.read_row(dep_time, v1, v2)) {
    //std::array<abm::graph::vertex_t, 2> od = {v1, v2};
    if (dep_time >= startSimulationH * 3600 && dep_time < endSimulationH * 3600){
      std::array<abm::graph::vertex_t, 2> od = {v1, v2};
      od_pairs.emplace_back(od);
      RoadGraphB2018::demandB2018.push_back(DemandB2018(1, v1, v2)); //there is only one person for each OD pair
    } else {
      count_outside_filter++;
    }
  }
  if (count_outside_filter > 0){
    std::cout << "WARNING: Filtering " << count_outside_filter << " trips outside the input time range." << std::endl;
  }
  RoadGraphB2018::totalNumPeople = RoadGraphB2018::demandB2018.size();
  if (nagents != std::numeric_limits<int>::max())
    od_pairs.resize(nagents);

  B18TrafficSP::numberOfPeopleRouted = 0;
  return od_pairs;
}

uint B18TrafficSP::numberOfPeopleRouted = 0;

// Read OD pairs file format
std::vector<float> B18TrafficSP::read_dep_times(
  const std::string& filename,
  const float startSimulationH,
  const float endSimulationH) {
  std::vector<float> dep_time_vec;
  csvio::CSVReader<1> in(filename);
  in.read_header(csvio::ignore_extra_column, "dep_time");
  float dep_time;
  int count_outside_filter = 0;
  while (in.read_row(dep_time)) {
    //printf("dep time %f\n", dep_time);
    if (dep_time >= startSimulationH * 3600 && dep_time < endSimulationH * 3600){
      dep_time_vec.emplace_back(dep_time);
    } else {
      count_outside_filter++;
    }
  }
  if (count_outside_filter > 0) {
    std::cout << "WARNING: Filtering " << count_outside_filter << " trips outside the input time range.";
  }
  return dep_time_vec;
}

void B18TrafficSP::filterODByTimeRange(
    const std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs,
    const std::vector<float> dep_times_in_seconds,
    const float start_time_mins,
    const float end_time_mins,
    std::vector<abm::graph::vertex_t>& filtered_od_pairs_sources_,
    std::vector<abm::graph::vertex_t>& filtered_od_pairs_targets_,
    std::vector<float>& filtered_dep_times_,
    std::vector<uint>& indexPathVecOrder,
    std::vector<uint>& currentBatchIndexPathVecOrder) {
  
  if (od_pairs.size() != dep_times_in_seconds.size()){
    throw std::runtime_error("Input from trips should match in size.");
  }

  filtered_od_pairs_sources_.clear();
  filtered_od_pairs_targets_.clear();
  currentBatchIndexPathVecOrder.clear();

  filtered_dep_times_.clear();
  for (uint person_id = 0; person_id < od_pairs.size(); person_id++) {
    if ((dep_times_in_seconds[person_id] >= start_time_mins * 60) && (dep_times_in_seconds[person_id] < end_time_mins * 60)) {
      filtered_od_pairs_sources_.push_back(od_pairs[person_id][0]);
      filtered_od_pairs_targets_.push_back(od_pairs[person_id][1]);
      filtered_dep_times_.push_back(dep_times_in_seconds[person_id]);
      indexPathVecOrder.push_back(person_id);
      currentBatchIndexPathVecOrder.push_back(person_id);
    }
  }
}

std::vector<abm::graph::edge_id_t> B18TrafficSP::loadPrevPathsFromFile(
  const std::string networkPathSP){

  std::vector<abm::graph::edge_id_t> paths_SP;
  // open file    
  const std::string& pathsFileName = networkPathSP + "all_paths_ch.txt";
  std::cout << "Loading " << pathsFileName << " as paths file" << std::endl;
  std::ifstream inputFile(pathsFileName);
  // test file open   
  if (inputFile) {        
    abm::graph::vertex_t value;
    // read the elements in the file into a vector  
    while (inputFile >> value) {
      paths_SP.push_back(value);
    }
  } else {
    throw std::runtime_error("Could not load previous paths file.");
  }

  return paths_SP;
}

template<typename T>
bool allValuesAreDifferent(std::vector<T> v){
  std::sort(v.begin(), v.end());
  for (int i = 0; i < v.size()-1; i++){
    if (v[i] == v[i+1]){
      return false;
    }
  }
  return true;
}

std::vector<abm::graph::edge_id_t> B18TrafficSP::RoutingWrapper (
  const std::vector<std::array<abm::graph::vertex_t, 2>> all_od_pairs_,
  const std::shared_ptr<abm::Graph>& street_graph,
  const std::vector<float>& dep_times,
  const float start_time_mins,
  const float end_time_mins,
  int reroute_batch_number,
  std::vector<uint>& indexPathVecOrder,
  const bool savePaths,
  const std::string networkPathSP,
  std::vector<LC::B18TrafficPerson>& trafficPersonVec,
  std::vector<uint> &indexPathVec) {

  std::cout << "numberOfPeopleRouted " << B18TrafficSP::numberOfPeopleRouted << std::endl;

  // --------------------------- preprocessing ---------------------------
  std::vector<abm::graph::vertex_t> filtered_od_pairs_sources_;
  std::vector<abm::graph::vertex_t> filtered_od_pairs_targets_;
  std::vector<float> filtered_dep_times_;
  std::vector<uint> currentBatchIndexPathVecOrder;

  //filter the next set of od pair/departures in the next increment
  B18TrafficSP::filterODByTimeRange(all_od_pairs_,
                                    dep_times,
                                    start_time_mins,
                                    end_time_mins,
                                    filtered_od_pairs_sources_,
                                    filtered_od_pairs_targets_,
                                    filtered_dep_times_,
                                    indexPathVecOrder,
                                    currentBatchIndexPathVecOrder);
  

  std::cout << "Simulating trips with dep_time between "
    << int(start_time_mins/60) << ":" << int(start_time_mins) % 60
    << "(" << start_time_mins << " in minutes)"
    << " and " << int(end_time_mins/60) << ":" << int(end_time_mins) % 60
    << "(" << end_time_mins << " in minutes)" << std::flush;
  std::cout << ". Trips in this time range: " << filtered_od_pairs_sources_.size() << "/" << dep_times.size() << std::endl;

  std::vector<std::vector<long>> edges_routing;
  std::vector<std::vector<double>> edge_weights_routing;
  edges_routing.reserve(street_graph->edges_.size());
  edge_weights_routing.reserve(street_graph->edges_.size());
  std::vector<double> edge_weights_routing_inside_vec;
  for (auto const& x : street_graph->edges_) {
    // build node routes vector
    auto nodes = x.first;
    auto nodeFrom = std::get<0>(nodes);
    auto nodeTo = std::get<1>(nodes);
    std::vector<long> edge_nodes = {nodeFrom, nodeTo};
    edges_routing.emplace_back(edge_nodes);

    // build weight vector
    std::shared_ptr<abm::Graph::Edge> edge = x.second;
    abm::Edge_vals edge_vals = edge->second;
    double edge_weight = edge_vals.weight;
    assert(edge_weight > 0);
    edge_weights_routing_inside_vec.emplace_back(edge_weight);
  }
  edge_weights_routing.emplace_back(edge_weights_routing_inside_vec);
  std::cout << "# nodes = " << street_graph->vertices_data_.size() << std::endl;

  // --------------------------- routing ---------------------------

  Benchmarker routingCH("Routing_CH_batch_" + std::to_string(reroute_batch_number), true);
  routingCH.startMeasuring();
  MTC::accessibility::Accessibility *graph_ch = new MTC::accessibility::Accessibility((int) street_graph->vertices_data_.size(), edges_routing, edge_weights_routing, false);
  std::vector<std::vector<abm::graph::edge_id_t> > all_paths_ch = graph_ch->Routes(filtered_od_pairs_sources_, filtered_od_pairs_targets_, 0);
  routingCH.stopAndEndBenchmark();

  std::cout << "# of paths = " << all_paths_ch.size() << std::endl;

  // --------------------------- postprocessing ---------------------------

  Benchmarker CHoutputNodesToEdgesConversion("CH_output_nodes_to_edges_conversion_batch_" +
                                              std::to_string(reroute_batch_number), true);
  CHoutputNodesToEdgesConversion.startMeasuring();
  //convert from nodes to edges
  std::vector<abm::graph::edge_id_t> all_paths;
  for (int i=0; i < all_paths_ch.size(); i++) {
    for (int j=0; j < all_paths_ch[i].size()-1; j++) {
      auto vertex_from = all_paths_ch[i][j];
      auto vertex_to = all_paths_ch[i][j+1];
      auto one_edge = street_graph->edge_ids_[vertex_from][vertex_to];
      all_paths.emplace_back(one_edge);
    }
    all_paths.emplace_back(-1);
  }
  CHoutputNodesToEdgesConversion.stopAndEndBenchmark();

  if (savePaths) {
    const std::string& pathsFileName = networkPathSP + "all_paths_ch.txt";
    std::cout << "Save " << pathsFileName << " as paths file\n";
    std::ofstream output_file(pathsFileName);
    std::ostream_iterator<abm::graph::vertex_t> output_iterator(output_file, "\n");
    std::copy(all_paths.begin(), all_paths.end(), output_iterator);
  }

  std::cout
    << "currentBatchIndexPathVecOrder.size() " << currentBatchIndexPathVecOrder.size()
    << " | all_paths_ch.size() " << all_paths_ch.size() << std::endl;
  assert(currentBatchIndexPathVecOrder.size() == all_paths_ch.size());


  std::cout << "B18TrafficSP::numberOfPeopleRouted  " << B18TrafficSP::numberOfPeopleRouted << std::endl;
  // map person to their initial edge
  bool next_edge_is_init_edge = true;
  int person_index_in_all_paths = 0;
  int count = 0;
  for (int i = 0; i < all_paths.size(); i++) {
    if (next_edge_is_init_edge) {
      uint converted_person_id = currentBatchIndexPathVecOrder[person_index_in_all_paths];
      if (converted_person_id == 1060 || converted_person_id == 1899 || converted_person_id == 3898 ||
          converted_person_id == 3978 || converted_person_id == 1014384 || converted_person_id == 1020327 ||
          converted_person_id == 1023659 || converted_person_id == 1024049 || converted_person_id == 1024868) {

        std::cout << "person_index_in_all_paths " << person_index_in_all_paths << std::endl;
        std::cout << "currentBatchIndexPathVecOrder[person_index_in_all_paths] " << currentBatchIndexPathVecOrder[person_index_in_all_paths] << std::endl;
        std::cout << "i " << i << std::endl;
        std::cout << "i + indexPathVec.size() " << i + indexPathVec.size() << std::endl;
      }
      trafficPersonVec[currentBatchIndexPathVecOrder[person_index_in_all_paths]].indexPathInit = i + indexPathVec.size();
      person_index_in_all_paths++;
      next_edge_is_init_edge = false;
    }
    if (all_paths[i] == -1){
      next_edge_is_init_edge = true;
    }
  }
  assert(all_paths_ch.size() == person_index_in_all_paths);
  std::cout << "person_index_in_all_paths is " << person_index_in_all_paths << std::endl;
  std::cout << "all_paths_ch has size " << all_paths_ch.size();
  B18TrafficSP::numberOfPeopleRouted += all_paths_ch.size();

  return all_paths;
}


void B18TrafficSP::convertVector(std::vector<abm::graph::edge_id_t> paths_SP,
  std::vector<uint>& indexPathVec,  std::vector<uint> &edgeIdToLaneMapNum,
  const std::shared_ptr<abm::Graph>& graph_) {
  for (abm::graph::edge_id_t& edge_in_path: paths_SP) {
    if (edge_in_path != -1) {
      indexPathVec.emplace_back(edgeIdToLaneMapNum[edge_in_path]);
    } else {
      indexPathVec.emplace_back(-1);
    }
  }
	//std::cout << "indexPathVec size = " << indexPathVec.size() << "\n";
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
    //const auto sp = graph_->dijkstra_edges(od_pairs[i][0], od_pairs[i][1]);
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

