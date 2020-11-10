#include <QString>
#include <string>

#include "b18CommandLineVersion.h"

#include "src/benchmarker.h"

#include "roadGraphB2018Loader.h"
#include "qcoreapplication.h"

#include "sp/graph.h"
#include "accessibility.h"
#include "traffic/b18TrafficSP.h"
#include "../roadGraphB2018Loader.h"

#ifdef B18_RUN_WITH_GUI
#include "b18TestSimpleRoadAndOD.h"
#endif

namespace LC {

using namespace std::chrono;

void B18CommandLineVersion::runB18Simulation() {
  QSettings settings(QCoreApplication::applicationDirPath() + "/command_line_options.ini",
      QSettings::IniFormat);
  bool useCPU = settings.value("USE_CPU", false).toBool();
  bool useJohnsonRouting = settings.value("USE_JOHNSON_ROUTING", false).toBool();
  bool useSP = settings.value("USE_SP_ROUTING", false).toBool();
  bool usePrevPaths = settings.value("USE_PREV_PATHS", false).toBool();

  QString networkPath = settings.value("NETWORK_PATH").toString();
  std::string networkPathSP = networkPath.toStdString();

  bool addRandomPeople = settings.value("ADD_RANDOM_PEOPLE", true).toBool();
  int limitNumPeople = settings.value("LIMIT_NUM_PEOPLE", -1).toInt(); // -1
  int numOfPasses = settings.value("NUM_PASSES", 1).toInt();
  const float deltaTime = settings.value("TIME_STEP", .5).toFloat();
  const float startDemandH = settings.value("START_HR", 5).toFloat();
  const float endDemandH = settings.value("END_HR", 12).toFloat();
  const bool showBenchmarks = settings.value("SHOW_BENCHMARKS", false).toBool();

  float startSimulationH = startDemandH;
  float endSimulationH = endDemandH;

  // new benchmarks
  if (showBenchmarks){ Benchmarker::enableShowBenchmarks(); }
  Benchmarker loadNetwork("Load_network", true);
  Benchmarker loadODDemandData("Load_OD_demand_data", true);
  Benchmarker routingCH("Routing_CH", true);
  Benchmarker CHoutputNodesToEdgesConversion("CH_output_nodes_to_edges_conversion", true);
  
  // ----

  Benchmarker initBench("Initialize traffic simulator task");
  Benchmarker peopleBench("People creation task");
  Benchmarker simulationBench("Simulation task");


  ClientGeometry cg;
  B18TrafficSimulator b18TrafficSimulator(deltaTime, &cg.roadGraph);
  
  std::vector<abm::graph::edge_id_t> all_paths;
  vector<vector<int>> all_paths_ch;
  const bool directed = true;
  auto street_graph = std::make_shared<abm::Graph>(directed, networkPathSP);
  if (useSP) {
	  //make the graph from edges file and load the OD demand from od file
    loadNetwork.startMeasuring();
	  std::string odFileName = RoadGraphB2018::loadABMGraph(networkPathSP, street_graph, (int) startDemandH, (int) endDemandH);
    loadNetwork.stopAndEndBenchmark();
    loadODDemandData.startMeasuring();
	  const auto all_od_pairs_ = B18TrafficSP::read_od_pairs(odFileName, std::numeric_limits<int>::max());
	  const auto dep_times = B18TrafficSP::read_dep_times(odFileName);
    loadODDemandData.stopAndEndBenchmark();

	  printf("# of OD pairs = %d\n", all_od_pairs_.size());

	  //compute the routes for every OD pair
	  int mpi_rank = 0;
	  int mpi_size = 1;
    if (usePrevPaths) {
      // open file    
      const std::string& pathsFileName = networkPathSP + "all_paths_ch.txt";
      std::cout << "Loading " << pathsFileName << " as paths file\n";
      std::ifstream inputFile(pathsFileName);
      // test file open   
      if (inputFile) {        
        abm::graph::vertex_t value;
        // read the elements in the file into a vector  
        while (inputFile >> value) {
          all_paths.push_back(value);
        }
      }
    } else {
      std::vector<std::vector<long>> edge_vals;

      std::vector<std::vector<double>> edge_weights;
      edge_weights.reserve(street_graph->edges_.size());

      std::vector<long> sources;
      std::vector<long> targets;
      sources.reserve(street_graph->edges_.size());
      targets.reserve(street_graph->edges_.size());


      for (int x = 0; x < all_od_pairs_.size(); x++) {
        sources.emplace_back(all_od_pairs_[x][0]);
        targets.emplace_back(all_od_pairs_[x][1]);
        //std::cout << "origin = " << sources[x] << " \n";
        //std::cout << "dest = " << targets[x] << " \n";
      }

      //get the edge lengths
      std::vector<double> edge_weights_inside_vec;
      for (auto const& x : street_graph->edges_) {
        double metersLength = std::get<1>(x)->second[0];
        edge_weights_inside_vec.emplace_back(metersLength);
        //std::cout << "edge length = " << metersLength << " \n";

        std::vector<long> edge_nodes = {std::get<0>(x.first), std::get<1>(x.first)};
        edge_vals.emplace_back(edge_nodes);
        //edge_vals.emplace_back((long) street_graph->edge_ids_[x.first]);
        //std::cout << "origin = " << sources[x] << " \n";

        //std::cout << "Start node id = " << edge_nodes[0] << " End node id = " << edge_nodes[1] << "\n";
      }
      edge_weights.emplace_back(edge_weights_inside_vec);
      std::cout << "# nodes = " << street_graph->vertices_data_.size() << "\n";

      routingCH.startMeasuring();
      MTC::accessibility::Accessibility *graph_ch = new MTC::accessibility::Accessibility((int) street_graph->vertices_data_.size(), edge_vals, edge_weights, false);

      all_paths_ch = graph_ch->Routes(sources, targets, 0);
      routingCH.stopAndEndBenchmark();
	    std::cout << "# of paths = " << all_paths_ch.size() << " \n";

      CHoutputNodesToEdgesConversion.startMeasuring();
      //convert from nodes to edges
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

      /*
      //add the -1 after every single route (for use in simulator)
      for (int i=0; i < all_paths_ch.size(); i++) {
          for (int j=0; j < all_paths_ch[i].size(); j++) {
              all_paths.emplace_back(all_paths_ch[i][j]);
          }
          all_paths.emplace_back(-1);
      }
      */
      /*
	    auto start = high_resolution_clock::now();
	    all_paths = B18TrafficSP::compute_routes(mpi_rank, mpi_size, street_graph, all_od_pairs_);
	    auto stop = high_resolution_clock::now();
	    auto duration = duration_cast<milliseconds>(stop - start);
        std::cout << "Shortest path time = " << duration.count() << " ms \n";
	    std::cout << "# of paths = " << all_paths.size() << "\n";
      */
      //write paths to file so that we can just load them instead
      const std::string& pathsFileName = networkPathSP + "all_paths_ch.txt";
      std::cout << "Save " << pathsFileName << " as paths file\n";
      std::ofstream output_file(pathsFileName);
      std::ostream_iterator<abm::graph::vertex_t> output_iterator(output_file, "\n");
      std::copy(all_paths.begin(), all_paths.end(), output_iterator);
    }
    //map person to their initial edge
    int count = 0;
    /*
    for (int i = 0; i < all_paths.size(); i++) {
      if ((all_paths[i] == -1) && (i == 0)) {
        street_graph->person_to_init_edge_[count] = -1;
        count++;
      } else if ((all_paths[i] == -1) && (all_paths[i+1] == -1)) {
        street_graph->person_to_init_edge_[count] = -1;
        count++;
      } else if ((all_paths[i] != -1) && (all_paths[i-1] == -1)) {
        street_graph->person_to_init_edge_[count] = all_paths[i];
        count++;
      } else if ((all_paths[i] == -1) && (i == (all_paths.size() - 1))) {
        break;
      }
    }
    */
    for (int i = 0; i < all_paths.size(); i++) {
      //if ((all_paths[i] == -1) && (i == 0)) {
      if (i == 0) { //first one that doesn't contain a -1 for logic
        street_graph->person_to_init_edge_[count] = i; 
        count++;
      } else if ((all_paths[i] == -1) && (all_paths[i+1] == -1)) { //if current is -1 and next is -1, increment (will result in nan)
        street_graph->person_to_init_edge_[count] = i;
        count++;
      } else if ((all_paths[i] != -1) && (all_paths[i-1] == -1)) { //if previous is -1, use this as first edge for p
        street_graph->person_to_init_edge_[count] = i;
        count++;
      } else if ((all_paths[i] == -1) && (i == (all_paths.size() - 1))) { //reach the end
        break;
      }
    }
    //std::cout << "person_to_init_edge size " << street_graph->person_to_init_edge_.size() << "\n";
    //std::cout << "# of paths = " << all_paths.size() << "\n";
    //std::cout << "Shortest path time = " << duration.count() << " ms \n";

    //create a set of people for simulation (trafficPersonVec)
    b18TrafficSimulator.createB2018PeopleSP(startDemandH, endDemandH, limitNumPeople, addRandomPeople, street_graph, dep_times);
  } else {
    RoadGraphB2018::loadB2018RoadGraph(cg.roadGraph, networkPath);
    b18TrafficSimulator.createB2018People(startDemandH, endDemandH, limitNumPeople, addRandomPeople, useSP);
  }

  
  if (useCPU) {
    b18TrafficSimulator.simulateInCPU_MultiPass(numOfPasses, startSimulationH, endSimulationH,
        useJohnsonRouting);
  } else {
	  //if useSP, convert all_paths to indexPathVec format and run simulation
    b18TrafficSimulator.simulateInGPU(numOfPasses, startSimulationH, endSimulationH,
        useJohnsonRouting, useSP, street_graph, all_paths);
  }

}



}  // LC
