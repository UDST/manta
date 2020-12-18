#pragma once
#include <QString>
#include <string>

#include "b18CommandLineVersion.h"

#include "src/benchmarker.h"

#include "roadGraphB2018Loader.h"
#include "qcoreapplication.h"

#include "sp/graph.h"
#include "traffic/b18TrafficSP.h"
#include "../roadGraphB2018Loader.h"
#include "accessibility.h"

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
  const std::string networkPathSP = networkPath.toStdString();

  bool addRandomPeople = settings.value("ADD_RANDOM_PEOPLE", true).toBool();
  int limitNumPeople = settings.value("LIMIT_NUM_PEOPLE", -1).toInt(); // -1
  int numOfPasses = settings.value("NUM_PASSES", 1).toInt();
  const float deltaTime = settings.value("TIME_STEP", .5).toFloat();
  const float startDemandH = settings.value("START_HR", 5).toFloat();
  const float endDemandH = settings.value("END_HR", 12).toFloat();
  const bool showBenchmarks = settings.value("SHOW_BENCHMARKS", false).toBool();
  int increment = settings.value("REROUTE_INCREMENT", 60).toInt(); //in minutes
  const parameters simParameters {
      settings.value("a",0.557040909258405).toDouble(),
      settings.value("b",2.9020578588167).toDouble(),
      settings.value("T",0.5433027817144876).toDouble(),
      settings.value("s_0",1.3807498735425845).toDouble()};


  std::string odDemandPath = settings.value("OD_DEMAND_FILENAME", "od_demand_5to12.csv").toString().toStdString();

  std::cout << "b18CommandLineVersion received the parameters "
            << "[a: " << simParameters.a 
            << ", b: " << simParameters.b
            << ", T: " << simParameters.T
            << ", s_0: " << simParameters.s_0
            << "]" << std::endl;

  float startSimulationH = startDemandH;
  float endSimulationH = endDemandH;

  // new benchmarks
  if (showBenchmarks){ Benchmarker::enableShowBenchmarks(); }
  Benchmarker loadNetwork("Load_network", true);
  Benchmarker loadODDemandData("Load_OD_demand_data", true);
  Benchmarker routingCH("Routing_CH", true);
  Benchmarker CHoutputNodesToEdgesConversion("CH_output_nodes_to_edges_conversion", true);
  Benchmarker initBench("Initialize traffic simulator task");
  Benchmarker peopleBench("People creation task");
  Benchmarker simulationBench("Simulation task");


  ClientGeometry cg;
  B18TrafficSimulator b18TrafficSimulator(deltaTime, &cg.roadGraph, simParameters);
  
  const bool directed = true;
  const std::shared_ptr<abm::Graph>& street_graph = std::make_shared<abm::Graph>(directed, networkPathSP);
  std::vector<abm::graph::edge_id_t> all_paths;
  std::vector<std::vector<int>> all_paths_ch;
  std::string odFileName = RoadGraphB2018::loadABMGraph(networkPathSP, odDemandPath, street_graph, (int) startDemandH, (int) endDemandH);
  const auto all_od_pairs_ = B18TrafficSP::read_od_pairs(odFileName, std::numeric_limits<int>::max());
  const auto dep_times = B18TrafficSP::read_dep_times(odFileName);
  std::vector<std::array<abm::graph::vertex_t, 2>> filtered_od_pairs_;
  std::vector<float> filtered_dep_times_;
  float newEndTimeH = ( startSimulationH + (increment / 60)) * 3600; 
  if (useSP) {
	  //make the graph from edges file and load the OD demand from od file
    loadNetwork.startMeasuring();
    loadNetwork.stopAndEndBenchmark();
    loadODDemandData.startMeasuring();
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

      B18TrafficSP::filterODByHour(all_od_pairs_,
                                  dep_times,
                                  startSimulationH * 3600,
                                  newEndTimeH,
                                  filtered_od_pairs_,
                                  filtered_dep_times_);

      for (int x = 0; x < filtered_od_pairs_.size(); x++) {
        sources.emplace_back(filtered_od_pairs_[x][0]);
        targets.emplace_back(filtered_od_pairs_[x][1]);
        //std::cout << "origin = " << sources[x] << " \n";
        //std::cout << "dest = " << targets[x] << " \n";
      }

      //get the edge lengths
      std::vector<double> edge_weights_inside_vec;
      for (auto const& x : street_graph->edges_) {
        double edge_impedance = std::get<1>(x)->second[0];
        //TODO(pavan, juli): also retrieve the speed limit of the edge and then we can divide to get time
        edge_weights_inside_vec.emplace_back(edge_impedance);
        //std::cout << "edge length = " << edge_impedance << " \n";

        std::vector<long> edge_nodes = {std::get<0>(x.first), std::get<1>(x.first)};
        edge_vals.emplace_back(edge_nodes);
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

      //write paths to file so that we can just load them instead
      const std::string& pathsFileName = networkPathSP + "all_paths_ch.txt";
      std::cout << "Save " << pathsFileName << " as paths file\n";
      std::ofstream output_file(pathsFileName);
      std::ostream_iterator<abm::graph::vertex_t> output_iterator(output_file, "\n");
      std::copy(all_paths.begin(), all_paths.end(), output_iterator);
    }
    //map person to their initial edge
    int count = 0;
    for (int i = 0; i < all_paths.size(); i++) {
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
        useJohnsonRouting, useSP, street_graph, all_paths, simParameters, increment, all_od_pairs_, dep_times);
  }

}



}  // LC
