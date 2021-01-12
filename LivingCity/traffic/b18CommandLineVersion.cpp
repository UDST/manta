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
#include <stdexcept>

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
  const float startSimulationH = settings.value("START_HR", 5).toFloat();
  const float endSimulationH = settings.value("END_HR", 12).toFloat();
  const bool showBenchmarks = settings.value("SHOW_BENCHMARKS", false).toBool();
  int rerouteIncrementMins = settings.value("REROUTE_INCREMENT", 60).toInt(); //in minutes

  if (rerouteIncrementMins <= 0){
    throw std::invalid_argument("Invalid reroute increment value.");
  }
    
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

  // new benchmarks
  if (showBenchmarks){ Benchmarker::enableShowBenchmarks(); }
  Benchmarker loadNetwork("Load_network", true);
  Benchmarker loadODDemandData("Load_OD_demand_data", true);
  Benchmarker initBench("Initialize traffic simulator task");
  Benchmarker peopleBench("People creation task");
  Benchmarker simulationBench("Simulation task");


  ClientGeometry cg;
  B18TrafficSimulator b18TrafficSimulator(deltaTime, &cg.roadGraph, simParameters);
  
  const bool directed = true;
  const std::shared_ptr<abm::Graph>& street_graph = std::make_shared<abm::Graph>(directed, networkPathSP);
  std::vector<abm::graph::edge_id_t> all_paths;
  std::vector<std::vector<int>> all_paths_ch;
  loadNetwork.startMeasuring();
  std::string odFileName = RoadGraphB2018::loadABMGraph(networkPathSP, odDemandPath, street_graph, (int) startSimulationH, (int) endSimulationH);
  loadNetwork.stopAndEndBenchmark();

  loadODDemandData.startMeasuring();
  const std::vector<std::array<abm::graph::vertex_t, 2>> all_od_pairs_ = B18TrafficSP::read_od_pairs(odFileName, std::numeric_limits<int>::max());
  const std::vector<float> dep_times = B18TrafficSP::read_dep_times(odFileName);
  loadODDemandData.stopAndEndBenchmark();
  std::vector<std::array<abm::graph::vertex_t, 2>> filtered_od_pairs_;
  std::vector<float> filtered_dep_times_;
  if (useSP) {
	  //make the graph from edges file and load the OD demand from od file
	  printf("# of OD pairs = %d\n", all_od_pairs_.size());

	  //compute the routes for every OD pair
	  int mpi_rank = 0;
	  int mpi_size = 1;
    if (usePrevPaths) {
      throw std::invalid_argument("prev_paths with reroute increment currently disabled.");
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
      const float startTimeMins = startSimulationH * 60;
      const float endTimeMins = startTimeMins + rerouteIncrementMins;
      cout << "startTime: " << startTimeMins << ", endTime: " << endTimeMins << endl;
      const int initialBatchNumber = 0;
      all_paths = B18TrafficSP::RoutingWrapper(all_od_pairs_, street_graph, dep_times, startTimeMins, endTimeMins, initialBatchNumber);
      std::cout << "paths size for batch #0: " << all_paths.size() << std::endl;
    }

    //std::cout << "person_to_init_edge size " << street_graph->person_to_init_edge_.size() << "\n";
    //std::cout << "# of paths = " << all_paths.size() << "\n";
    //std::cout << "Shortest path time = " << duration.count() << " ms \n";

    //create a set of people for simulation (trafficPersonVec)
    b18TrafficSimulator.createB2018PeopleSP(startSimulationH, endSimulationH, limitNumPeople, addRandomPeople, street_graph, dep_times);
  } else {
    RoadGraphB2018::loadB2018RoadGraph(cg.roadGraph, networkPath);
    b18TrafficSimulator.createB2018People(startSimulationH, endSimulationH, limitNumPeople, addRandomPeople, useSP);
  }

  
  if (useCPU) {
    b18TrafficSimulator.simulateInCPU_MultiPass(numOfPasses, startSimulationH, endSimulationH,
        useJohnsonRouting);
  } else {
	  //if useSP, convert all_paths to indexPathVec format and run simulation
    b18TrafficSimulator.simulateInGPU(numOfPasses, startSimulationH, endSimulationH,
        useJohnsonRouting, useSP, street_graph, all_paths, simParameters, rerouteIncrementMins, all_od_pairs_, dep_times);
  }

}



}  // LC
