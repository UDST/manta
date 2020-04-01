#include <QString>
#include <string>

#include "b18CommandLineVersion.h"

#include "src/benchmarker.h"

#include "roadGraphB2018Loader.h"
#include "qcoreapplication.h"

#include "sp/graph.h"
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
  bool saveFiles = settings.value("SAVE_FILES", true).toBool();
  const float a = settings.value("A", .8).toFloat();
  const float b = settings.value("B", .8).toFloat();
  const float T = settings.value("T", 1.5).toFloat();
  //const float s_0 = settings.value("s_0", 7.0).toFloat();
  const float s_0 = 1.9108366323843835;

  //const float deltaTime = 0.5f; //Time step of .5 seconds
  //const float startDemandH = 5.00f; //Start time for the simulation (hour)
  //const float endDemandH = 12.00f; //End time for the simulation (hour)

  float startSimulationH = startDemandH;
  float endSimulationH = endDemandH;


  Benchmarker graphLoadBench("Load graph task", 1);
  Benchmarker initBench("Initialize traffic simulator task", 1);
  Benchmarker peopleBench("People creation task", 1);
  Benchmarker simulationBench("Simulation task", 1);

  ClientGeometry cg;
  B18TrafficSimulator b18TrafficSimulator(deltaTime, &cg.roadGraph);
  std::vector<abm::graph::vertex_t> all_paths;
  std::vector<std::array<abm::graph::vertex_t, 2>> all_od_pairs_;
  std::vector<float> dep_times_;

  const bool directed = true;
  auto street_graph = std::make_shared<abm::Graph>(directed);
  if (useSP) {
	  //make the graph from edges file and load the OD demand from od file
	  std::string odFileName = RoadGraphB2018::loadABMGraph(networkPathSP, street_graph, (int) startDemandH, (int) endDemandH);
	  all_od_pairs_ = B18TrafficSP::read_od_pairs(odFileName, std::numeric_limits<int>::max());
	  dep_times_ = B18TrafficSP::read_dep_times(odFileName);

      //std::cout << "first pair second vertex = " << all_od_pairs_[0][2] << "\n";
      //std::cout << "first dep time = " << dep_times_[0] << "\n";
	  //printf("# of OD pairs = %d\n", all_od_pairs_.size());

	  //compute the routes for every OD pair
	  int mpi_rank = 0;
	  int mpi_size = 1;
	  //auto start = high_resolution_clock::now();
      QTime timer_shortest_path;
      timer_shortest_path.start();

	  //create a set of people for simulation (trafficPersonVec)
	  b18TrafficSimulator.createB2018PeopleSP(startDemandH, endDemandH, limitNumPeople, addRandomPeople, street_graph, a, b, T, all_od_pairs_);

  } else {
	  graphLoadBench.startMeasuring();
	  RoadGraphB2018::loadB2018RoadGraph(cg.roadGraph, networkPath);
	  graphLoadBench.stopAndEndBenchmark();
  
	  peopleBench.startMeasuring();
	  b18TrafficSimulator.createB2018People(startDemandH, endDemandH, limitNumPeople, addRandomPeople, useSP);
	  peopleBench.stopAndEndBenchmark();

  }

  simulationBench.startMeasuring();
  if (useCPU) {
    b18TrafficSimulator.simulateInCPU_MultiPass(numOfPasses, startSimulationH, endSimulationH,
        useJohnsonRouting);
  } else {
	  //if useSP, convert all_paths to indexPathVec format and run simulation
    b18TrafficSimulator.simulateInGPU(numOfPasses, startSimulationH, endSimulationH,
        useJohnsonRouting, useSP, street_graph, all_paths, saveFiles, s_0, all_od_pairs_, dep_times_, usePrevPaths, networkPathSP);
  }
  simulationBench.stopAndEndBenchmark();
}



}  // LC
