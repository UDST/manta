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

  QString networkPath = settings.value("NETWORK_PATH").toString();
  std::string networkPathSP = networkPath.toStdString();

  bool addRandomPeople = settings.value("ADD_RANDOM_PEOPLE", true).toBool();
  int limitNumPeople = settings.value("LIMIT_NUM_PEOPLE", -1).toInt(); // -1
  int numOfPasses = settings.value("NUM_PASSES", 1).toInt();

  //const float deltaTime = 0.5f; //Time step of 30 minutes
  //const float deltaTime = .0167f; //Time step of 1 minute
  const float deltaTime = .00028f; //Time step of 1 minute
  const float startDemandH = 5.00f; //Start time for the simulation (hour)
  const float endDemandH = 12.00f; //End time for the simulation (hour)

  float startSimulationH = startDemandH;
  float endSimulationH = endDemandH;


  Benchmarker graphLoadBench("Load graph task", 1);
  Benchmarker initBench("Initialize traffic simulator task", 1);
  Benchmarker peopleBench("People creation task", 1);
  Benchmarker simulationBench("Simulation task", 1);

  ClientGeometry cg;
  B18TrafficSimulator b18TrafficSimulator(deltaTime, &cg.roadGraph);
  //auto all_paths = std::vector<abm::graph::vertex_t>;
  std::vector<abm::graph::vertex_t> all_paths;
  const bool directed = true;
  auto street_graph = std::make_shared<abm::Graph>(directed);
  if (useSP) {
	  //make the graph from edges file and load the OD demand from od file
	  std::string odFileName = RoadGraphB2018::loadABMGraph(networkPathSP, street_graph);
	  const auto all_od_pairs_ = B18TrafficSP::read_od_pairs(odFileName, std::numeric_limits<int>::max());
	  printf("# of OD pairs = %d\n", all_od_pairs_.size());

	  //compute the routes for every OD pair
	  int mpi_rank = 0;
	  int mpi_size = 1;
	  auto start = high_resolution_clock::now();
	  all_paths = B18TrafficSP::compute_routes(mpi_rank, mpi_size, street_graph, all_od_pairs_);
	  auto stop = high_resolution_clock::now();
	  auto duration = duration_cast<milliseconds>(stop - start);
	  std::cout << "# of paths = " << all_paths.size() << "\n";
	  
      std::cout << "total time to compute shortest paths = " << duration.count() << "ms \n";

	  //create a set of people for simulation (trafficPersonVec)
	  b18TrafficSimulator.createB2018PeopleSP(startDemandH, endDemandH, limitNumPeople, addRandomPeople, street_graph);

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
        useJohnsonRouting, useSP, street_graph, all_paths);
  }
  simulationBench.stopAndEndBenchmark();
}



}  // LC
