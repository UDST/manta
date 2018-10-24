#include <QString>

#include "b18CommandLineVersion.h"

#include "benchmarker.h"

#include "roadGraphB2018Loader.h"
#include "qcoreapplication.h"

#ifdef B18_RUN_WITH_GUI
#include "b18TestSimpleRoadAndOD.h"
#endif

namespace LC {
void B18CommandLineVersion::runB18Simulation() {
  QSettings settings(QCoreApplication::applicationDirPath() + "/command_line_options.ini",
      QSettings::IniFormat);
  QString networkPath = settings.value("NETWORK_PATH").toString();
  bool addRandomPeople = settings.value("ADD_RANDOM_PEOPLE", true).toBool();
  bool useCPU = settings.value("USE_CPU", false).toBool();
  bool useJohnsonRouting = settings.value("USE_JOHNSON_ROUTING", false).toBool();
  int limitNumPeople = settings.value("LIMIT_NUM_PEOPLE", -1).toInt(); // -1
  int numOfPasses = settings.value("NUM_PASSES", 1).toInt();

  const float deltaTime = 0.5f;
  const float startDemandH = 5.00f;
  const float endDemandH = 12.00f;

  float startSimulationH = startDemandH;
  float endSimulationH = endDemandH;


  Benchmarker graphLoadBench("Load graph task", 1);
  Benchmarker initBench("Initialize traffic simulator task", 1);
  Benchmarker peopleBench("People creation task", 1);
  Benchmarker simulationBench("Simulation task", 1);

  graphLoadBench.startMeasuring();
  ClientGeometry cg;
  RoadGraphB2018::loadB2018RoadGraph(cg.roadGraph, networkPath);
  graphLoadBench.stopAndEndBenchmark();

  initBench.startMeasuring();
  B18TrafficSimulator b18TrafficSimulator(deltaTime, &cg.roadGraph);
  initBench.stopAndEndBenchmark();

  peopleBench.startMeasuring();
  b18TrafficSimulator.createB2018People(startDemandH, endDemandH, limitNumPeople, addRandomPeople);
  peopleBench.stopAndEndBenchmark();

  simulationBench.startMeasuring();
  if (useCPU) {
    b18TrafficSimulator.simulateInCPU_MultiPass(numOfPasses, startSimulationH, endSimulationH,
        useJohnsonRouting);
  } else {
    b18TrafficSimulator.simulateInGPU(numOfPasses, startSimulationH, endSimulationH,
        useJohnsonRouting);
  }
  simulationBench.stopAndEndBenchmark();
}



}  // LC
