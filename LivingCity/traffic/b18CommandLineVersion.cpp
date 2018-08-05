#include "b18CommandLineVersion.h"

#include "roadGraphB2018Loader.h"
#include "qcoreapplication.h"

#ifdef B18_RUN_WITH_GUI
#include "b18TestSimpleRoadAndOD.h"
#endif

namespace LC {
void B18CommandLineVersion::runB18Simulation() {
  QSettings settings(QCoreApplication::applicationDirPath() +
                     "/command_line_options.ini",
                     QSettings::IniFormat);
  bool useBasicTest = settings.value("USE_BASIC_TEST",
                                     false).toBool(); // false = B2018; true = basic intersection
  bool useCPU = settings.value("USE_CPU",
                               false).toBool(); // false = GPU; true = CPU
  bool useFullB18Network = settings.value("USE_FULL_B2018_NETWORK",
                                          false).toBool(); // false = GPU; true = CPU
  bool useJohnsonRouting = settings.value("USE_JOHNSON_ROUTING",
                                          false).toBool(); // false = Disjktra; true = Johnson
  int limitNumPeople = settings.value("LIMIT_NUM_PEOPLE", -1).toInt(); // -1
  int numOfPasses = settings.value("NUM_PASSES", 1).toInt();

  const float deltaTime = 0.5f;
  const float startDemandH = 5.00f;
  const float endDemandH = 12.00f;

  float startSimulationH = startDemandH;
  float endSimulationH = endDemandH;

  ClientGeometry cg;

  /*
#ifdef B18_RUN_WITH_GUI

  if (useBasicTest) {
    // SIMPLE TEST: Create a basic road and basic OD to run a simulation.
    printf("useBasicTest == true");
    B18TestSimpleRoadAndOD::generateTest(cg.roadGraph,
        b18TrafficSimulator.trafficPersonVec, startDemandH, endDemandH, nullptr);
    b18TrafficSimulator.initSimulator(deltaTime, &cg.roadGraph);
  } else
#endif
  */
  // B18 CODE: Normal Simulation
  printf("1. Load RoadGraph\n");
  RoadGraphB2018::loadB2018RoadGraph(cg.roadGraph, useFullB18Network);

  printf("2. Init\n");
  B18TrafficSimulator b18TrafficSimulator(deltaTime, &cg.roadGraph);

  printf("3. Create people\n");
  b18TrafficSimulator.createB2018People(startDemandH, endDemandH, limitNumPeople);

  printf("4. Simulate\n");

  if (useCPU) {
    b18TrafficSimulator.simulateInCPU_MultiPass(numOfPasses, startSimulationH,
        endSimulationH, useJohnsonRouting);
  } else {
    b18TrafficSimulator.simulateInGPU(numOfPasses, startSimulationH, endSimulationH,
                                      useJohnsonRouting);
  }

  printf(">>++ runB18Simulation\n");
}



}  // LC
