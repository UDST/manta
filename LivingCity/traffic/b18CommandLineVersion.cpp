#include "b18CommandLineVersion.h"

#include "roadGraphB2018Loader.h"
#include "b18TestSimpleRoadAndOD.h"
#include <QApplication>

namespace LC {
  void B18CommandLineVersion::runB18Simulation() {
    printf("<<+ runB18Simulation\n");

    QSettings settings(QApplication::applicationDirPath() + "/command_line_options.ini", QSettings::IniFormat);
    bool useBasicTest = settings.value("USE_BASIC_TEST", false).toBool(); // false = B2018; true = basic intersection
    bool useCPU = settings.value("USE_CPU", false).toBool(); // false = GPU; true = CPU
    bool useFullB18Network = settings.value("USE_FULL_B2018_NETWORK", false).toBool(); // false = GPU; true = CPU
    bool useJohnsonRouting = settings.value("USE_JOHNSON_ROUTING", false).toBool(); // false = Disjktra; true = Johnson
    
    const float deltaTime = 0.5f;
    const float startDemandH = 7.30f;
    const float endDemandH = 9.00f;

    float startSimulationH = startDemandH;
    float endSimulationH = endDemandH;

    if (useBasicTest) {
      // SIMPLE TEST: Create a basic road and basic OD to run a simulation.
      printf("useBasicTest == true");
      B18TestSimpleRoadAndOD::generateTest(cg.roadGraph, b18TrafficSimulator.trafficPersonVec, startDemandH, endDemandH, nullptr);
      b18TrafficSimulator.initSimulator(deltaTime, &cg.roadGraph);
    } else {
      // B18 CODE: Normal Simulation
      printf("1. Load RoadGraph");
      RoadGraphB2018::loadB2018RoadGraph(cg.roadGraph, useFullB18Network);

      printf("2. Init Simulator");
      b18TrafficSimulator.initSimulator(deltaTime, &cg.roadGraph);

      printf("3. Create people");
      const int limitNumPeople = 1; // -1 to not limit
      b18TrafficSimulator.createB2018People(startDemandH, endDemandH, limitNumPeople);
    }

    printf("4. Simulate");
    const int numPasses = 1;
    if (useCPU) {
      b18TrafficSimulator.simulateInCPU_MultiPass(numPasses, startSimulationH, endSimulationH, useJohnsonRouting);
    } else {
      b18TrafficSimulator.simulateInGPU(startSimulationH, endSimulationH, useJohnsonRouting);
    }

    printf(">>++ runB18Simulation\n");
  }



}  // LC