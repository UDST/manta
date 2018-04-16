#include "b18CommandLineVersion.h"

#include "roadGraphB2018Loader.h"
#include "b18TestSimpleRoadAndOD.h"

namespace LC {
  void B18CommandLineVersion::runB18Simulation() {
    printf("<<+ runB18Simulation\n");

    const bool kUseBasicTest = true;

    const float deltaTime = 0.5f;
    const float startDemandH = 7.30f;
    const float endDemandH = 9.00f;

    float startSimulationH = startDemandH;
    float endSimulationH = endDemandH;

    if (kUseBasicTest) {
      // SIMPLE TEST: Create a basic road and basic OD to run a simulation.
      printf("kUseBasicTest == true");
      B18TestSimpleRoadAndOD::generateTest(cg.roadGraph, b18TrafficSimulator.trafficPersonVec, startDemandH, endDemandH, nullptr);
      b18TrafficSimulator.initSimulator(deltaTime, &cg.roadGraph);
    } else {
      // B18 CODE: Normal Simulation
      printf("1. Load RoadGraph");
      RoadGraphB2018::loadB2018RoadGraph(cg.roadGraph);

      printf("2. Init Simulator");
      b18TrafficSimulator.initSimulator(deltaTime, &cg.roadGraph);

      printf("3. Create people");
      const int limitNumPeople = 1; // -1 to not limit
      b18TrafficSimulator.createB2018People(startDemandH, endDemandH, limitNumPeople);
    }

    printf("4. Simulate");
    const int numPasses = 1;
    const bool useJohnsonRouting = false;
    b18TrafficSimulator.simulateInCPU_MultiPass(numPasses, startSimulationH, endSimulationH, useJohnsonRouting);

    printf(">>++ runB18Simulation\n");
  }



}  // LC