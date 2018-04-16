#include "b18CommandLineVersion.h"

#include "roadGraphB2018Loader.h"

namespace LC {
  void B18CommandLineVersion::runB18Simulation() {
    printf("<<+ runB18Simulation\n");

    // 1. Load RoadGraph
    printf("1. Load RoadGraph");
    RoadGraphB2018::loadB2018RoadGraph(cg.roadGraph);

    printf("2. Init Simulator");
    const float deltaTime = 0.5f;
    b18TrafficSimulator.initSimulator(deltaTime, &cg.roadGraph);
    
    printf("3. Create people"); 
    const int limitNumPeople = 1; // -1 to not limit
    const float startDemandH = 7.30f;
    const float endDemandH = 9.00f;
    b18TrafficSimulator.createB2018People(startDemandH, endDemandH, limitNumPeople);

    printf("4. Simulate");
    float startTimeH = 7.30f;
    float endTimeH = endDemandH;

    int numPasses = 1;
    bool useJohnsonRouting = false;
    b18TrafficSimulator.simulateInCPU_MultiPass(numPasses, startTimeH, endTimeH, useJohnsonRouting);
    printf(">>++ runB18Simulation\n");
  }



}  // LC