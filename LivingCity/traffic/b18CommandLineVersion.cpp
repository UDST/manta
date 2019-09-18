#include <iostream>
#include <cassert>

#include <qcoreapplication.h>
#include <string>

#include "b18CommandLineVersion.h"

#include "roadGraphB2018Loader.h"

#include "sp/graph.h"
#include "traffic/b18TrafficSP.h"

#ifdef B18_RUN_WITH_GUI
#include "./LivingCity/b18TestSimpleRoadAndOD.h"
#endif


namespace LC {



void B18CommandLineVersion::runB18Simulation() {
  std::cerr << "[Log] Loading configuration." << std::endl;
  SimulatorConfiguration simulatorConfiguration("./command_line_options.ini");

  std::cerr << "[Log] Initializing simulator." << std::endl;
  B18TrafficSimulator b18TrafficSimulator(simulatorConfiguration);

  std::cerr << "[Log] Starting simulation." << std::endl;
  if (simulatorConfiguration.UseCPU()) {
    b18TrafficSimulator.simulateInCPU_MultiPass();
  } else {
    b18TrafficSimulator.simulateInGPU();
  }
}


}  // namespace LC

