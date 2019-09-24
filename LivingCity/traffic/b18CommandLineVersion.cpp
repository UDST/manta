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
  DataExporter dataExporter;

  std::cerr << "[Log] Loading configuration." << std::endl;
  SimulatorConfiguration simulatorConfiguration("./command_line_options.ini");

  std::cerr << "[Log] Initializing simulator." << std::endl;
  B18TrafficSimulator b18TrafficSimulator(simulatorConfiguration, dataExporter);

  std::cerr << "[Log] Starting simulation." << std::endl;
  b18TrafficSimulator.simulateInGPU();

  dataExporter.ExportTimes();
}


}  // namespace LC

