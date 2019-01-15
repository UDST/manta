#include <iostream>
#include <cassert>

#include <QString>
#include <qcoreapplication.h>

#include "b18CommandLineVersion.h"
#include "../roadGraphB2018Loader.h"

#ifdef B18_RUN_WITH_GUI
#include "./LivingCity/b18TestSimpleRoadAndOD.h"
#endif


namespace LC {


void B18CommandLineVersion::runB18Simulation() {
  QSettings settings("./command_line_options.ini", QSettings::IniFormat);

  const QString networkPath = settings.value("NETWORK_PATH").toString();
  const bool addRandomPeople = settings.value("ADD_RANDOM_PEOPLE", true).toBool();
  const bool useCPU = settings.value("USE_CPU", false).toBool();
  const bool useJohnsonRouting = settings.value("USE_JOHNSON_ROUTING", false).toBool();
  const float deltaTime = 0.5f;
  const float endDemandH = 12.00f;
  const float endSimulationH = endDemandH;
  const float startDemandH = 5.00f;
  const float startSimulationH = startDemandH;
  const int limitNumPeople = settings.value("LIMIT_NUM_PEOPLE", -1).toInt();
  const int numOfPasses = settings.value("NUM_PASSES", 1).toInt();

  ClientGeometry cg;
  RoadGraphB2018::loadB2018RoadGraph(cg.roadGraph, networkPath);

  B18TrafficSimulator b18TrafficSimulator(deltaTime, &cg.roadGraph);

  b18TrafficSimulator.createB2018People(startDemandH, endDemandH, limitNumPeople, addRandomPeople);

  if (useCPU) {
    b18TrafficSimulator.simulateInCPU_MultiPass(
      numOfPasses, startSimulationH, endSimulationH, useJohnsonRouting);
  } else {
    b18TrafficSimulator.simulateInGPU(
      numOfPasses, startSimulationH, endSimulationH, useJohnsonRouting);
  }
}


}  // namespace LC

