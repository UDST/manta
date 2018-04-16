/************************************************************************************************
*
*		Command Line Version Of Traffic Simulator
*
*		@desc Class to execute the simulator for command line without opening the GUI.
*		@author igaciad
*
************************************************************************************************/

#pragma once

#include "traffic/b18TrafficSimulator.h"
#include "Geometry/client_geometry.h"

namespace LC {

  class B18CommandLineVersion{
  public:
    ClientGeometry cg;
    B18TrafficSimulator b18TrafficSimulator;

    void runB18Simulation();

  };

}  // namespace LC