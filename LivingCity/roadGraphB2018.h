
/************************************************************************************************
*		@desc Class to load the B2018 road graph
*		@author igarciad
************************************************************************************************/
#pragma once

#include "RoadGraph/roadGraph.h"

namespace LC {

  class LCGLWidget3D;

  /**
  * RoadGraph.
  **/
  class RoadGraphB2018 {

  public:

    /**
    * Load
    **/
    static void loadB2018RoadGraph(RoadGraph& inRoadGraph, LCGLWidget3D* glWidget3D);
  private:

  };

}
