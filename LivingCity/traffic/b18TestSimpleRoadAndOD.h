
/************************************************************************************************
*		@desc Class to create simple road and OD demmand (for simple testing)
*		@author igarciad
************************************************************************************************/
#pragma once

#include "../RoadGraph/roadGraph.h"
#include "cudaTrafficPerson.h"

namespace LC {

class LCGLWidget3D;

/**
* RoadGraph.
**/
class B18TestSimpleRoadAndOD {

 public:

  /**
  * Generate test: Road+People+OD
  **/
   static void generateTest(RoadGraph &inRoadGraph, std::vector<CUDATrafficPerson> &trafficPersonVec, 
     float startTimeH, float endTimeH, LCGLWidget3D *glWidget3D);
 private:

};

}
