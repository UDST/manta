
/************************************************************************************************
*		@desc Class to load the B2018 road graph
*		@author igarciad
************************************************************************************************/
#pragma once

#include "RoadGraph/roadGraph.h"

namespace LC {


class DemandB2018 {
public:
  int num_people;
  int src_vertex;
  int tgt_vertex;
  DemandB2018(int num_people, int src_vertex, int tgt_vertex): 
    num_people(num_people), src_vertex(src_vertex), tgt_vertex(tgt_vertex) {}
};

/**
* RoadGraph.
**/
class RoadGraphB2018 {

 public:

  /**
  * Load
  **/
  static void loadB2018RoadGraph(RoadGraph &inRoadGraph, bool loadFullNetwork); // select big or small network.
  
  static std::vector<DemandB2018> demandB2018;
  static int totalNumPeople;
 private:

};

}
