
/************************************************************************************************
*		@desc Class to load the B2018 road graph
*		@author igarciad
************************************************************************************************/
#pragma once

#include "RoadGraph/roadGraph.h"
#include <QString>
#include "traffic/b18TrafficSP.h"
#include "traffic/sp/graph.h"

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

  static void loadB2018RoadGraph(
      std::shared_ptr<RoadGraph> inRoadGraph,
      const QString & networkPath,
      std::map<RoadGraph::roadGraphVertexDesc, uchar> & intersection_types);

  static std::string loadABMGraph(const std::string& networkPath, const std::shared_ptr<abm::Graph>& graph_);
  static std::vector<DemandB2018> demandB2018;
  static int totalNumPeople;
  static QHash<int, uint64_t> indToOsmid;

};

}
