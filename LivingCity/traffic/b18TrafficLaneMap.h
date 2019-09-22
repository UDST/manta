/************************************************************************************************
*
*		LC Project - B18 Traffic lane map
*
*
*		@desc Class that contains the structure of the lane maps
*		@author igaciad
*
************************************************************************************************/

#ifndef LC_B18_TRAFFIC_LANEMAP_H
#define LC_B18_TRAFFIC_LANEMAP_H


#include "RoadGraph/roadGraph.h"
#include "./b18EdgeData.h"
#include "traffic/sp/graph.h"
#include "./routing.h"
#include "./simulatorConfiguration.h"

namespace LC {


class SimulatorDataInitializer {

public:
    SimulatorDataInitializer(
        const std::shared_ptr<RoadGraph> & boost_street_graph_shared_ptr,
        const std::shared_ptr<abm::Graph> & abm_street_graph_shared_ptr,
        const SimulatorConfiguration & configuration);

    void initializeDataStructures(
        std::vector<uchar> &laneMap,
        std::vector<B18EdgeData> &edgesData,
        std::vector<uchar> &trafficLights,
        std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> & edgeDescToLaneMapNum,
        std::map<uint, RoadGraph::roadGraphEdgeDesc_BI> & laneMapNumToEdgeDesc,
        std::map<std::shared_ptr<abm::Graph::Edge>, uint> & edgeDescToLaneMapNumSP,
        std::map<uint, std::shared_ptr<abm::Graph::Edge>> & laneMapNumToEdgeDescSP,
        std::vector<LC::Connection> &connections,
        std::vector<uint> &connectionsBlocking,
        std::vector<LC::Intersection> &updatedIntersections,
        std::vector<TrafficLightScheduleEntry> &trafficLightSchedules,
        std::vector<uint> &inLanesIndexes) const;
    void resetIntersections(
        std::vector<uchar> &trafficLights) const;

private:
    const std::shared_ptr<RoadGraph> & boost_street_graph_shared_ptr_;
    const std::shared_ptr<abm::Graph> & abm_street_graph_shared_ptr_;
    const bool use_boost_graph_;
};

}

#endif  // LC_B18_TRAFFIC_LANEMAP_H

