/************************************************************************************************
*
*		LC Project - B18 Traffic lane map
*
*
*		@desc Class that contains the structure of the lane maps
*		@author igaciad
*
************************************************************************************************/

#include <cassert>
#include <iomanip>
#include <ios>
#include <cmath>
#include <map>
#include "b18TrafficLaneMap.h"
#include "sp/graph.h"

#include "./b18TrafficLaneMap.h"
#include "./laneCoordinatesComputer.h"
#include "OSMConstants.h"

namespace LC {

void addBlockingToIntersection(
    const uint vertexIdx,
    const std::vector<LC::Intersection> & updatedIntersections,
    std::vector<LC::Connection> & connections,
    std::vector<uint> & connectionsBlocking,
    LaneCoordinatesComputer & laneCoordinatesComputer)
{
  // First compute each lanes coordinates
  const std::unordered_map<uint, BoostPoint> lanesCoordinates =
    laneCoordinatesComputer.computeLanesCoordinatesFor(vertexIdx);

  if (lanesCoordinates.size() == 0)
    return;

  const Intersection & intersection = updatedIntersections.at(vertexIdx);
  for (
      uint mainConnectionIdx = intersection.connectionGraphStart;
      mainConnectionIdx < intersection.connectionGraphEnd;
      mainConnectionIdx++) {
    Connection & mainConnection = connections.at(mainConnectionIdx);
    mainConnection.connectionsBlockingStart = connectionsBlocking.size();
    const BoostSegment mainSegment{
      lanesCoordinates.at(mainConnection.inLaneLcId),
      lanesCoordinates.at(mainConnection.outLaneLcId)
    };

    for (
        uint secondConnectionIdx = intersection.connectionGraphStart;
        secondConnectionIdx < intersection.connectionGraphEnd;
        secondConnectionIdx++) {
      if (mainConnectionIdx == secondConnectionIdx)
        continue;

      const Connection & secondConnection = connections.at(secondConnectionIdx);

      if (mainConnection.inLaneLcId == secondConnection.inLaneLcId
          || mainConnection.inLaneLcId == secondConnection.outLaneLcId
          || mainConnection.outLaneLcId == secondConnection.inLaneLcId
          || mainConnection.outLaneLcId == secondConnection.outLaneLcId)
        continue;

      const BoostSegment secondSegment{
        lanesCoordinates.at(secondConnection.inLaneLcId),
        lanesCoordinates.at(secondConnection.outLaneLcId)
      };

      std::vector<BoostPoint> intersectionPoints;
      boost::geometry::intersection(mainSegment, secondSegment, intersectionPoints);
      if (intersectionPoints.size() == 0)
        continue;

      // If there is an intersection between the main conneciton and the second connection,
      // then add the second connection index to the list of connections blocked by the main one.
      connectionsBlocking.push_back(secondConnectionIdx);
    }
    mainConnection.connectionsBlockingEnd = connectionsBlocking.size();
  }
}


void addTrafficLightScheduleToIntersection(
    Intersection & tgtIntersection,
    const long long vertexIdx,
    std::vector<TrafficLightScheduleEntry> & trafficLightSchedules,
    const std::vector<LC::Connection> & connections,
    const std::vector<B18EdgeData> & edges)
{
  tgtIntersection.trafficLightSchedulesStart = trafficLightSchedules.size();

  std::unordered_set<uint> indexesOfNotYetScheduledConnections;
  for (
      uint connectionIdx = tgtIntersection.connectionGraphStart;
      connectionIdx < tgtIntersection.connectionGraphEnd;
      ++connectionIdx) {
    indexesOfNotYetScheduledConnections.insert(connectionIdx);
  }

  const auto sourceVertex = [&connections, &edges] (const uint connectionIdx) {
    return edges.at(connections.at(connectionIdx).inEdgeLcId).sourceVertexLcId;
  };

  const auto targetVertex = [&connections, &edges] (const uint connectionIdx) {
    return edges.at(connections.at(connectionIdx).outEdgeLcId).targetVertexLcId;
  };

  const uint scheduledTime = 10;
  uint scheduleGroup = 0;
  while (!indexesOfNotYetScheduledConnections.empty()) {
    // Choose one not yet scheduled connection
    const auto currentIdxIt = indexesOfNotYetScheduledConnections.cbegin();
    const uint currentIdx = *currentIdxIt;
    indexesOfNotYetScheduledConnections.erase(currentIdxIt);

    // Find all the connections which are compatible with this connection
    std::unordered_set<uint> indexesOfCompatibleConnections{currentIdx};
    for (
        uint otherConnectionIdx = tgtIntersection.connectionGraphStart;
        otherConnectionIdx < tgtIntersection.connectionGraphEnd;
        ++otherConnectionIdx) {
      const bool isSameConnection = otherConnectionIdx == currentIdx;
      const bool goesInSameDirection =
        sourceVertex(otherConnectionIdx) == sourceVertex(currentIdx)
        && targetVertex(otherConnectionIdx) == targetVertex(currentIdx);
      const bool goesInOppositeDirection =
        sourceVertex(otherConnectionIdx) == targetVertex(currentIdx)
        && targetVertex(otherConnectionIdx) == sourceVertex(currentIdx);

      const bool isCompatible =
        !isSameConnection && (goesInOppositeDirection || goesInSameDirection);
      if (!isCompatible) continue;

      indexesOfCompatibleConnections.insert(otherConnectionIdx);
      indexesOfNotYetScheduledConnections.erase(otherConnectionIdx);
    }

    // Create a new schedule group with all the compatible connections
    for (const uint connectionIdx : indexesOfCompatibleConnections) {
      trafficLightSchedules.emplace(
        trafficLightSchedules.end(),
        vertexIdx,
        connectionIdx,
        scheduleGroup,
        scheduledTime);
    }
    scheduleGroup++;
  }

  tgtIntersection.timeOfNextUpdate = 0;
  tgtIntersection.scheduleIdx = tgtIntersection.trafficLightSchedulesStart;
  tgtIntersection.currentScheduleGroup = 0;
  tgtIntersection.trafficLightSchedulesEnd = trafficLightSchedules.size();
}


SimulatorDataInitializer::SimulatorDataInitializer(
        const std::shared_ptr<RoadGraph> & boost_street_graph_shared_ptr,
        const std::shared_ptr<abm::Graph> & abm_street_graph_shared_ptr,
        const SimulatorConfiguration & configuration) :
    boost_street_graph_shared_ptr_(boost_street_graph_shared_ptr),
    abm_street_graph_shared_ptr_(abm_street_graph_shared_ptr),
    use_boost_graph_(configuration.SimulationRouting() != Routing::SP) {}

void SimulatorDataInitializer::initializeDataStructures(
    std::vector<uchar> &laneMap,
    std::vector<B18EdgeData> &edgesData,
    std::vector<B18IntersectionData> &intersections,
    std::vector<uchar> &trafficLights,
    std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> & edgeDescToLaneMapNum,
    std::map<uint, RoadGraph::roadGraphEdgeDesc_BI> & laneMapNumToEdgeDesc,
    std::map<std::shared_ptr<abm::Graph::Edge>, uint> & edgeDescToLaneMapNumSP,
    std::map<uint, std::shared_ptr<abm::Graph::Edge>> & laneMapNumToEdgeDescSP,
    std::vector<LC::Connection> &connections,
    std::vector<uint> &connectionsBlocking,
    std::vector<LC::Intersection> &updatedIntersections,
    std::vector<TrafficLightScheduleEntry> &trafficLightSchedules,
    std::vector<uint> &inLanesIndexes) const
{
  const auto & boost_input_graph = boost_street_graph_shared_ptr_->myRoadGraph_BI;
  const int amount_of_edges = use_boost_graph_
    ? boost::num_edges(boost_input_graph)
    : abm_street_graph_shared_ptr_->nedges();
  const int amount_of_vertices = use_boost_graph_
    ? boost::num_vertices(boost_input_graph)
    : abm_street_graph_shared_ptr_->amount_of_vertices_;

  edgesData.resize(amount_of_edges * 4);  //4 to make sure it fits

  edgeDescToLaneMapNum.clear();
  laneMapNumToEdgeDesc.clear();
  edgeDescToLaneMapNumSP.clear();
  laneMapNumToEdgeDescSP.clear();
  connections.clear();
  connectionsBlocking.clear();
  updatedIntersections.clear();
  trafficLightSchedules.clear();
  inLanesIndexes.clear();

  RoadGraph::roadGraphEdgeIter_BI ei, ei_end;

  updatedIntersections.resize(amount_of_vertices);

  int totalLaneMapChunks = 0;
  const auto initialize_edge_data = [&] (
      const uint roadAmountOfLanes,
      const float edgeLength,
      const uint sourceVertexLcId,
      const uint targetVertexLcId,
      const float maxSpeedMperSec,
      const bool startsAtHighway
    ) {

    if (roadAmountOfLanes == 0) {
      return;
    }

    const int numWidthNeeded = static_cast<int>(std::ceil(edgeLength / kMaxMapWidthM));

    edgesData[totalLaneMapChunks].length = edgeLength;
    edgesData[totalLaneMapChunks].sourceVertexLcId = sourceVertexLcId;
    edgesData[totalLaneMapChunks].targetVertexLcId = targetVertexLcId;
    edgesData[totalLaneMapChunks].maxSpeedMperSec = maxSpeedMperSec;
    edgesData[totalLaneMapChunks].nextInters = targetVertexLcId;
    edgesData[totalLaneMapChunks].numLines = roadAmountOfLanes;
    edgesData[totalLaneMapChunks].valid = true;
    edgesData[totalLaneMapChunks].startsAtHighway = startsAtHighway;

    totalLaneMapChunks += roadAmountOfLanes * numWidthNeeded;
  };

  std::cerr << "[Log] Initializing edge data." << std::endl;
  if (use_boost_graph_) {
    for (boost::tie(ei, ei_end) = boost::edges(boost_input_graph); ei != ei_end; ++ei) {
      edgeDescToLaneMapNum.insert(std::make_pair(*ei, totalLaneMapChunks));
      laneMapNumToEdgeDesc.insert(std::make_pair(totalLaneMapChunks, *ei));

      initialize_edge_data(
        boost_input_graph[*ei].numberOfLanes,
        boost_input_graph[*ei].edgeLength,
        source(*ei, boost_input_graph),
        target(*ei, boost_input_graph),
        boost_input_graph[*ei].maxSpeedMperSec,
        boost_input_graph[source(*ei, boost_input_graph)].intersectionType == OSMConstant::MotorwayJunction);
    }
  } else {
    for (auto const& x : abm_street_graph_shared_ptr_->edges_) {
      edgeDescToLaneMapNumSP.insert(std::make_pair(x.second, totalLaneMapChunks));
      laneMapNumToEdgeDescSP.insert(std::make_pair(totalLaneMapChunks, x.second));

      initialize_edge_data(
        std::get<1>(x)->second[1],
        std::get<1>(x)->second[0],
        abm_street_graph_shared_ptr_->vertex_osm_ids_to_lc_ids_.at(std::get<0>(std::get<0>(x))),
        abm_street_graph_shared_ptr_->vertex_osm_ids_to_lc_ids_.at(std::get<1>(std::get<0>(x))),
        std::get<1>(x)->second[2],
        abm_street_graph_shared_ptr_->vertex_OSM_type_.at(std::get<0>(std::get<0>(x))) == OSMConstant::MotorwayJunction);
    }
  }
  edgesData.resize(totalLaneMapChunks);

  const auto initialize_connections_between = [
    &connections
  ] (
      uint & connectionsCount,
      const uint & vertexLcId,
      const uint & inEdgeLcId,
      const B18EdgeData & inEdgeLcData,
      const uint & outEdgeLcId,
      const B18EdgeData & outEdgeLcData) {
    if (inEdgeLcData.sourceVertexLcId == outEdgeLcData.targetVertexLcId) {
      // Avoid U-turns
      return;
    }

    assert(inEdgeLcData.targetVertexLcId == outEdgeLcData.sourceVertexLcId);
    for (uint inIdx = 0; inIdx < inEdgeLcData.numLines; inIdx++) {
      for (uint outIdx = 0; outIdx < outEdgeLcData.numLines; outIdx++) {
        Connection connection;
        connection.vertexLcId = vertexLcId;
        connection.inEdgeLcId = inEdgeLcId;
        connection.outEdgeLcId = outEdgeLcId;
        connection.inLaneLcId = inEdgeLcId + inIdx;
        connection.outLaneLcId = outEdgeLcId + outIdx;
        connection.enabled = false;
        connections.push_back(connection);
        ++connectionsCount;
      }
    }
  };

  const CoordinatesRetriever coordinatesRetriever = [&] (const uint & vertexLcId) {
    const auto vertexGraphId = use_boost_graph_
      ? vertexLcId
      : abm_street_graph_shared_ptr_->vertex_lc_ids_to_osm_ids_.at(vertexLcId);

    const double x = use_boost_graph_
      ? boost_input_graph[vertexGraphId].x
      : abm_street_graph_shared_ptr_->vertices_positions_.at(vertexGraphId).x();

    const double y = use_boost_graph_
      ? boost_input_graph[vertexGraphId].y
      : abm_street_graph_shared_ptr_->vertices_positions_.at(vertexGraphId).y();

    return std::make_pair(x, y);
  };

  LaneCoordinatesComputer laneCoordinatesComputer{
    coordinatesRetriever,
    edgesData,
    connections,
    updatedIntersections};

  const auto initialize_updated_intersection = [&] (uint & connectionsCount, const uint vertexLcId) {
    const uint vertexGraphId = use_boost_graph_
      ? vertexLcId
      : abm_street_graph_shared_ptr_->vertex_lc_ids_to_osm_ids_.at(vertexLcId);

    Intersection & intersection = updatedIntersections.at(vertexLcId);

    // Check whether this intersection is a stop junction
    const auto intersectionOsmType = use_boost_graph_
      ? boost_input_graph[vertexGraphId].intersectionType
      : abm_street_graph_shared_ptr_->vertex_OSM_type_.at(vertexGraphId);
    intersection.isStopIntersection = intersectionOsmType == OSMConstant::StopJunction;

    // Create connections information
    intersection.connectionGraphStart = connectionsCount;

    if (use_boost_graph_) {
      const auto in_edges_pair = boost::in_edges(vertexGraphId, boost_input_graph);
      const auto in_edges_begin = in_edges_pair.first;
      const auto in_edges_end = in_edges_pair.second;
      for (auto in_edges_it = in_edges_begin; in_edges_it != in_edges_end; ++in_edges_it) {
        const auto & inEdgeLcId = edgeDescToLaneMapNum.at(*in_edges_it);
        const auto & inEdgeLcData = edgesData.at(inEdgeLcId);

        const auto p2 = boost::out_edges(vertexGraphId, boost_input_graph);
        const auto outEdgesBegin = p2.first;
        const auto outEdgesEnd = p2.second;
        for (auto outEdgesIt = outEdgesBegin; outEdgesIt != outEdgesEnd; ++outEdgesIt) {
          const auto & outEdgeLcId = edgeDescToLaneMapNum.at(*outEdgesIt);
          const auto & outEdgeLcData = edgesData.at(outEdgeLcId);
          initialize_connections_between(
              connectionsCount, vertexLcId, inEdgeLcId, inEdgeLcData, outEdgeLcId, outEdgeLcData);
        }
      }
    } else {
      // TODO(ffigari): This could be optimized by using the adjacency lists
      for (const auto & inEdgeIt : abm_street_graph_shared_ptr_->edges_) {
        const auto & inEdgeGraphData = inEdgeIt.second;
        const auto inEdgeTargetGraphId = inEdgeGraphData->first.second;
        if (inEdgeTargetGraphId != vertexGraphId)
          continue;

        const auto & inEdgeLcId = edgeDescToLaneMapNumSP.at(inEdgeGraphData);
        const auto & inEdgeLcData = edgesData.at(inEdgeLcId);
        for (const auto & outEdgeIt : abm_street_graph_shared_ptr_->edges_) {
          const auto & outEdgeGraphData = outEdgeIt.second;
          const auto outEdgeSourceGraphId = outEdgeGraphData->first.first;
          if (outEdgeSourceGraphId != vertexGraphId)
            continue;

          const auto & outEdgeLcId = edgeDescToLaneMapNumSP.at(outEdgeGraphData);
          const auto & outEdgeLcData = edgesData.at(outEdgeLcId);

          initialize_connections_between(
              connectionsCount, vertexLcId, inEdgeLcId, inEdgeLcData, outEdgeLcId, outEdgeLcData);
        }
      }
    }
    intersection.connectionGraphEnd = connectionsCount;

    addTrafficLightScheduleToIntersection(
      intersection,
      vertexLcId,
      trafficLightSchedules,
      connections,
      edgesData);
    addBlockingToIntersection(
      vertexLcId,
      updatedIntersections,
      connections,
      connectionsBlocking,
      laneCoordinatesComputer);

    intersection.trafficControl = mapOSMConstantToTrafficControl(intersectionOsmType);

    intersection.inLanesIndexesStart = inLanesIndexes.size();
    std::unordered_set<uint> intersectionInLanesIndexes;
    for (
        uint connectionIdx = intersection.connectionGraphStart;
        connectionIdx < intersection.connectionGraphEnd;
        connectionIdx++) {
      const Connection & connection = connections.at(connectionIdx);
      intersectionInLanesIndexes.insert(connection.inLaneLcId);
    }
    for (const uint & index : intersectionInLanesIndexes) {
      inLanesIndexes.push_back(index);
    }
    intersection.inLanesIndexesEnd = inLanesIndexes.size();
  };

  std::cerr << "[Log] Initializing intersections data." << std::endl;
  uint connectionsCount = 0;
  if (use_boost_graph_) {
    auto p = boost::vertices(boost_input_graph);
    const auto verticesBegin = p.first;
    const auto verticesEnd = p.second;
    for (auto vertices_it = verticesBegin; vertices_it != verticesEnd; ++vertices_it) {
      const uint vertexLcId = *vertices_it;
      initialize_updated_intersection(connectionsCount, vertexLcId);
    }
  } else {
    for (const auto & vertex_mapping : abm_street_graph_shared_ptr_->vertex_osm_ids_to_lc_ids_) {
      const uint vertexLcId = vertex_mapping.second;
      initialize_updated_intersection(connectionsCount, vertexLcId);
    }
  }

  for (uint idx = 0; idx < connections.size(); ++idx) {
    const Connection & connection = connections.at(idx);
    for (
        uint i = connection.connectionsBlockingStart;
        i < connection.connectionsBlockingEnd;
        ++i) {
      const uint blockedConnectionIdx = connectionsBlocking.at(i);
      const Connection & blockedConnection = connections.at(blockedConnectionIdx);
      if (blockedConnection.connectionsBlockingStart == blockedConnection.connectionsBlockingEnd)
        continue;

      bool found = false;
      for (
          uint j = blockedConnection.connectionsBlockingStart;
          !found && j < blockedConnection.connectionsBlockingEnd;
          ++j) {
        const uint secondBlockedConnectionIdx = connectionsBlocking.at(j);
        found = secondBlockedConnectionIdx == idx;
      }
      assert(found && "Blocked connections should be symmetric.");
    }
  }

  // Instantiate lane map
  laneMap.resize(kMaxMapWidthM * totalLaneMapChunks * 2); // 2: to have two maps.
  memset(laneMap.data(), -1, laneMap.size()*sizeof(unsigned char)); //

  intersections.resize(amount_of_vertices);
  trafficLights.assign(totalLaneMapChunks, 0);

  std::cerr << "[Log] Initializing old intersections data." << std::endl;
  RoadGraph::roadGraphVertexIter_BI vi, viEnd;
  RoadGraph::in_roadGraphEdgeIter_BI Iei, Iei_end;
  RoadGraph::out_roadGraphEdgeIter_BI Oei, Oei_end;
  // TODO: Implement this so that it works with SP routing too
  return;
  if (!use_boost_graph_) throw std::runtime_error("Not yet implemented. #7");
  for (boost::tie(vi, viEnd) = boost::vertices(boost_input_graph); vi != viEnd; ++vi) {
    intersections.at(*vi).totalInOutEdges = boost::degree(*vi, boost_input_graph);

    if (intersections.at(*vi).totalInOutEdges <= 0) {
      printf("Vertex without in/out edges\n");
      continue;
    }

    if (intersections.at(*vi).totalInOutEdges >= 20) {
      printf("Vertex with more than 20 in/out edges\n");
      continue;
    }

    int numOutEdges = 0;
    std::vector<LC::RoadGraph::roadGraphEdgeDesc_BI> edgeAngleOut;
    for (boost::tie(Oei, Oei_end) = boost::out_edges(*vi, boost_input_graph); Oei != Oei_end; ++Oei) {
      if (boost_input_graph[*Oei].numberOfLanes == 0) {
        continue;
      }
      edgeAngleOut.push_back(*Oei);
      numOutEdges++;
    }

    int numInEdges = 0;
    std::vector<LC::RoadGraph::roadGraphEdgeDesc_BI> edgeAngleIn;
    for (boost::tie(Iei, Iei_end) = boost::in_edges(*vi, boost_input_graph); Iei != Iei_end; ++Iei) {
      if (boost_input_graph[*Iei].numberOfLanes == 0) {
        continue;
      }
      edgeAngleIn.push_back(*Iei);
      numInEdges++;
    }

    intersections.at(*vi).totalInOutEdges = numOutEdges + numInEdges;

    // Intersection data:
    //  Store the edges that go in or out of this intersection
    //  Said edges will be sorted by angle
    //
    //      0xFF00 0000 Num lines
    //      0x0080 0000 in out (one bit)
    //      0x007F FFFF Edge number
    size_t totalCount = 0;
    for (const auto & outEdgeIdx : edgeAngleOut) {
        assert(edgeDescToLaneMapNum[outEdgeIdx] < 0x007fffff && "Edge number is too high");
        intersections.at(*vi).edge[totalCount] = edgeDescToLaneMapNum[outEdgeIdx];
        intersections.at(*vi).edge[totalCount] |= (edgesData[intersections.at(*vi).edge[totalCount]].numLines << 24); //put the number of lines in each edge
        intersections.at(*vi).edge[totalCount] |= kMaskOutEdge; // 0x000000 mask to define out edge
        totalCount++;
    }
    for (const auto & inEdgeIdx : edgeAngleIn) {
        assert(edgeDescToLaneMapNum[inEdgeIdx] < 0x007fffff && "Edge number is too high");
        intersections.at(*vi).edge[totalCount] = edgeDescToLaneMapNum[inEdgeIdx];
        intersections.at(*vi).edge[totalCount] |= (edgesData[intersections.at(*vi).edge[totalCount]].numLines << 24); //put the number of lines in each edge
        intersections.at(*vi).edge[totalCount] |= kMaskInEdge; // 0x800000 mask to define in edge
        totalCount++;
    }
  }
}

void SimulatorDataInitializer::resetIntersections(
    std::vector<B18IntersectionData> & intersections,
    std::vector<uchar> & trafficLights) const {
  if (trafficLights.size() > 0) {
    memset(trafficLights.data(), 0, trafficLights.size()*sizeof(
             uchar));  //to make the system to repeat same execution
  }
}//

}//
