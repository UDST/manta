/************************************************************************************************
*
*		LC Project - B18 Traffic lane map
*
*
*		@desc Class that contains the structure of the lane maps
*		@author igaciad
*
************************************************************************************************/

#include <ios>
#include <cassert>
#include <map>
#include "./b18TrafficLaneMap.h"
#include <thrust/device_vector.h>

namespace LC {

B18TrafficLaneMap::B18TrafficLaneMap() {
}//

B18TrafficLaneMap::~B18TrafficLaneMap() {
}//

namespace {
  bool compareSecondPartTupleC(const
    std::pair<LC::RoadGraph::roadGraphEdgeDesc_BI, float> &i,
    const std::pair<LC::RoadGraph::roadGraphEdgeDesc_BI, float> &j) {
    return (i.second < j.second);
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
    return edges.at(connections.at(connectionIdx).inEdgeNumber).originalSourceVertexIndex;
  };

  const auto targetVertex = [&connections, &edges] (const uint connectionIdx) {
    return edges.at(connections.at(connectionIdx).outEdgeNumber).originalTargetVertexIndex;
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

void B18TrafficLaneMap::createLaneMap(
    const RoadGraph &inRoadGraph,
    std::vector<uchar> &laneMap,
    std::vector<B18EdgeData> &edgesData,
    std::vector<B18IntersectionData> &intersections,
    std::vector<uchar> &trafficLights,
    std::map<uint, RoadGraph::roadGraphEdgeDesc_BI> &laneMapNumToEdgeDesc,
    std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum,
    std::vector<LC::Connection> &connections,
    std::vector<LC::Intersection> &updatedIntersections,
    std::vector<TrafficLightScheduleEntry> &trafficLightSchedules) {
  edgesData.resize(boost::num_edges(inRoadGraph.myRoadGraph_BI) * 4);  //4 to make sure it fits

  edgeDescToLaneMapNum.clear();
  laneMapNumToEdgeDesc.clear();
  connections.clear();
  updatedIntersections.clear();
  trafficLightSchedules.clear();

  auto & inputGraph = inRoadGraph.myRoadGraph_BI;
  RoadGraph::roadGraphEdgeIter_BI ei, ei_end;

  updatedIntersections.resize(boost::num_vertices(inputGraph));

  int totalLaneMapChunks = 0;
  for (boost::tie(ei, ei_end) = boost::edges(inRoadGraph.myRoadGraph_BI); ei != ei_end; ++ei) {
    const int roadAmountOfLanes = inRoadGraph.myRoadGraph_BI[*ei].numberOfLanes;
    if (roadAmountOfLanes == 0) { continue; }

    const auto edgeLength = inRoadGraph.myRoadGraph_BI[*ei].edgeLength;
    const int numWidthNeeded = static_cast<int>(std::ceil(edgeLength / kMaxMapWidthM));

    edgesData[totalLaneMapChunks].originalSourceVertexIndex = source(*ei, inputGraph);
    edgesData[totalLaneMapChunks].originalTargetVertexIndex = target(*ei, inputGraph);
    edgesData[totalLaneMapChunks].length = edgeLength;
    edgesData[totalLaneMapChunks].maxSpeedMperSec = inRoadGraph.myRoadGraph_BI[*ei].maxSpeedMperSec;
    edgesData[totalLaneMapChunks].nextInters = boost::target(*ei, inRoadGraph.myRoadGraph_BI);
    edgesData[totalLaneMapChunks].numLines = roadAmountOfLanes;
    edgesData[totalLaneMapChunks].valid = true;

    edgeDescToLaneMapNum.insert(std::make_pair(*ei, totalLaneMapChunks));
    laneMapNumToEdgeDesc.insert(std::make_pair(totalLaneMapChunks, *ei));

    totalLaneMapChunks += roadAmountOfLanes * numWidthNeeded;
  }
  edgesData.resize(totalLaneMapChunks);

  auto p = boost::vertices(inputGraph);
  const auto verticesBegin = p.first;
  const auto verticesEnd = p.second;
  uint connectionsCount = 0;
  for (auto vertices_it = verticesBegin; vertices_it != verticesEnd; ++vertices_it) {
    const auto in_edges_pair = boost::in_edges(*vertices_it, inputGraph);
    const auto in_edges_begin = in_edges_pair.first;
    const auto in_edges_end = in_edges_pair.second;
    const auto vertexIdx = *vertices_it;

    Intersection & intersection = updatedIntersections.at(vertexIdx);

    // Create connections information
    intersection.connectionGraphStart = connectionsCount;
    for (auto in_edges_it = in_edges_begin; in_edges_it != in_edges_end; ++in_edges_it) {
      const auto & inEdgeNumber = edgeDescToLaneMapNum.at(*in_edges_it);
      const auto & inEdgeData = edgesData[inEdgeNumber];
      const auto p2 = boost::out_edges(*vertices_it, inputGraph);
      const auto outEdgesBegin = p2.first;
      const auto outEdgesEnd = p2.second;
      for (auto outEdgesIt = outEdgesBegin; outEdgesIt != outEdgesEnd; ++outEdgesIt) {
        const auto & outEdgeNumber = edgeDescToLaneMapNum.at(*outEdgesIt);
        const auto & outEdgeData = edgesData[outEdgeNumber];
        if (inEdgeData.originalSourceVertexIndex == outEdgeData.originalTargetVertexIndex) {
          // Avoid U-turns
          continue;
        }

        assert(inEdgeData.originalTargetVertexIndex == outEdgeData.originalSourceVertexIndex);
        for (uint inIdx = 0; inIdx < inEdgeData.numLines; inIdx++) {
          for (uint outIdx = 0; outIdx < outEdgeData.numLines; outIdx++) {
            Connection connection;
            connection.vertexNumber = vertexIdx;
            connection.inEdgeNumber = inEdgeNumber;
            connection.outEdgeNumber = outEdgeNumber;
            connection.inLaneNumber = inEdgeNumber + inIdx;
            connection.outLaneNumber = outEdgeNumber + outIdx;
            connection.enabled = false;
            connections.push_back(connection);
            ++connectionsCount;
          }
        }
      }
    }
    intersection.connectionGraphEnd = connectionsCount;

    addTrafficLightScheduleToIntersection(
      intersection,
      vertexIdx,
      trafficLightSchedules,
      connections,
      edgesData);
  }

  std::cout << "\ntrafficLightSchedules" << std::endl;
  for (const auto & i : trafficLightSchedules) {
    std::cout
      << "vertexIdx: " << i.vertexIdx << " "
      << "connectionIdx: " << i.connectionIdx << " "
      << "scheduleGroup: " << i.scheduleGroup << " "
      << "scheduledTime: " << i.scheduledTime << std::endl;
  }

  std::cout << "\nintersections" << std::endl;
  for (const auto & i : updatedIntersections) {
    std::cout
      << "connectionGraphStart: " << i.connectionGraphStart << " "
      << "connectionGraphEnd: " << i.connectionGraphEnd << " "
      << "trafficLightSchedulesStart: " << i.trafficLightSchedulesStart << " "
      << "trafficLightSchedulesEnd: " << i.trafficLightSchedulesEnd << " "
      << "scheduleIdx: " << i.scheduleIdx << " "
      << "currentScheduleGroup: " << i.currentScheduleGroup << " "
      << "timeOfNextUpdate: " << i.timeOfNextUpdate << std::endl;
  }

  std::cout << "\nconnections" << std::endl;
  for (const auto & i : connections) {
    std::cout
      << "inLaneNumber: " << i.inLaneNumber << " "
      << "outLaneNumber: " << i.outLaneNumber << " "
      << "enabled: " << i.enabled << " "
      << "vertexNumber: " << i.vertexNumber << " "
      << "inEdgeNumber: " << i.inEdgeNumber << " "
      << "outEdgeNumber: " << i.outEdgeNumber << std::endl;
  }

  // Instantiate lane map
  laneMap.resize(kMaxMapWidthM * totalLaneMapChunks * 2); // 2: to have two maps.
  memset(laneMap.data(), -1, laneMap.size()*sizeof(unsigned char)); //

  RoadGraph::roadGraphVertexIter_BI vi, viEnd;
  RoadGraph::in_roadGraphEdgeIter_BI Iei, Iei_end;
  RoadGraph::out_roadGraphEdgeIter_BI Oei, Oei_end;
  intersections.resize(boost::num_vertices(inputGraph));
  trafficLights.assign(totalLaneMapChunks, 0);

  for (boost::tie(vi, viEnd) = boost::vertices(inRoadGraph.myRoadGraph_BI); vi != viEnd; ++vi) {
    intersections.at(*vi).state = 0;
    intersections.at(*vi).nextEvent = 0.0f;
    intersections.at(*vi).totalInOutEdges = boost::degree(*vi, inRoadGraph.myRoadGraph_BI);

    if (intersections.at(*vi).totalInOutEdges <= 0) {
      printf("Vertex without in/out edges\n");
      continue;
    }

    if (intersections.at(*vi).totalInOutEdges >= 20) {
      printf("Vertex with more than 20 in/out edges\n");
      continue;
    }

    //sort by angle
    QVector3D referenceVector(0, 1, 0);
    QVector3D p0, p1;
    std::vector<std::pair<LC::RoadGraph::roadGraphEdgeDesc_BI, float>> edgeAngleOut;
    int numOutEdges = 0;

    float angleRef = atan2(referenceVector.y(), referenceVector.x());

    for (boost::tie(Oei, Oei_end) = boost::out_edges(*vi, inRoadGraph.myRoadGraph_BI);
        Oei != Oei_end; ++Oei) {
      if (inRoadGraph.myRoadGraph_BI[*Oei].numberOfLanes == 0) { continue; }
      p0 = inRoadGraph.myRoadGraph_BI[boost::source(*Oei, inRoadGraph.myRoadGraph_BI)].pt;
      p1 = inRoadGraph.myRoadGraph_BI[boost::target(*Oei, inRoadGraph.myRoadGraph_BI)].pt;
      QVector3D edgeDir = (p1 - p0).normalized();
      float angle = angleRef - atan2(edgeDir.y(), edgeDir.x());
      /*LC::RoadGraph::roadGraphVertexDesc_BI sV=boost::source(*Oei, inRoadGraph.myRoadGraph_BI);
      LC::RoadGraph::roadGraphVertexDesc_BI tV=boost::target(*Oei, inRoadGraph.myRoadGraph_BI);
      std::pair<RoadGraph::roadGraphEdgeDesc,bool> edge_pair =
        boost::edge(sV,tV,inRoadGraph.myRoadGraph_BI);*/
      edgeAngleOut.push_back(std::make_pair(*Oei, angle));

      if (edgeDescToLaneMapNum.find(*Oei) == edgeDescToLaneMapNum.end()) {
        printf("->ERROR OUT\n");//edge desc not found in map
      }

      numOutEdges++;
      //edgeAngleOut.push_back(std::make_pair(edge_pair.first,angle));
    }

    std::vector<std::pair<LC::RoadGraph::roadGraphEdgeDesc_BI, float>> edgeAngleIn;
    int numInEdges = 0;

    for (boost::tie(Iei, Iei_end) = boost::in_edges(*vi, inRoadGraph.myRoadGraph_BI);
        Iei != Iei_end; ++Iei) {
      if (inRoadGraph.myRoadGraph_BI[*Iei].numberOfLanes == 0) { continue; }
      p0 = inRoadGraph.myRoadGraph_BI[boost::source(*Iei, inRoadGraph.myRoadGraph_BI)].pt;
      p1 = inRoadGraph.myRoadGraph_BI[boost::target(*Iei, inRoadGraph.myRoadGraph_BI)].pt;
      QVector3D edgeDir = (p0 - p1).normalized();
      float angle = angleRef - atan2(edgeDir.y(), edgeDir.x());
      //std::pair<RoadGraph::roadGraphEdgeDesc, bool> edge_pair = boost::edge(boost::source(*Iei, inRoadGraph.myRoadGraph_BI),boost::target(*Iei, inRoadGraph.myRoadGraph_BI),inRoadGraph.myRoadGraph_BI);
      edgeAngleIn.push_back(std::make_pair(*Iei, angle));

      if (edgeDescToLaneMapNum.find(*Iei) == edgeDescToLaneMapNum.end()) {
        printf("->ERROR IN\n");//edge desc not found in map
        continue;
      }

      numInEdges++;
      //edgeAngleIn.push_back(std::make_pair(edge_pair.first,angle));
    }

    intersections.at(*vi).totalInOutEdges = numOutEdges + numInEdges;

    //save in sorterd way as lane number
    if (edgeAngleOut.size() > 0) {
      std::sort(edgeAngleOut.begin(), edgeAngleOut.end(), compareSecondPartTupleC);
    }

    if (edgeAngleIn.size() > 0) {
      std::sort(edgeAngleIn.begin(), edgeAngleIn.end(), compareSecondPartTupleC);
    }

    //!!!!
    size_t outCount = 0;
    size_t inCount = 0;
    size_t totalCount = 0;

    // Intersection data:
    //  Store the edges that go in or out of this intersection
    //  Said edges will be sorted by angle
    //
    //      0xFF00 0000 Num lines
    //      0x0080 0000 in out (one bit)
    //      0x007F FFFF Edge number
    for (size_t iter = 0; iter < edgeAngleOut.size() + edgeAngleIn.size(); iter++) {
      if ((outCount < edgeAngleOut.size() && inCount < edgeAngleIn.size() &&
           edgeAngleOut[outCount].second <= edgeAngleIn[inCount].second) ||
          (outCount < edgeAngleOut.size() && inCount >= edgeAngleIn.size())) {
        assert(edgeDescToLaneMapNum[edgeAngleOut[outCount].first] < 0x007fffff
            && "Edge number is too high");
        intersections.at(*vi).edge[totalCount] = edgeDescToLaneMapNum[edgeAngleOut[outCount].first];
        intersections.at(*vi).edge[totalCount] |= (edgesData[intersections.at(*vi).edge[totalCount]].numLines << 24); //put the number of lines in each edge
        intersections.at(*vi).edge[totalCount] |= kMaskOutEdge; // 0x000000 mask to define out edge
        outCount++;
      } else {
        assert(edgeDescToLaneMapNum[edgeAngleIn[inCount].first] < 0x007fffff
            && "Edge number is too high");
        intersections.at(*vi).edge[totalCount] = edgeDescToLaneMapNum[edgeAngleIn[inCount].first];
        intersections.at(*vi).edge[totalCount] |= (edgesData[intersections.at(*vi).edge[totalCount]].numLines << 24); //put the number of lines in each edge
        intersections.at(*vi).edge[totalCount] |= kMaskInEdge; // 0x800000 mask to define in edge
        inCount++;
      }

      totalCount++;
    }
  }
}

void B18TrafficLaneMap::resetIntersections(std::vector<B18IntersectionData>
    &intersections, std::vector<uchar> &trafficLights) {
  for (size_t idx = 0; idx < intersections.size(); idx++) {
    intersections[idx].nextEvent = 0.0f; //otherwise they not change until reach time again
    intersections[idx].state = 0; //to make the system to repeat same execution
  }

  if (trafficLights.size() > 0) {
    memset(trafficLights.data(), 0, trafficLights.size()*sizeof(
             uchar));  //to make the system to repeat same execution
  }
}//

}
