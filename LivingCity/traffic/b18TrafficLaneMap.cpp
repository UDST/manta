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


// TODO: Add traffic light schedules
void B18TrafficLaneMap::createLaneMap(
    const RoadGraph &inRoadGraph,
    std::vector<uchar> &laneMap,
    std::vector<B18EdgeData> &edgesData,
    std::vector<B18IntersectionData> &intersections,
    std::vector<uchar> &trafficLights,
    std::map<uint, RoadGraph::roadGraphEdgeDesc_BI> &laneMapNumToEdgeDesc,
    std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum,
    std::vector<LC::Connection> &connections,
    std::vector<LC::Intersection> &updatedIntersections) {
  edgesData.resize(boost::num_edges(inRoadGraph.myRoadGraph_BI) * 4);  //4 to make sure it fits

  edgeDescToLaneMapNum.clear();
  laneMapNumToEdgeDesc.clear();
  connections.clear();
  updatedIntersections.clear();

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
  auto vertices_begin = p.first;
  const auto vertices_end = p.second;
  int connectionsIndex = 0;
  // TODO: Add traffic light schedules
  for (auto vertices_it = vertices_begin; vertices_it != vertices_end; ++vertices_it) {
    auto in_edges_pair = boost::in_edges(*vertices_it, inputGraph);
    auto in_edges_begin = in_edges_pair.first;
    const auto in_edges_end = in_edges_pair.second;
    const auto vertexNumber = *vertices_it;
    updatedIntersections.at(vertexNumber).connectionGraphStart = connectionsIndex;
    for (auto in_edges_it = in_edges_begin; in_edges_it != in_edges_end; ++in_edges_it) {
      const auto & inEdgeNumber = edgeDescToLaneMapNum.at(*in_edges_it);
      const auto & inEdgeData = edgesData[inEdgeNumber];
      auto p2 = boost::out_edges(*vertices_it, inputGraph);
      auto outEdgesBegin = p2.first;
      const auto outEdgesEnd = p2.second;
      for (auto outEdgesIt = outEdgesBegin; outEdgesIt != outEdgesEnd; ++outEdgesIt) {
        const auto & outEdgeNumber = edgeDescToLaneMapNum.at(*outEdgesIt);
        const auto & outEdgeData = edgesData[outEdgeNumber];
        if (inEdgeData.originalSourceVertexIndex == outEdgeData.originalTargetVertexIndex) {
          // Avoid U-turns
          continue;
        }

        assert(inEdgeData.originalTargetVertexIndex == outEdgeData.originalSourceVertexIndex);
        for (int inIdx = 0; inIdx < inEdgeData.numLines; inIdx++) {
          for (int outIdx = 0; outIdx < outEdgeData.numLines; outIdx++) {
            Connection connection;
            connection.vertexNumber = vertexNumber;
            connection.inEdgeNumber = inEdgeNumber;
            connection.outEdgeNumber = outEdgeNumber;
            connection.inLaneNumber = inEdgeNumber + inIdx;
            connection.outLaneNumber = outEdgeNumber + outIdx;
            // TODO: Set this default value to false once the traffic light schedules are ready
            connection.enabled = true;
            connections.push_back(connection);
            ++connectionsIndex;
          }
        }
      }
    }
    updatedIntersections.at(vertexNumber).connectionGraphEnd = connectionsIndex;
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
  // TODO: Check if this needs to be modified
  for (size_t idx = 0; idx < intersections.size(); idx++) {
    intersections[idx].nextEvent =
      0.0f; //otherwise they not change until reach time again
    intersections[idx].state = 0; //to make the system to repeat same execution
  }

  if (trafficLights.size() > 0) {
    memset(trafficLights.data(), 0, trafficLights.size()*sizeof(
             uchar));  //to make the system to repeat same execution
  }
}//

}
