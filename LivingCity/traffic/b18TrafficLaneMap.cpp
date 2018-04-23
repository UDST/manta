/************************************************************************************************
*
*		LC Project - B18 Traffic lane map
*
*
*		@desc Class that contains the structure of the lane maps
*		@author igaciad
*
************************************************************************************************/

#include "b18TrafficLaneMap.h"

#define LANE_DEBUG 1

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


void B18TrafficLaneMap::createLaneMap(
  RoadGraph &inRoadGraph,
  std::vector<uchar> &laneMap,
  std::vector<B18EdgeData> &edgesData,
  std::vector<B18IntersectionData> &intersections,
  std::vector<uchar> &trafficLights,
  std::map<uint, RoadGraph::roadGraphEdgeDesc_BI> &laneMapNumToEdgeDesc,
  std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum) {
  //////////////////////////////////////////////////////////
  // GENERATE LANE MAP
  if (LANE_DEBUG) {
    printf("  >> createLaneMap\n");
  }

  // 1. Cretae edgesData and find requires sizes.
  RoadGraph::roadGraphEdgeIter_BI ei, ei_end;
  int edge_count = 0;
  int tNumMapWidth = 0;
  int tNumLanes = 0;
  float maxLength = 0;
  int maxNumLanes = 0;

  //printf("edgesData %d\n",edgesData);
  edgesData.resize(boost::num_edges(inRoadGraph.myRoadGraph_BI) * 4); //4 to make sure it fits

  edgeDescToLaneMapNum.clear();
  laneMapNumToEdgeDesc.clear();

  ////////////////////////////////
  // Check distribution of street length
  float binLength = 1.0f;//1km
  int numBins = 31 / binLength;//maxlength is 26km
  std::vector<int> bins(numBins);
  std::fill(bins.begin(), bins.end(), 0);
  for (boost::tie(ei, ei_end) = boost::edges(inRoadGraph.myRoadGraph_BI); ei != ei_end; ++ei) {
    float length = inRoadGraph.myRoadGraph_BI[*ei].edgeLength;
    int binN = (length / 1000.0f) / binLength;
    //printf("l %.2f binN %d\n", length, binN);
    if (binN < 0 || binN >= numBins) {
      printf("ERROR: Bin out of range: %d of %d\n", binN, numBins);
    }
    bins[binN]++;
  }
  for (int binN = 0; binN < bins.size(); binN++) {
    printf("%.0fkm, %d\n", (binN * binLength+1.0f), bins[binN]);
  }

  /////////////////////////////////
  // Create EdgeData
  // Instead of having maxWidth (for b18 would have been 26km), we define a max width for the map and we wrap down.
  for (boost::tie(ei, ei_end) = boost::edges(inRoadGraph.myRoadGraph_BI); ei != ei_end; ++ei) {

    // num lanes
    int numLanes = inRoadGraph.myRoadGraph_BI[*ei].numberOfLanes;

    if (numLanes == 0) {
      continue;//edges with zero lines just skip
    }

    // length
    edgesData[tNumMapWidth].length = inRoadGraph.myRoadGraph_BI[*ei].edgeLength;
    edgesData[tNumMapWidth].maxSpeedMperSec = inRoadGraph.myRoadGraph_BI[*ei].maxSpeedMperSec;

    if (maxLength < edgesData[tNumMapWidth].length) {
      maxLength = edgesData[tNumMapWidth].length;
    }
    if (maxNumLanes < numLanes) {
      maxNumLanes = numLanes;
    }

    int numWidthNeeded = ceil(edgesData[tNumMapWidth].length / kMaxMapWidthM); // number of width needed if > than kMaxMapWidthM
    edgesData[tNumMapWidth].numLines = numLanes * numWidthNeeded;
    // next intersection
    edgesData[tNumMapWidth].nextInters = boost::target(*ei,
                                      inRoadGraph.myRoadGraph_BI);

    edgeDescToLaneMapNum.insert(std::make_pair(*ei, tNumMapWidth));
    laneMapNumToEdgeDesc.insert(std::make_pair(tNumMapWidth, *ei));

    tNumMapWidth += numLanes * numWidthNeeded;
    tNumLanes += numLanes;
    edge_count++;
  }

  edgesData.resize(tNumMapWidth);

  if (LANE_DEBUG) {
    printf("Num edges %d Num Lanes %d Max Leng %f Max num lanes %d\n", edge_count, tNumMapWidth, maxLength, maxNumLanes);
  }

  // 2. RESIZE LANE MAP

  printf("Total Memory %d\n", kMaxMapWidthM * tNumMapWidth * 2);
  laneMap.resize(kMaxMapWidthM * tNumMapWidth * 2); // 2: to have two maps.
  memset(laneMap.data(), -1, laneMap.size()*sizeof(unsigned char)); //

  //printf("aa\n");
  //////////////////////////////////////////////////////////
  // GENERATE INTERSECTION INFO
  printf("Start intersection info\n");
  RoadGraph::roadGraphVertexIter_BI vi, viEnd;
  RoadGraph::in_roadGraphEdgeIter_BI Iei, Iei_end;
  RoadGraph::out_roadGraphEdgeIter_BI Oei, Oei_end;
  intersections.resize(boost::num_vertices(inRoadGraph.myRoadGraph_BI));//as many as vertices
  trafficLights.resize(tNumLanes);
  memset(trafficLights.data(), 0, trafficLights.size()*sizeof(uchar));

  if (LANE_DEBUG) {
    printf(">>Generate Intersection Info\n");
  }

  for (boost::tie(vi, viEnd) = boost::vertices(inRoadGraph.myRoadGraph_BI); vi != viEnd; ++vi) {
    intersections[*vi].state = 0;
    intersections[*vi].nextEvent = 0.0f;
    
    intersections[*vi].totalInOutEdges = boost::degree(*vi, inRoadGraph.myRoadGraph_BI);

    if (intersections[*vi].totalInOutEdges <= 0) {
      //printf("Vertex without in/out edges\n");
      continue;
    }

    if (intersections[*vi].totalInOutEdges >= 20) {
      printf("Vertex with more than 20 in/out edges\n");
      continue;
    }

    //printf("Total: %d\n",intersections[*vi].totalInOutEdges);
    //sort by angle
    QVector3D referenceVector(0, 1, 0);
    QVector3D p0, p1;
    std::vector<std::pair<LC::RoadGraph::roadGraphEdgeDesc_BI, float>> edgeAngleOut;
    //printf("Out\n");
    int numOutEdges = 0;

    float angleRef = atan2(referenceVector.y(), referenceVector.x());

    for (boost::tie(Oei, Oei_end) = boost::out_edges(*vi,
                                    inRoadGraph.myRoadGraph_BI); Oei != Oei_end; ++Oei) {
      if (inRoadGraph.myRoadGraph_BI[*Oei].numberOfLanes == 0) {
        continue;
      }

      p0 = inRoadGraph.myRoadGraph_BI[boost::source(*Oei,
                                      inRoadGraph.myRoadGraph_BI)].pt;
      p1 = inRoadGraph.myRoadGraph_BI[boost::target(*Oei,
                                      inRoadGraph.myRoadGraph_BI)].pt;
      QVector3D edgeDir = (p1 - p0).normalized(); // NOTE p1-p0
      //float angle=QVector3D::dotProduct(referenceVector,edgeDir);
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

    //printf("Out %d\n",numOutEdges);
    std::vector<std::pair<LC::RoadGraph::roadGraphEdgeDesc_BI, float>> edgeAngleIn;
    //printf("In\n");
    int numInEdges = 0;

    for (boost::tie(Iei, Iei_end) = boost::in_edges(*vi,
                                    inRoadGraph.myRoadGraph_BI); Iei != Iei_end; ++Iei) {
      if (inRoadGraph.myRoadGraph_BI[*Iei].numberOfLanes == 0) {
        continue;
      }

      p0 = inRoadGraph.myRoadGraph_BI[boost::source(*Iei,
                                      inRoadGraph.myRoadGraph_BI)].pt;
      p1 = inRoadGraph.myRoadGraph_BI[boost::target(*Iei,
                                      inRoadGraph.myRoadGraph_BI)].pt;
      QVector3D edgeDir = (p0 - p1).normalized(); // NOTE p0-p1
      //float angle=QVector3D::dotProduct(referenceVector,edgeDir);
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

    intersections[*vi].totalInOutEdges = numOutEdges + numInEdges;

    //printf("In %d\n",numInEdges);
    //printf("Sort\n");
    //save in sorterd way as lane number
    if (edgeAngleOut.size() > 0) {
      std::sort(edgeAngleOut.begin(), edgeAngleOut.end(), compareSecondPartTupleC);
    }

    if (edgeAngleIn.size() > 0) {
      std::sort(edgeAngleIn.begin(), edgeAngleIn.end(), compareSecondPartTupleC);
    }

    //!!!!
    int outCount = 0;
    int inCount = 0;
    int totalCount = 0;

    //printf("count %d\n",);
    // INTERSECTION
    //      0xFF00 0000 Num lines
    //      0x0080 0000 in out (one bit)
    //      0x007F FFFF Edge number
    for (int iter = 0; iter < edgeAngleOut.size() + edgeAngleIn.size(); iter++) {
      if ((outCount < edgeAngleOut.size() && inCount < edgeAngleIn.size() && 
           edgeAngleOut[outCount].second <= edgeAngleIn[inCount].second) ||
          (outCount < edgeAngleOut.size() && inCount >= edgeAngleIn.size())) {
        intersections[*vi].edge[totalCount] = edgeDescToLaneMapNum[edgeAngleOut[outCount].first];
        intersections[*vi].edge[totalCount] |= (edgesData[intersections[*vi].edge[totalCount]].numLines << 24); //put the number of lines in each edge
        intersections[*vi].edge[totalCount] |= kMaskOutEdge; // 0x000000 mask to define out edge
        outCount++;
      } else {
        intersections[*vi].edge[totalCount] = edgeDescToLaneMapNum[edgeAngleIn[inCount].first];
        intersections[*vi].edge[totalCount] |= (edgesData[intersections[*vi].edge[totalCount]].numLines << 24); //put the number of lines in each edge
        intersections[*vi].edge[totalCount] |= kMaskInEdge; // 0x800000 mask to define in edge
        inCount++;
      }

      totalCount++;
    }

    if (totalCount != intersections[*vi].totalInOutEdges) {
      printf("Error totalCount!=intersections[*vi].totalInOutEdges %d %d\n",
             totalCount, intersections[*vi].totalInOutEdges);
    }

  }

  if (LANE_DEBUG) {
    printf("<<Generate Intersection Info\n");
  }

  if (LANE_DEBUG) {
    printf("  << createLaneMap\n");
  }
}//

void B18TrafficLaneMap::resetIntersections(std::vector<B18IntersectionData>
    &intersections, std::vector<uchar> &trafficLights) {
  for (int i = 0; i < intersections.size(); i++) {
    intersections[i].nextEvent =
      0.0f; //otherwise they not change until reach time again
    intersections[i].state = 0; //to make the system to repeat same execution
  }

  if (trafficLights.size() > 0) {
    memset(trafficLights.data(), 0, trafficLights.size()*sizeof(
             uchar));  //to make the system to repeat same execution
  }
}//

}
