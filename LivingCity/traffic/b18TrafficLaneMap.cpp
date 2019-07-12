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
#include "b18TrafficLaneMap.h"
#include "sp/graph.h"

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
  bool compareSecondPartTupleCNew(const std::pair<std::shared_ptr<abm::Graph::Edge>, float> &i,
		  	       const std::pair<std::shared_ptr<abm::Graph::Edge>, float> &j) {
    return (i.second < j.second);
  }

	/*
  bool compareSecondPartTupleCNew(float i,
    float j) {
    return (i < j);
  }
*/
}

void B18TrafficLaneMap::createLaneMapSP(const std::shared_ptr<abm::Graph>& graph_, std::vector<uchar> &laneMap,
      std::vector<B18EdgeData> &edgesData, std::vector<B18IntersectionData> &intersections,
      std::vector<uchar> &trafficLights, 
      //std::map<uint,
      //RoadGraph::roadGraphEdgeDesc_BI> &laneMapNumToEdgeDesc,
      //std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum) {
      std::map<uint, std::shared_ptr<abm::Graph::Edge>> &laneMapNumToEdgeDescSP,
      std::map<std::shared_ptr<abm::Graph::Edge>, uint> &edgeDescToLaneMapNumSP) {
      //std::tuple<graph::vertex_t, graph::vertex_t>, std::shared_ptr<Edge>> &laneMapNumToEdgeDesc,
      //std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum) {
  	printf("edgesData size = %d\n", edgesData.size());
	/* FOR DATA ACCESS SANITY
	for (auto const& x : graph_->edges_) {
		//std::cout << "vertex = " << std::get<1>(std::get<0>(x)) << "\n"; //gets the second vertex value of edge
		//std::cout << "weight = " << std::get<1>(x)->second[0]<< "\n"; //gets weight of edge
		std::cout << "weight = " << std::get<1>(x)->second[1]<< "\n"; //gets weight of edge
	}
	*/
  // GENERATE LANE MAP
  if (LANE_DEBUG) {
    printf("  >> createLaneMap\n");
  }

  // 1. Cretae edgesData and find requires sizes.
  RoadGraph::roadGraphEdgeIter_BI ei, ei_end;
  int edge_count = 0;
  int tNumMapWidth = 0;
  int tNumLanes = 0;
  
  edgesData.resize(graph_->nedges() * 4); //4 to make sure it fits
  
  edgeDescToLaneMapNumSP.clear();
  laneMapNumToEdgeDescSP.clear();

  // Check distribution of street length
  float binLength = 1.0f;//1km
  int numBins = 31 / binLength;//maxlength is 26km
  std::vector<int> bins(numBins, 0);
  for (auto const& x : graph_->edges_) {
	  const float metersLength = std::get<1>(x)->second[0];
	  const int binN = (metersLength / 1000.0f) / binLength;
	  assert(0 <= binN && binN < numBins && "Edge over max length");
	  bins[binN]++;
  }
  printf("bin size = %d\n", bins.size());
  for (int binN = 0; binN < bins.size(); binN++) {
    printf("%.0fkm, %d\n", (binN * binLength+1.0f), bins[binN]);
  }

  /////////////////////////////////
  // Create EdgeData
  // Instead of having maxWidth (for b18 would have been 26km), we define a max width for the map and we wrap down.
  float maxLength = 0;
  int maxNumLanes = 0;
  for (auto const& x : graph_->edges_) {
	  const int numLanes = std::get<1>(x)->second[1];

    if (numLanes == 0) { continue; }

    edgesData[tNumMapWidth].length = std::get<1>(x)->second[0];
    edgesData[tNumMapWidth].maxSpeedMperSec = std::get<1>(x)->second[2];

    if (maxLength < edgesData[tNumMapWidth].length) { maxLength = edgesData[tNumMapWidth].length; }
    if (maxNumLanes < numLanes) { maxNumLanes = numLanes; }

    const int numWidthNeeded = ceil(edgesData[tNumMapWidth].length / kMaxMapWidthM);
    edgesData[tNumMapWidth].numLines = numLanes;
    edgesData[tNumMapWidth].nextInters = std::get<1>(std::get<0>(x));

    //edgeDescToLaneMapNum.insert(std::make_pair(*ei, tNumMapWidth));
    //laneMapNumToEdgeDesc.insert(std::make_pair(tNumMapWidth, *ei));
    edgeDescToLaneMapNumSP.insert(std::make_pair(std::get<1>(x), tNumMapWidth));
    laneMapNumToEdgeDescSP.insert(std::make_pair(tNumMapWidth, std::get<1>(x)));

    tNumMapWidth += numLanes * numWidthNeeded;
    tNumLanes += numLanes;
    edge_count++;
  }

  edgesData.resize(tNumMapWidth);

  if (LANE_DEBUG) {
    printf("Num edges %d Num Lanes %d Num Lanes Width %d Max Leng %f Max num lanes %d\n",
        edge_count, tNumLanes, tNumMapWidth, maxLength, maxNumLanes);
  }

  // 2. RESIZE LANE MAP

  printf("Total Memory %d\n", kMaxMapWidthM * tNumMapWidth * 2);
  laneMap.resize(kMaxMapWidthM * tNumMapWidth * 2); // 2: to have two maps.
  memset(laneMap.data(), -1, laneMap.size()*sizeof(unsigned char)); //

  //////////////////////////////////////////////////////////
  // GENERATE INTERSECTION INFO
  printf("Start intersection info\n");
  RoadGraph::roadGraphVertexIter_BI vi, viEnd;
  RoadGraph::in_roadGraphEdgeIter_BI Iei, Iei_end;
  RoadGraph::out_roadGraphEdgeIter_BI Oei, Oei_end;
  //intersections.resize(boost::num_vertices(inRoadGraph.myRoadGraph_BI));//as many as vertices
  intersections.resize(graph_->nvertices_);//as many as vertices
  trafficLights.assign(tNumMapWidth, 0);
  //trafficLights.resize(tNumMapWidth); // we could use tNumLanes but then the edge number would not match and we would need to add logic.
  //memset(trafficLights.data(), 0, trafficLights.size()*sizeof(uchar));

  //for (boost::tie(vi, viEnd) = boost::vertices(inRoadGraph.myRoadGraph_BI); vi != viEnd; ++vi) {}
  int index = 0;
  for (const auto& vertex : graph_->vertex_edges_) {
    intersections[std::get<0>(vertex)].state = 0;
    intersections[std::get<0>(vertex)].nextEvent = 0.0f;
    intersections[std::get<0>(vertex)].totalInOutEdges = vertex.second.size();
    if (intersections[std::get<0>(vertex)].totalInOutEdges <= 0) {
      printf("Vertex without in/out edges\n");
      continue;
    }

  if (intersections[std::get<0>(vertex)].totalInOutEdges >= 20) {
      printf("Vertex with more than 20 in/out edges\n");
      continue;
    }
    index++;

    
  //sort by angle
    QVector3D referenceVector(0, 1, 0);
    QVector3D p0, p1;
    //std::vector<std::pair<LC::RoadGraph::roadGraphEdgeDesc_BI, float>> edgeAngleOut;
    std::vector<std::pair<std::shared_ptr<abm::Graph::Edge>, float>> edgeAngleOut;
    int numOutEdges = 0;

    float angleRef = atan2(referenceVector.y(), referenceVector.x());

    for (const auto& edge : graph_->vertex_out_edges_[std::get<0>(vertex)]) {
      //if (inRoadGraph.myRoadGraph_BI[*Oei].numberOfLanes == 0) { continue; }
      if (edge->second[1] == 0) { continue; }

      p0 = graph_->vertices_data_[std::get<0>(vertex)];
      p1 = graph_->vertices_data_[edge->first.second];

      QVector3D edgeDir = (p1 - p0).normalized();
      float angle = angleRef - atan2(edgeDir.y(), edgeDir.x());

      //LC::RoadGraph::roadGraphVertexDesc_BI sV=boost::source(*Oei, inRoadGraph.myRoadGraph_BI);
      //LC::RoadGraph::roadGraphVertexDesc_BI tV=boost::target(*Oei, inRoadGraph.myRoadGraph_BI);
      //LC::RoadGraph::roadGraphVertexDesc_BI sV = std::get<0>(std::get<0>(x));
      //LC::RoadGraph::roadGraphVertexDesc_BI tV = std::get<1>(std::get<0>(x));
      //std::pair<RoadGraph::roadGraphEdgeDesc,bool> edge_pair =
      //  boost::edge(sV,tV,inRoadGraph.myRoadGraph_BI);
      //edgeAngleOut.push_back(std::make_pair(*Oei, angle));
      edgeAngleOut.push_back(std::make_pair(edge, angle));

      if (edgeDescToLaneMapNumSP.find(edge) == edgeDescToLaneMapNumSP.end()) {
        printf("->ERROR OUT\n");//edge desc not found in map
      }

      numOutEdges++;
      //edgeAngleOut.push_back(std::make_pair(edge_pair.first,angle));
    }
    
  printf("Out %d\n",numOutEdges);
   
  //std::vector<std::pair<LC::RoadGraph::roadGraphEdgeDesc_BI, float>> edgeAngleIn;
    //std::vector<std::pair<abm::graph::vertex_t, std::vector<std::shared_ptr<Edge>>, float>> edgeAngleIn;
    std::vector<std::pair<std::shared_ptr<abm::Graph::Edge>, float>> edgeAngleIn;
    int numInEdges = 0;

    //for (boost::tie(Iei, Iei_end) = boost::in_edges(*vi, inRoadGraph.myRoadGraph_BI);
    //    Iei != Iei_end; ++Iei) {}
    for (const auto& edge : graph_->vertex_in_edges_[std::get<0>(vertex)]) {
      //if (inRoadGraph.myRoadGraph_BI[*Iei].numberOfLanes == 0) { continue; }
      if (edge->second[1] == 0) { continue; }

      p0 = graph_->vertices_data_[std::get<0>(vertex)];
      p1 = graph_->vertices_data_[edge->first.second];

      QVector3D edgeDir = (p0 - p1).normalized();
      float angle = angleRef - atan2(edgeDir.y(), edgeDir.x());
      
      //edgeAngleIn.push_back(std::make_pair(*Iei, angle));
      edgeAngleIn.push_back(std::make_pair(edge, angle));

      
      if (edgeDescToLaneMapNumSP.find(edge) == edgeDescToLaneMapNumSP.end()) {
        printf("->ERROR IN\n");//edge desc not found in map
        continue;
      }

      numInEdges++;
      //edgeAngleIn.push_back(std::make_pair(edge_pair.first,angle));
    }
    printf("In %d\n",numInEdges);

    intersections[std::get<0>(vertex)].totalInOutEdges = numOutEdges + numInEdges;
    //save in sorterd way as lane number
    if (edgeAngleOut.size() > 0) {
      std::sort(edgeAngleOut.begin(), edgeAngleOut.end(), compareSecondPartTupleCNew);
    }

    if (edgeAngleIn.size() > 0) {
      std::sort(edgeAngleIn.begin(), edgeAngleIn.end(), compareSecondPartTupleCNew);
    }
    //!!!!
    int outCount = 0;
    int inCount = 0;
    int totalCount = 0;

    // Use intersections[std::get<0>(vertex)] = blah blah to set those values for each vertex
    // Intersection data:
    //  Store the edges that go in or out of this intersection
    //  Said edges will be sorted by angle
    //  
    //      0xFF00 0000 Num lines
    //      0x0080 0000 in out (one bit)
    //      0x007F FFFF Edge number
    for (int iter = 0; iter < edgeAngleOut.size() + edgeAngleIn.size(); iter++) {
      if ((outCount < edgeAngleOut.size() && inCount < edgeAngleIn.size() && 
           edgeAngleOut[outCount] <= edgeAngleIn[inCount]) ||
          (outCount < edgeAngleOut.size() && inCount >= edgeAngleIn.size())) {
        assert(edgeDescToLaneMapNumSP[edgeAngleOut[outCount].first] < 0x007fffff && "Edge number is too high");
        intersections[std::get<0>(vertex)].edge[totalCount] = edgeDescToLaneMapNumSP[edgeAngleOut[outCount].first];
        intersections[std::get<0>(vertex)].edge[totalCount] |= (edgesData[intersections[std::get<0>(vertex)].edge[totalCount]].numLines << 24); //put the number of lines in each edge
        intersections[std::get<0>(vertex)].edge[totalCount] |= kMaskOutEdge; // 0x000000 mask to define out edge
        outCount++;
      } else {
        assert(edgeDescToLaneMapNumSP[edgeAngleIn[inCount].first] < 0x007fffff && "Edge number is too high");
        intersections[std::get<0>(vertex)].edge[totalCount] = edgeDescToLaneMapNumSP[edgeAngleIn[inCount].first];
        intersections[std::get<0>(vertex)].edge[totalCount] |= (edgesData[intersections[std::get<0>(vertex)].edge[totalCount]].numLines << 24); //put the number of lines in each edge
        intersections[std::get<0>(vertex)].edge[totalCount] |= kMaskInEdge; // 0x800000 mask to define in edge
        inCount++;
      }

      totalCount++;
      std::cout << "intersections data " << intersections[std::get<0>(vertex)].edge[totalCount] << "\n";
    }


    if (totalCount != intersections[*vi].totalInOutEdges) {
      printf("Error totalCount!=intersections[std::get<0>(vertex)].totalInOutEdges %d %d\n",
             totalCount, intersections[*vi].totalInOutEdges);
    }
  }
}

void B18TrafficLaneMap::createLaneMap(
    const RoadGraph &inRoadGraph,
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

  edgesData.resize(boost::num_edges(inRoadGraph.myRoadGraph_BI) * 4); //4 to make sure it fits

  edgeDescToLaneMapNum.clear();
  laneMapNumToEdgeDesc.clear();

  ////////////////////////////////
  // Check distribution of street length
  float binLength = 1.0f;//1km
  int numBins = 31 / binLength;//maxlength is 26km
  std::vector<int> bins(numBins, 0);
  for (boost::tie(ei, ei_end) = boost::edges(inRoadGraph.myRoadGraph_BI); ei != ei_end; ++ei) {
    const float metersLength = inRoadGraph.myRoadGraph_BI[*ei].edgeLength;
    const int binN = (metersLength / 1000.0f) / binLength;
    assert(0 <= binN && binN < numBins && "Edge over max length");
    bins[binN]++;
  }
  for (int binN = 0; binN < bins.size(); binN++) {
    printf("%.0fkm, %d\n", (binN * binLength+1.0f), bins[binN]);
  }

  /////////////////////////////////
  // Create EdgeData
  // Instead of having maxWidth (for b18 would have been 26km), we define a max width for the map and we wrap down.
  float maxLength = 0;
  int maxNumLanes = 0;
  for (boost::tie(ei, ei_end) = boost::edges(inRoadGraph.myRoadGraph_BI); ei != ei_end; ++ei) {
    const int numLanes = inRoadGraph.myRoadGraph_BI[*ei].numberOfLanes;

    if (numLanes == 0) { continue; }

    edgesData[tNumMapWidth].length = inRoadGraph.myRoadGraph_BI[*ei].edgeLength;
    edgesData[tNumMapWidth].maxSpeedMperSec = inRoadGraph.myRoadGraph_BI[*ei].maxSpeedMperSec;

    if (maxLength < edgesData[tNumMapWidth].length) { maxLength = edgesData[tNumMapWidth].length; }
    if (maxNumLanes < numLanes) { maxNumLanes = numLanes; }

    const int numWidthNeeded = ceil(edgesData[tNumMapWidth].length / kMaxMapWidthM);
    edgesData[tNumMapWidth].numLines = numLanes;
    edgesData[tNumMapWidth].nextInters = boost::target(*ei, inRoadGraph.myRoadGraph_BI);

    edgeDescToLaneMapNum.insert(std::make_pair(*ei, tNumMapWidth));
    laneMapNumToEdgeDesc.insert(std::make_pair(tNumMapWidth, *ei));

    tNumMapWidth += numLanes * numWidthNeeded;
    tNumLanes += numLanes;
    edge_count++;
  }

  edgesData.resize(tNumMapWidth);

  if (LANE_DEBUG) {
    printf("Num edges %d Num Lanes %d Num Lanes Width %d Max Leng %f Max num lanes %d\n",
        edge_count, tNumLanes, tNumMapWidth, maxLength, maxNumLanes);
  }

  // 2. RESIZE LANE MAP

  printf("Total Memory %d\n", kMaxMapWidthM * tNumMapWidth * 2);
  laneMap.resize(kMaxMapWidthM * tNumMapWidth * 2); // 2: to have two maps.
  memset(laneMap.data(), -1, laneMap.size()*sizeof(unsigned char)); //

  //////////////////////////////////////////////////////////
  // GENERATE INTERSECTION INFO
  printf("Start intersection info\n");
  RoadGraph::roadGraphVertexIter_BI vi, viEnd;
  RoadGraph::in_roadGraphEdgeIter_BI Iei, Iei_end;
  RoadGraph::out_roadGraphEdgeIter_BI Oei, Oei_end;
  intersections.resize(boost::num_vertices(inRoadGraph.myRoadGraph_BI));//as many as vertices
  trafficLights.assign(tNumMapWidth, 0);
  //trafficLights.resize(tNumMapWidth); // we could use tNumLanes but then the edge number would not match and we would need to add logic.
  //memset(trafficLights.data(), 0, trafficLights.size()*sizeof(uchar));

  for (boost::tie(vi, viEnd) = boost::vertices(inRoadGraph.myRoadGraph_BI); vi != viEnd; ++vi) {
    intersections[*vi].state = 0;
    intersections[*vi].nextEvent = 0.0f;
    intersections[*vi].totalInOutEdges = boost::degree(*vi, inRoadGraph.myRoadGraph_BI);
    if (intersections[*vi].totalInOutEdges <= 0) {
      printf("Vertex without in/out edges\n");
      continue;
    }

    if (intersections[*vi].totalInOutEdges >= 20) {
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

    //printf("Out %d\n",numOutEdges);
    std::vector<std::pair<LC::RoadGraph::roadGraphEdgeDesc_BI, float>> edgeAngleIn;
    //printf("In\n");
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

    intersections[*vi].totalInOutEdges = numOutEdges + numInEdges;

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

    // Intersection data:
    //  Store the edges that go in or out of this intersection
    //  Said edges will be sorted by angle
    //  
    //      0xFF00 0000 Num lines
    //      0x0080 0000 in out (one bit)
    //      0x007F FFFF Edge number
    for (int iter = 0; iter < edgeAngleOut.size() + edgeAngleIn.size(); iter++) {
      if ((outCount < edgeAngleOut.size() && inCount < edgeAngleIn.size() && 
           edgeAngleOut[outCount].second <= edgeAngleIn[inCount].second) ||
          (outCount < edgeAngleOut.size() && inCount >= edgeAngleIn.size())) {
        assert(edgeDescToLaneMapNum[edgeAngleOut[outCount].first] < 0x007fffff && "Edge number is too high");
        intersections[*vi].edge[totalCount] = edgeDescToLaneMapNum[edgeAngleOut[outCount].first];
        intersections[*vi].edge[totalCount] |= (edgesData[intersections[*vi].edge[totalCount]].numLines << 24); //put the number of lines in each edge
        intersections[*vi].edge[totalCount] |= kMaskOutEdge; // 0x000000 mask to define out edge
        outCount++;
      } else {
        assert(edgeDescToLaneMapNum[edgeAngleIn[inCount].first] < 0x007fffff && "Edge number is too high");
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
