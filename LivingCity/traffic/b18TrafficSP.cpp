#include "b18TrafficSP.h"

#include <boost/graph/exterior_property.hpp>
#include <iostream>
#include <fstream>
#include "src/linux_host_memory_logger.h"

#define ROUTE_DEBUG 0
//#define DEBUG_JOHNSON 0

namespace LC {


////////////////
/////////////////////////////
using namespace boost;

inline bool fileExists(const std::string& fileName) {
  std::ifstream f(fileName.c_str());
  return f.good();
}

typedef exterior_vertex_property<RoadGraph::roadBGLGraph_BI, float>
DistanceProperty;
typedef DistanceProperty::matrix_type DistanceMatrix;
typedef DistanceProperty::matrix_map_type DistanceMatrixMap;

void B18TrafficSP::generateRoutesSP(
    LC::RoadGraph::roadBGLGraph_BI &roadGraph,
    std::vector<B18TrafficPerson> &trafficPersonVec,
    std::vector<uint>& indexPathVec,
    std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum,
    int weigthMode, float sample) {
  if (trafficPersonVec.empty()) {
    printf("ERROR generateRoutesSP: people empty");
    return;
  }

  printf(">> generatePathRoutes\n");
  QTime timer;
  timer.start();
  
  uint currIndexPath = 0;
  std::vector<uint> oldIndexPathVec = indexPathVec; // std::move(indexPathVec); // avoid copying
  indexPathVec.clear(); // = std::vector<uint>();

  // 1. Update weight edges
  printf(">> generateRoutes Update weight edges\n");
  int numEdges = 0;
  property_map<RoadGraph::roadBGLGraph_BI, float RoadGraphEdge::*>::type
    weight_pmap = boost::get(&RoadGraphEdge::edge_weight, roadGraph);

  RoadGraph::roadGraphEdgeIter_BI ei, eiEnd;
  float minTravelTime = FLT_MAX;
  float maxTravelTime = -FLT_MAX;
  float minLength = FLT_MAX;
  float maxLength = -FLT_MAX;
  float minSpeed = FLT_MAX;
  float maxSpeed = -FLT_MAX;
  for (boost::tie(ei, eiEnd) = boost::edges(roadGraph); ei != eiEnd; ++ei) {
    numEdges++;
    if ((edgeDescToLaneMapNum.size() > 20) && (numEdges % (edgeDescToLaneMapNum.size() / 20)) == 0) {
      printf("Edge %d of %d (%2.0f%%)\n", numEdges, edgeDescToLaneMapNum.size(), (100.0f * numEdges) / edgeDescToLaneMapNum.size());
    }
    if (roadGraph[*ei].numberOfLanes > 0) {
      float speed;
      if (weigthMode == 0 || roadGraph[*ei].averageSpeed.size() <= 0) {
        speed = roadGraph[*ei].maxSpeedMperSec;//(p0-p1).length();

      } else {
        // Use the avarage speed sampled in a former step
        float avOfAv = 0;
        for (int avOfAvN = 0; avOfAvN < roadGraph[*ei].averageSpeed.size(); avOfAvN++) {
          avOfAv += roadGraph[*ei].averageSpeed[avOfAvN];
        }
        speed = avOfAv / roadGraph[*ei].averageSpeed.size();
      }
      float travelTime = roadGraph[*ei].edgeLength / speed;
      roadGraph[*ei].edge_weight = travelTime;
      weight_pmap[*ei] = travelTime;

      minTravelTime = minTravelTime>travelTime ? travelTime : minTravelTime;
      maxTravelTime = maxTravelTime < travelTime ? travelTime : maxTravelTime;

      minLength = minLength > roadGraph[*ei].edgeLength ? roadGraph[*ei].edgeLength : minLength;
      maxLength = maxLength < roadGraph[*ei].edgeLength ? roadGraph[*ei].edgeLength : maxLength;

      minSpeed = minSpeed > speed ? speed : minSpeed;
      maxSpeed = maxSpeed < speed ? speed : maxSpeed;
      printf("Edge ID: %s, Vertex IDs: %hu %hu, Edge weight: %f\n", roadGraph[*ei].label, roadGraph[*ei].inNum, roadGraph[*ei].outNum, roadGraph[*ei].edge_weight);
    } else {
      roadGraph[*ei].edge_weight =
        100000000.0; //FLT_MAX;// if it has not lines, then the weight is inf
    }
  }
  printf("Travel time Min: %f Max: %f\n", minTravelTime, maxTravelTime);
  printf("Length Min: %f Max: %f\n", minLength, maxLength);
  printf("Speed Min: %f Max: %f\n", minSpeed, maxSpeed);

  //2. Generate route for each person
  printf(">> Printing vertex IDs, edge IDs, and weights\n");
  int numVertex = boost::num_vertices(roadGraph);
  //vertex
  typedef LC::RoadGraph::roadGraphVertexDesc_BI VertexDescriptor;
  typedef std::map<VertexDescriptor, size_t> VertexIndexMap;
  VertexIndexMap mapVertexIndex;
  boost::associative_property_map<VertexIndexMap> pmVertexIndex(mapVertexIndex);

  DistanceMatrix distances(numVertex);
  DistanceMatrixMap dm(distances, roadGraph);



  #ifdef DEBUG_JOHNSON
  std::cerr << "indexPathVec: " << std::endl;
  int i = 0;
  for (const auto x : indexPathVec) {
    std::cerr << i++ << " " << x << " " << std::endl;
  }
  std::cerr << "trafficPersonVec: " << std::endl;
  for (const auto p : trafficPersonVec) {
    std::cerr << p.indexPathInit << " " << p.indexPathCurr << std::endl;
  }
  #endif

}


}  // Closing namespace LC

