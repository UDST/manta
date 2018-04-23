
#include "b18TrafficJohnson.h"

#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/exterior_property.hpp>

#define ROUTE_DEBUG 0

namespace LC {


////////////////
/////////////////////////////
using namespace boost;

typedef exterior_vertex_property<RoadGraph::roadBGLGraph_BI, float>
DistanceProperty;
typedef DistanceProperty::matrix_type DistanceMatrix;
typedef DistanceProperty::matrix_map_type DistanceMatrixMap;

void B18TrafficJohnson::generateRoutes(
  LC::RoadGraph::roadBGLGraph_BI &roadGraph,
  std::vector<B18TrafficPerson> &trafficPersonVec,
  std::vector<uint>& indexPathVec,
  std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum,
  int weigthMode, float sample) {

  if (trafficPersonVec.empty()) {
    printf("ERROR generateRoutes: people empty");
    return;
  }

  printf(">> generatePathRoutes\n");
  QTime timer;
  timer.start();
  
  uint currIndexPath = 0; // counter to keep track where to put more
  indexPathVec.clear();
  indexPathVec.resize(trafficPersonVec.size() * 140); // initial allocation (so we do not add)
  indexPathVec[currIndexPath++] = -1; // first path is empty

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
    } else {
      roadGraph[*ei].edge_weight =
        100000000.0; //FLT_MAX;// if it has not lines, then the weight is inf
    }
  }
  printf("Travel time Min: %f Max: %f\n", minTravelTime, maxTravelTime);
  printf("Length Min: %f Max: %f\n", minLength, maxLength);
  printf("Speed Min: %f Max: %f\n", minSpeed, maxSpeed);

  //2. Generate route for each person
  printf(">> generateRoutesGenerate\n");
  int numVertex = boost::num_vertices(roadGraph);
  //vertex
  typedef LC::RoadGraph::roadGraphVertexDesc_BI VertexDescriptor;
  typedef std::map<VertexDescriptor, size_t> VertexIndexMap;
  VertexIndexMap mapVertexIndex;
  boost::associative_property_map<VertexIndexMap> pmVertexIndex(mapVertexIndex);

  DistanceMatrix distances(numVertex);
  DistanceMatrixMap dm(distances, roadGraph);


  ////////////////////////
  // CALL JOHNSON
  printf("Call Johnson\n");
  boost::johnson_all_pairs_shortest_paths(roadGraph, dm, weight_map(weight_pmap));
  
  // check maxDist
  printf("numVertex %d\n", numVertex);
  float maxDist = -1.0f;
  for (int vN = 0; vN < numVertex; vN++) {
    for (int vN2 = 0; vN2 < numVertex; vN2++) {
      maxDist = maxDist < dm[vN][vN2] ? dm[vN][vN2] : maxDist;
    }
    //printf("\n");
  }

  printf("maxDist %f\n", maxDist);
  ////////////////////////
  // Create routes
  uint noAccesible = 0;
  uint sameSrcDst = 0;
  QTime timer2;
  timer2.start();
  std::vector<int> pathHistogram(250, 0); // this should be around 140.

  for (int p = 0; p < trafficPersonVec.size(); p++) {
    if (trafficPersonVec.size() > 200) {
      if ((p % (trafficPersonVec.size() / 20)) == 0) {
        printf("Route %d of %d (%2.0f%%)\n", p, trafficPersonVec.size(), (100.0f * p) / trafficPersonVec.size());
      }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    trafficPersonVec[p].indexPathInit = currIndexPath;

    LC::RoadGraph::roadGraphVertexDesc_BI srcvertex = trafficPersonVec[p].init_intersection;
    LC::RoadGraph::roadGraphVertexDesc_BI tgtvertex = trafficPersonVec[p].end_intersection;

    // check whether source same than target (we have arrived)
    if (tgtvertex == srcvertex) {
      trafficPersonVec[p].indexPathInit = 0; // that index points to -1
      sameSrcDst++;
      continue;
    }

    // check if accesible
    if (dm[srcvertex][tgtvertex] == (std::numeric_limits < float >::max)()) {
      trafficPersonVec[p].indexPathInit = 0; // that index points to -1
      noAccesible++;
      continue;
    }

    // find path
    LC::RoadGraph::roadGraphVertexDesc_BI currvertex = srcvertex;//init
    RoadGraph::out_roadGraphEdgeIter_BI Oei, Oei_end;
    LC::RoadGraph::roadGraphVertexDesc_BI srcPosEdgeV, tgtPosEdgeV;

    int currIndex = 0;

    while (currvertex != tgtvertex) {
      
      // check all outedges from currvertex which one continues to shortest path
      bool cont = false;

      for (boost::tie(Oei, Oei_end) = boost::out_edges(currvertex, roadGraph); Oei != Oei_end; ++Oei) {

        srcPosEdgeV = boost::source(*Oei, roadGraph);
        tgtPosEdgeV = boost::target(*Oei, roadGraph);

        if (std::abs<float>(dm[currvertex][tgtPosEdgeV] + dm[tgtPosEdgeV][tgtvertex] - dm[currvertex][tgtvertex]) < 0.1f) {  // link found
          std::pair<RoadGraph::roadGraphEdgeDesc_BI, bool> edge_pair = boost::edge(srcPosEdgeV, tgtPosEdgeV, roadGraph);

          if (edge_pair.second == false) {
            printf("****edge not found\n");//this should not happen
            currvertex = tgtPosEdgeV;
            break;//break for
          } else {
            if (edgeDescToLaneMapNum.find(edge_pair.first) == edgeDescToLaneMapNum.end()) {
              printf("****Unknown edge\n");//this should not happen
              currvertex = tgtvertex;//end loop
              break;//break for
            }

            uint lane = edgeDescToLaneMapNum[edge_pair.first];
            if (indexPathVec.size() < (currIndexPath - 2)) {  // almost full
              printf("p %d *****indexPathVec too small -> size: %d currIndexPath: %d\n", p, indexPathVec.size(), currIndexPath);
              indexPathVec.resize((int) (indexPathVec.size()*1.4f));// Expand 40%.
            }
            indexPathVec[currIndexPath++] = lane;
            currIndex++;  // for histogram.

            if (currIndex > 250) {
              printf("currIndex > 250\n");
              currvertex = tgtvertex;//end loop
              break;//break for
            }

            currvertex = tgtPosEdgeV;
            cont = true;
            break;//break for
          }
        }
      }

      if (cont == true) {
        continue;
      }

      // not foudn edge
      printf("****none edge works\n");//this should not happen
      //exit(0);//!!! REMOVE
      break;
    }//while find tgt

    indexPathVec[currIndexPath++] = -1; // end path with -1
    pathHistogram[currIndex]++;
    ////////////////////////////////////////////////////////////////////////////////////////////
  }

  // Resize tight indexPathVec and set everyone to the first edge (maybe do in sim?)
  printf("Final Path Size %u\n", currIndexPath);
  indexPathVec.resize(currIndexPath);
  for (int p = 0; p < trafficPersonVec.size(); p++) {
    trafficPersonVec[p].indexPathCurr = trafficPersonVec[p].indexPathInit;
  }

  { // Print histogram
    
    int maxNumSteps = INT_MAX;
    for (int h = pathHistogram.size()-1; h >=0; h--) {
      if (pathHistogram[h] != 0) {
        maxNumSteps = h;
        break;
      }
    }
    for (int h = 0; h < maxNumSteps +1; h++) {
      printf("pathHistogram,%d,%d\n", h, pathHistogram[h]);
    }
  }
  ////////////////////
  // Debug
  const bool kDebugFirstPerson = true;
  if (trafficPersonVec.size()>0 && kDebugFirstPerson) {
    int currIndex = 0;
    while (true) {
      uint laneMap = indexPathVec[trafficPersonVec[0].indexPathInit + currIndex];
      if (laneMap != -1) {
        printf("-> %u ", laneMap);
        currIndex++;
      } else {
        break;
      }
    }
    printf("\n");
  }
  printf("<< generateRoutePathsJohnson: individual routes time %d ms --> numPeople %d (No Acc %d sameSrcDst %d)\n",
    timer2.elapsed(), trafficPersonVec.size(), noAccesible, sameSrcDst);
}//
}


