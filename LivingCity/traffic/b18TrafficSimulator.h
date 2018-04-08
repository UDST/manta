/************************************************************************************************
*		@desc Class that contains the traffic simulator b2018.
*		@author igaciad
************************************************************************************************/

#ifndef LC_B18_TRAFFIC_SIMULATOR_H
#define LC_B18_TRAFFIC_SIMULATOR_H

#include "../misctools/misctools.h"

#include "b18TrafficOD.h"
#include "cudaTrafficLaneMap.h"
#include "cudaTrafficPersonShortestPath.h"

#include "../VBOPeopleJobInfoLayer.h"
#include "../VBORenderManager.h"

#include "cudaGridPollution.h"


namespace LC {

class LCUrbanMain;
//class VBORenderManager;

class B18TrafficSimulatorRender {
 public:
  std::vector<QVector3D> positions;
  std::vector<QVector3D> directions;
  std::vector<bool> goodPoints;
  int indexToRead;

  bool getInterpolated(bool goodPoint, QVector3D newPos, QVector3D newDir,
                       QVector3D &interPos, QVector3D &interDir);
};//

class B18TrafficLightRender {
 public:
  std::vector<uchar> trafficLight;
  int indexToRead;

  void getInterpolated(uchar newTrafficLight, uchar &interTrafficLight);
};//


class B18TrafficSimulator {

 public:
  B18TrafficSimulator();
  ~B18TrafficSimulator();

  // init data
  RoadGraph *simRoadGraph;
  LCUrbanMain *clientMain;

  float deltaTime;
  float cellSize;
  ushort maxWidth;
  int threadNumber;
  float avgTravelTime;


  bool initialized;
  void initSimulator(float deltaTime, float cellSize, RoadGraph *geoRoadGraph,
                     LCUrbanMain *urbanMain);

  // int numberPeople,

  //PM
  B18TrafficOD b18TrafficOD;
  CUDATrafficLaneMap cudaTrafficLaneMap;
  CUDATrafficPersonShortestPath cudaTrafficPersonShortestPath;

  void simulateInCPU_MultiPass(int numOfPasses,
                               float startTimeH, float endTimeH);
  void simulateInCPU_Onepass(float startTimeH, float endTimeH);


  void simulateInCPU(float startTimeH, float endTimeH);
  void simulateInGPU(float startTimeH, float endTimeH);

  // Lanes
  std::vector<uchar> laneMap;
  std::vector<edgeData> edgesData;
  std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> edgeDescToLaneMapNum;
  std::map<uint, RoadGraph::roadGraphEdgeDesc_BI> laneMapNumToEdgeDesc;
  void createLaneMap();

  // car path
  void generateCarPaths();

  // People
  std::vector<CUDATrafficPerson> trafficPersonVec;
  void createRandomPeople(float startTime, float endTime, int numberPeople,
                          PeopleJobInfoLayers &peopleJobInfoLayers);
  void createB2018People(float startTime, float endTime);

  void resetPeopleJobANDintersections();
  void saveODToFile() {}; // TODO
  void loadODFromFile() {};

  // Traffic lights
  std::vector<uchar> trafficLights;
  std::vector<intersectionData> intersections;

  // measurements
  std::vector<float> accSpeedPerLinePerTimeInterval;
  std::vector<float> numVehPerLinePerTimeInterval;

  void calculateAndDisplayTrafficDensity();
  void calculateAndDisplayTrafficDensity(std::vector<float>
                                         &accSpeedPerLinePerTimeInterval,
                                         std::vector<float> &numVehPerLinePerTimeInterval,
                                         std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum,
                                         std::map<uint, RoadGraph::roadGraphEdgeDesc_BI> &laneMapNumToEdgeDesc,
                                         int tNumLanes);

  void render(VBORenderManager &rendManager);
  std::vector<B18TrafficSimulatorRender> b18TrafficSimulatorRender;
  std::vector<B18TrafficLightRender> b18TrafficLightRender;
  // pollution
  CUDAGridPollution cudaGridPollution;
};
}

#endif  // LC_B18_TRAFFIC_SIMULATOR_H
