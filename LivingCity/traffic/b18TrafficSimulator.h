/************************************************************************************************
*		@desc Class that contains the traffic simulator b2018.
*		@author igaciad
************************************************************************************************/

#ifndef LC_B18_TRAFFIC_SIMULATOR_H
#define LC_B18_TRAFFIC_SIMULATOR_H

#include "../misctools/misctools.h"

#include "b18TrafficOD.h"
#include "b18TrafficLaneMap.h"

#ifdef B18_RUN_WITH_GUI
#include "../VBOPeopleJobInfoLayer.h"
#include "../VBORenderManager.h"
#endif

#include "b18GridPollution.h"


namespace LC {

class LCUrbanMain;

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
  B18TrafficSimulator(float deltaTime, RoadGraph *geoRoadGraph, LCUrbanMain *urbanMain = nullptr);
  ~B18TrafficSimulator();

  // init data
  RoadGraph *simRoadGraph;
  LCUrbanMain *clientMain;

  float deltaTime;
  int threadNumber;
  float avgTravelTime;

  //PM
  B18TrafficOD b18TrafficOD;
  B18TrafficLaneMap b18TrafficLaneMap;

  void simulateInCPU_MultiPass(int numOfPasses,
                               float startTimeH, float endTimeH, bool useJohnsonRouting);
  void simulateInCPU_Onepass(float startTimeH, float endTimeH,
                             bool useJohnsonRouting);
  void simulateInCPU(float startTimeH, float endTimeH);

  //void simulateInGPU(int numOfPasses, float startTimeH, float endTimeH,
  //                   bool useJohnsonRouting, bool useSP);
  
  void simulateInGPU(int numOfPasses, float startTimeH, float endTimeH,
    bool useJohnsonRouting, bool useSP, const std::shared_ptr<abm::Graph>& graph_, std::vector<abm::graph::vertex_t> paths_SP, bool saveFiles, float s_0, std::vector<std::array<abm::graph::vertex_t, 2>> all_od_pairs, std::vector<float> dep_times, bool usePrevPaths, std::string networkPathSP);

  // Lanes
  std::vector<uchar> laneMap;
  std::vector<B18EdgeData> edgesData;
  std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> edgeDescToLaneMapNum;
  std::map<uint, RoadGraph::roadGraphEdgeDesc_BI> laneMapNumToEdgeDesc;
  std::map<uint, std::shared_ptr<abm::Graph::Edge>> laneMapNumToEdgeDescSP;
  //std::map<std::shared_ptr<abm::Graph::Edge>, uint> edgeDescToLaneMapNumSP;
  tsl::robin_map<std::shared_ptr<abm::Graph::Edge>, uint> edgeDescToLaneMapNumSP;
  void createLaneMap();
  void createLaneMapSP(const std::shared_ptr<abm::Graph>& graph_);

  // car path
  void generateCarPaths(bool useJohnsonRouting);

  // People
  std::vector<B18TrafficPerson> trafficPersonVec;
  std::vector<uint> indexPathVec;
  std::vector<std::vector<float>> timeMatrix;
  std::vector<std::vector<float>> routeShareMatrix;
  std::vector<std::vector<abm::graph::vertex_t>> allPathsMatrix;


  std::vector<std::array<abm::graph::vertex_t, 2>> filtered_od_pairs_;
  std::vector<float> filtered_dep_times_;


#ifdef B18_RUN_WITH_GUI
  void createRandomPeople(float startTime, float endTime, int numberPeople,
                          PeopleJobInfoLayers &peopleJobInfoLayers);
#endif
  void createB2018People(float startTime, float endTime, int limitNumPeople, bool addRandomPeople, bool useSP);
  
  void createB2018PeopleSP(float startTime, float endTime, int limitNumPeople, bool addRandomPeople, const std::shared_ptr<abm::Graph>& graph_, float a, float b, float T, std::vector<std::array<abm::graph::vertex_t, 2>> all_od_pairs_);

  void resetPeopleJobANDintersections();
  void saveODToFile() {}; // TODO
  void loadODFromFile() {};

  // Traffic lights
  std::vector<uchar> trafficLights;
  std::vector<B18IntersectionData> intersections;

  // measurements
  std::vector<float> accSpeedPerLinePerTimeInterval;
  std::vector<float> numVehPerLinePerTimeInterval;

  void calculateAndDisplayTrafficDensity(int numOfPass);
  void savePeopleAndRoutes(int numOfPass);
  void savePeopleAndRoutesSP(int numOfPass, const std::shared_ptr<abm::Graph>& graph_, std::vector<abm::graph::vertex_t> paths_SP, int start_time, int end_time);
#ifdef B18_RUN_WITH_GUI
  void render(VBORenderManager &rendManager);
#endif
  std::vector<B18TrafficSimulatorRender> b18TrafficSimulatorRender;
  std::vector<B18TrafficLightRender> b18TrafficLightRender;

  // pollution
  B18GridPollution gridPollution;
};
}

#endif  // LC_B18_TRAFFIC_SIMULATOR_H
