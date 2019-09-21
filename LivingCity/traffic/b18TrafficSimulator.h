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
#include "./simulatorConfiguration.h"


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
  B18TrafficSimulator(const SimulatorConfiguration & simulatorConfiguration);

  // init data
  std::shared_ptr<RoadGraph> simRoadGraph_shared_ptr_;
  std::shared_ptr<abm::Graph> street_graph_shared_ptr_;
  const SimulatorConfiguration & configuration_;

  LCUrbanMain *clientMain;
  int threadNumber;
  float avgTravelTime;

  B18TrafficOD b18TrafficOD_;
  SimulatorDataInitializer simulatorDataInitializer_;

  void simulateInCPU_MultiPass(void);
  void simulateInCPU(void);
  void simulateInGPU(void);

  // Internal data structures
  std::vector<uchar> laneMap;
  std::vector<B18EdgeData> edgesData;
  std::vector<LC::Connection> connections;
  std::vector<uint> connectionsBlocking;
  std::vector<LC::Intersection> updatedIntersections;
  std::vector<TrafficLightScheduleEntry> trafficLightSchedules;
  std::vector<uint> inLanesIndexes;
  std::vector<B18TrafficPerson> trafficPersonVec;
  std::vector<uint> indexPathVec;

  // Mappings between edges of the graphs and positions in laneMap
  std::map<std::shared_ptr<abm::Graph::Edge>, uint> edgeDescToLaneMapNumSP;
  std::map<uint, std::shared_ptr<abm::Graph::Edge>> laneMapNumToEdgeDescSP;
  std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> edgeDescToLaneMapNum;
  std::map<uint, RoadGraph::roadGraphEdgeDesc_BI> laneMapNumToEdgeDesc;

  // SP routing paths. Elments of this vector are the edges of street_graph_shared_ptr_, not the
  // ones from edgesData.
  std::vector<abm::graph::vertex_t> all_paths_;

  void enerateCarPaths(bool useJohnsonRouting);

#ifdef B18_RUN_WITH_GUI
  void createRandomPeople(float startTime, float endTime, int numberPeople,
                          PeopleJobInfoLayers &peopleJobInfoLayers);
#endif
  void resetPeopleJobANDintersections();
  void saveODToFile() {}; // Todo
  void loadODFromFile() {};

  // Traffic lights
  std::vector<uchar> trafficLights;
  std::vector<B18IntersectionData> intersections;

  // measurements
  std::vector<float> accSpeedPerLinePerTimeInterval;
  std::vector<float> numVehPerLinePerTimeInterval;

  void calculateAndDisplayTrafficDensity(int numOfPass);
  void savePeopleAndRoutes(int numOfPass);
  void savePeopleAndRoutesSP(int numOfPass, const std::shared_ptr<abm::Graph>& street_graph_shared_ptr_);
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
