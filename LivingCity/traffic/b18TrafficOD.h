
/************************************************************************************************
*
*		Procedural Machine Traffic Person B2018
*
*		@desc Class that generates procedually a traffic person
*		@author igaciad
*
************************************************************************************************/

#ifndef LC_B18_PM_TRAFFIC_PERSON_H
#define LC_B18_PM_TRAFFIC_PERSON_H

#include "../misctools/misctools.h"

#include <QtGlobal>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include "b18TrafficPerson.h"
#include "RoadGraph/roadGraph.h"
#include "sp/graph.h"

#include <boost/graph/dijkstra_shortest_paths.hpp>

#ifdef B18_RUN_WITH_GUI
#include "../VBOPeopleJobInfoLayer.h"
#endif

#include <random>
#include "./routing.h"
#include "./simulatorConfiguration.h"

namespace LC {
class B18TrafficOD {

 public:
  B18TrafficOD(const SimulatorConfiguration & configuration) : configuration_(configuration) {};

  cv::Mat *peopleDistribution;
  cv::Mat *jobDistribution;

  void randomPerson(int p, B18TrafficPerson &person, QVector3D housePos3D,
                    QVector3D jobPos3D, float startTimeH,
                    LC::RoadGraph::roadBGLGraph_BI &roadGraph);
  void randomPerson(int p, B18TrafficPerson &person, uint srcvertex,
                    uint tgtvertex, float startTimeH);

#ifdef B18_RUN_WITH_GUI
  void sampleDistribution(int numberToSample, PeopleJobOneLayer &distribution,
                          std::vector<QVector2D> &samples, QString &name);
  // generate random
  void createRandomPeople(int numberPerGen,
                          float startTimeH, float endTimeH,
                          std::vector<B18TrafficPerson> &trafficPersonVec,
                          PeopleJobInfoLayers &simPeopleJobInfoLayersn,
                          LC::RoadGraph::roadBGLGraph_BI &roadGraph);
#endif

  void loadB18TrafficPeople(std::vector<B18TrafficPerson> &trafficPersonVec);

  void resetTrafficPersonJob(std::vector<B18TrafficPerson> &trafficPersonVec);

 private:
  const SimulatorConfiguration & configuration_;
};
}

#endif  // LC_B18_PM_TRAFFIC_PERSON_H
