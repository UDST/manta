
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
#include <QGLWidget>

#include <QtGlobal>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include "cudaTrafficPerson.h"
#include "RoadGraph/roadGraph.h"

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "../VBOPeopleJobInfoLayer.h"
#include <random>

namespace LC {
class B18TrafficOD {

 public:
  B18TrafficOD();
  ~B18TrafficOD();

  cv::Mat *peopleDistribution;
  cv::Mat *jobDistribution;

  void randomPerson(int p, CUDATrafficPerson &person, QVector3D housePos3D,
                    QVector3D jobPos3D, float startTimeH,
                    LC::RoadGraph::roadBGLGraph_BI &roadGraph);
  void randomPerson(int p, CUDATrafficPerson &person, uint srcvertex,
                    uint tgtvertex, float startTimeH);

  void sampleDistribution(int numberToSample, PeopleJobOneLayer &distribution,
                          std::vector<QVector2D> &samples, QString &name);
  // generate random
  void createRandomPeople(int numberPerGen,
                          float startTimeH, float endTimeH,
                          std::vector<CUDATrafficPerson> &trafficPersonVec,
                          PeopleJobInfoLayers &simPeopleJobInfoLayersn,
                          LC::RoadGraph::roadBGLGraph_BI &roadGraph);

  // generate from b18
  void loadB18TrafficPeople(float startTimeH, float endTimeH,
                                std::vector<CUDATrafficPerson> &trafficPersonVec,
                                RoadGraph::roadBGLGraph_BI &roadGraph);

  void resetTrafficPersonJob(std::vector<CUDATrafficPerson> &trafficPersonVec);
  std::vector<ushort> backUpInitEdge;
};
}

#endif  // LC_B18_PM_TRAFFIC_PERSON_H
