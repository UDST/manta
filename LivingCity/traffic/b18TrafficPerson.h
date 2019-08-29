/************************************************************************************************
*
*		LC Project - B18 Traffic Person
*
*
*		@desc Class that contains the info of a person
*		@author igaciad
*
************************************************************************************************/

#ifndef LC_B18_TRAFFIC_PERSON_H
#define LC_B18_TRAFFIC_PERSON_H

#include <climits>

namespace LC {

struct B18TrafficPerson {
  unsigned int init_intersection;
  unsigned int end_intersection;
  float time_departure;

  unsigned short active;//0 inactive 1 active 2 finished
  unsigned short numOfLaneInEdge;//number of lane in that edge

  float posInLaneM;

  //////////////////////////
  // current edge (from edgeData)
  unsigned short edgeNumLanes;  //total number of lanes in that edge
  unsigned int edgeNextInters;
  float length;
  float maxSpeedMperSec;
  /////////////////////////
  // to check next edge
  unsigned short nextEdgeNumLanes;
  unsigned short nextEdgeNextInters;
  float nextEdgeLength;
  float nextEdgemaxSpeedMperSec;
  ///////////////////////////
  unsigned int indexPathInit;
  unsigned int indexPathCurr;

  /*
   * The main idea behind stop junction is to handle them as obstacles in the IDM.
   * Cars however should only need to stop once before the stop signal, which is what the following
   * members are for.
   * */
  bool isApproachingStopJunction = false;
  int distanceUntilIntersection = INT_MAX >> 2;
  bool stoppedBeforeNextIntersection = false;

  // data
  unsigned short num_steps;
  unsigned int color;
  float co;
  float gas;
  // IDM
  float v;//current velocity
  float a;//acceleration
  float b;//break
  float T;// Time heading

  // lane changing
  unsigned short LC_initOKLanes;
  unsigned short LC_endOKLanes;
  unsigned short LC_stateofLaneChanging;
};

}

#endif  // LC_B18_TRAFFIC_PERSON_H
