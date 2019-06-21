/************************************************************************************************
*
*		LC Project - B18 Edge data
*
*		@author igaciad
*
************************************************************************************************/

#ifndef LIVINGCITY_TRAFFIC_B18EDGEDATA_H_
#define LIVINGCITY_TRAFFIC_B18EDGEDATA_H_

#include "stdint.h"

#ifndef ushort
#define ushort uint16_t
#endif
#ifndef uint
#define uint uint32_t
#endif
#ifndef uchar
#define uchar uint8_t
#endif

const int kMaxMapWidthM = 1024;
const uint kMaskOutEdge = 0x000000;
const uint kMaskInEdge = 0x800000;
const uint kMaskLaneMap = 0x007FFFFF;


namespace LC {


enum IntersectionType {TrafficLight, Unsupervised};

// Object to store intersections (ie vertices) information
struct Intersection {
  // Start and end indexes of connections assigned to this intersection
  uint connectionGraphStart;
  uint connectionGraphEnd;

  // Start and end indexes of traffic light schedules entries assigned to this intersection
  uint trafficLightSchedulesStart;
  uint trafficLightSchedulesEnd;

  // Start and end indexes of in-lanes indexes
  uint inLanesIndexesStart;
  uint inLanesIndexesEnd;

  // Indicates index in traffic lights entries to where to keep updating on the next iteration
  uint scheduleIdx;
  // Indicates schedule group on which we are now set
  uint currentScheduleGroup;

  // Time of this intersection's last update
  float timeOfNextUpdate;

  IntersectionType intersectionType;
};

// Object to abstract whether a connection between two lanes is enabled or not
struct Connection {
  // The lane numbers have a range of possible values depending on the edge:
  // edgeNumber <= laneNumber < edgeNumber + amountOfLanesInEdge
  // Note also that
  //    laneNumber == edgeNumber is for the left-most lane
  //    laneNumber == edgeNumber + amountOfLanesInEdge - 1 is for the right-most lane
  uint inLaneNumber;
  uint outLaneNumber;

  // Flag indicating if the connection can be used
  bool enabled;

  // Start and end of indexes of connections that are blocked by this connection
  uint connectionsBlockingStart;
  uint connectionsBlockingEnd;

  uint vertexNumber;
  uint inEdgeNumber;
  uint outEdgeNumber;
};

// Object to abstract traffic lights schedules
// Each traffic light will have many groups of entries, each group sharing the schedule position
// The position indicates the order in which the positions must be enabled
struct TrafficLightScheduleEntry {
  // Vertex to which this entry belongs
  long long vertexIdx;

  // Connection which must be enabled by this schedule entry
  uint connectionIdx;

  // Group to which this schedule entry belongs
  // Entries of the same vertex with same group indicate that they must be enabled at the same
  // time
  uint scheduleGroup;

  // Amount of time assigned to this entry
  // Entries of the same vertex with same position should all have the same scheduled time
  float scheduledTime;

  // Last time this intersection was updated
  float lastUpdate;

  TrafficLightScheduleEntry(
      long long _vertexIndex,
      uint _connectionIdx,
      uint _scheduleGroup,
      float _scheduledTime) :
    vertexIdx(_vertexIndex),
    connectionIdx(_connectionIdx),
    scheduleGroup(_scheduleGroup),
    scheduledTime(_scheduledTime),
    lastUpdate(0) {}
};

struct B18EdgeData {
  uint sourceVertexIndex;
  uint targetVertexIndex;
  ushort numLines;
  uint nextInters;
  float length;
  float maxSpeedMperSec;
  bool valid = 0;
  bool startsAtHighway = false;
};

struct B18IntersectionData {
  ushort state;
  ushort stateLine;
  ushort totalInOutEdges;
  uint edge[24];// up to six arms intersection
  float nextEvent;
};


}  // namespace LC


#endif  // LIVINGCITY_TRAFFIC_B18EDGEDATA_H_

