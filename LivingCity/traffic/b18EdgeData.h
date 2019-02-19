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


// Object to store intersections (ie vertices) information
struct Intersection {
  // Start and end indexes of connections assigned to this intersection
  uint connectionGraphStart;
  uint connectionGraphEnd;

  // Start and end indexes of traffic light schedules entries assigned to this intersection
  uint trafficLightSchedulesStart;
  uint trafficLightSchedulesEnd;
};

// Object to abstract whether a connection between two lanes is enabled or not
struct Connection {
  // The lane numbers are computed as the edge number + the position of that lane in said edge
  uint inLaneNumber;
  uint outLaneNumber;

  // Flag indicating if the connection can be used
  bool enabled;

  // Some extra-info for easier debugging
  uint vertexNumber;
  uint inEdgeNumber;
  uint outEdgeNumber;
};

// Object to abstract traffic lights schedules
// Each traffic light will have many groups of entries, each group sharing the schedule position
// The position indicates the order in which the positions must be enabled
struct TrafficLightScheduleEntry {
  // Vertex to which this entry belongs
  uint vertexIdx;

  // Connection which must be enabled by this schedule entry
  uint connectionIdx;

  // Position in schedule
  // Entries of the same vertex with same position indicate that they must be enabled at the same
  // time
  uint vertexSchedulePosition;

  // Amount of time assigned to this entry
  // Entries of the same vertex with same position should al have the same scheduled time
  float scheduledTime;
};

struct B18EdgeData {
  uint originalSourceVertexIndex;
  uint originalTargetVertexIndex;
  ushort numLines;
  uint nextInters;
  float length;
  float maxSpeedMperSec;
  bool valid = 0;
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

