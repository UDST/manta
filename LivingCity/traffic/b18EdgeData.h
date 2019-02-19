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


struct Intersection {
  /**
   * Object to store intersections information
   */
  uint connectionGraphStart;
  uint connectionGraphEnd;
};

struct Connection {
  /**
   * Object that represent whether a connection between lanes is enabled or not
   */
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

//struct TrafficLightScheduleEntry {

//};

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

