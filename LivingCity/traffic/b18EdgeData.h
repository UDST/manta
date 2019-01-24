/************************************************************************************************
*
*		LC Project - B18 Edge data
*
*		@author igaciad
*
************************************************************************************************/

#ifndef LIVINGCITY_TRAFFIC_B18EDGEDATA_H_
#define LIVINGCITY_TRAFFIC_B18EDGEDATA_H_

#include <map>

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


struct ConnectionsInfo {
  int in_lane_number;
  int out_lane_number;
  bool enabled;
};

struct B18EdgeData {
  int original_src_vertex_index;
  int original_target_vertex_index;
  ushort numLines;
  uint nextInters;
  float length;
  float maxSpeedMperSec;
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

