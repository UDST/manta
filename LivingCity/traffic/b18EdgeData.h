/************************************************************************************************
*
*		LC Project - B18 Edge data
*
*		@author igaciad
*
************************************************************************************************/

#ifndef LC_B18_EDGE_DATA_H
#define LC_B18_EDGE_DATA_H


#define ushort unsigned short
#define uint unsigned int
#define uchar unsigned char

namespace LC {

  const int kMaxMapWidthM = 1024;

struct B18EdgeData {
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
}

#endif // LC_B18_EDGE_DATA_H
