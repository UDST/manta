/************************************************************************************************
*		@desc Class to load the B2018 road graph
*		@author igarciad
************************************************************************************************/

#include "roadGraphB2018Loader.h"

#include "Geometry/client_geometry.h"
#include "LC_UrbanMain.h"
#include "LC_GLWidget3D.h"

#include "global.h"
#include "bTraffic/bTrafficIntersection.h"

#include <QHash>
#include <QVector2D>

namespace LC {

std::vector<DemandB2018> RoadGraphB2018::demandB2018;
int RoadGraphB2018::totalNumPeople;

void updateMinMax2(QVector2D &newPoint, QVector2D &minBox, QVector2D &maxBox) {
  if (newPoint.x() < minBox.x()) {
    minBox.setX(newPoint.x());
  }

  if (newPoint.y() < minBox.y()) {
    minBox.setY(newPoint.y());
  }

  if (newPoint.x() > maxBox.x()) {
    maxBox.setX(newPoint.x());
  }

  if (newPoint.y() > maxBox.y()) {
    maxBox.setY(newPoint.y());
  }
}//

///////////////////////////////////////////////////////////

static QVector2D projLatLonToWorldMercator(float lat, float lon,
    bool isDeg = false) {
  // Ellipsoid model constants (actual values here are for WGS84)
  const float sm_a = 6378137.0f;
  const float sm_b = 6356752.314f;

  QVector2D result;
  float lon0 = 0;

  if (isDeg == true) {
    lat = lat / 180.0f * M_PI;
    lon = lon / 180.0f * M_PI;
  }

  //printf("lat %.2f lon %.2f\n", lat, lon);
  //result.setX(sm_a*(lon - lon0));
  //result.setY(sm_a*log((sin(lat) + 1) / cos(lat)));

  result.setX(sm_a * (cos(lat) * cos(lon)));
  result.setY(sm_a * (cos(lat) * sin(lon)));

  //qDebug() << result;
  return  result;
}//

void saveSetToFile(QSet<uint64> &set, QString &filename) {
  QFile file(filename);

  if (file.open(QIODevice::ReadWrite)) {
    QTextStream stream(&file);
    QSetIterator<uint64> nA(set);

    while (nA.hasNext()) {
      stream << nA.next() << "\n";
    }
  }
}

//////////////////////////////////////////////////////////

void RoadGraphB2018::loadB2018RoadGraph(RoadGraph &inRoadGraph,
                                        LCGLWidget3D *glWidget3D) {

  printf(">>loadB2018RoadGraph\n");
  printf(">>Remove\n");
  inRoadGraph.myRoadGraph.clear();
  inRoadGraph.myRoadGraph_BI.clear();
  glWidget3D->cg.geoZone.blocks.clear();
  glWidget3D->vboRenderManager.removeAllStreetElementName("tree");
  glWidget3D->vboRenderManager.removeAllStreetElementName("streetLamp");
  printf("<<Remove\n");


  printf("loadB2018RoadGraph\n");

  /////////////////////////////////////////////////
  // READ NODES

  QString fileName = "berkeley_2018/bay_area_full_strongly_nodes.csv";
  //QString fileName="data/Dynameq/smallTestNet_base.dqt";
  QFile baseFile(fileName); // Create a file handle for the file named

  QString line;

  if (!baseFile.open(QIODevice::ReadOnly | QIODevice::Text)) { // Open the file
    printf("Can't open file '%s'\n", fileName.toUtf8().constData());
    return;
  }

  QTextStream stream(&baseFile); // Set the stream to read from myFile
  QTime timer;
  timer.start();
  QVector2D minBox(FLT_MAX, FLT_MAX);
  QVector2D maxBox(-FLT_MAX, -FLT_MAX);
  QHash<uint64, QVector2D> osmidToVertexLoc;
  QHash<uint64, uchar> osmidToBType; // node type

  QHash<QString, uchar> bTypeStringTobType;
  bTypeStringTobType[""] = 0;
  bTypeStringTobType["motorway_junction"] = 1;
  bTypeStringTobType["traffic_signals"] = 2;
  bTypeStringTobType["stop"] = 3;
  bTypeStringTobType["turning_circle"] = 4;

  QStringList headers = stream.readLine().split(",");
  const int indexOsmid = headers.indexOf("osmid");
  const int indexX = headers.indexOf("x");
  const int indexY = headers.indexOf("y");
  const int indexHigh = headers.indexOf("highway");
  printf("Node Index %d %d %d %d\n", indexOsmid, indexX, indexY, indexHigh);

  while (!stream.atEnd()) {
    line = stream.readLine();
    QStringList fields = line.split(',', QString::SkipEmptyParts);

    if (indexX >= fields.size() || indexY >= fields.size()) {
      qDebug() << "ERROR line " << line << " --> SKIP";
      continue;
    }

    float x = fields[indexX].toFloat();
    float y = fields[indexY].toFloat();
    //qDebug() << "x " << x << " y " << y;
    uint64 osmid = fields[indexOsmid].toLongLong();
    osmidToVertexLoc[osmid] = QVector2D(x, y);
    updateMinMax2(QVector2D(x, y), minBox, maxBox);

    if (indexHigh >= fields.size()) {
      osmidToBType[osmid] = 0;
    } else {
      QString bType = fields[indexHigh];
      osmidToBType[osmid] = (!bTypeStringTobType.contains(bType)) ? 0 :
                            bTypeStringTobType[bType];
    }

    //printf("MinBox %f %f MaxBox %f %f--> %f %f\n", minBox.x(), minBox.y(), maxBox.x(), maxBox.y(), maxBox.x() - minBox.x(), maxBox.y() - minBox.y());
  }

  printf(">>Degrees --> MinBox %f %f MaxBox %f %f--> %f %f\n", minBox.x(),
         minBox.y(), maxBox.x(), maxBox.y(), maxBox.x() - minBox.x(),
         maxBox.y() - minBox.y());

  // Update coordenades to East-North-Up coordinates;

  const float lat0 = (maxBox.x() + minBox.x()) * 0.5f;
  const float lon0 = (maxBox.y() + minBox.y()) * 0.5f;
  minBox = QVector2D(FLT_MAX, FLT_MAX);
  maxBox = QVector2D(-FLT_MAX, -FLT_MAX);
  QHash<uint64, QVector2D>::iterator i;

  for (i = osmidToVertexLoc.begin(); i != osmidToVertexLoc.end(); ++i) {
    //printf("1 ANG: %.2f %.2f MinBox %.2f %.2f MaxBox %.2f %.2f--> %.2f %.2f\n", i.value().x(), i.value().y(), minBox.x(), minBox.y(), maxBox.x(), maxBox.y(), maxBox.x() - minBox.x(), maxBox.y() - minBox.y());

    osmidToVertexLoc[i.key()] = projLatLonToWorldMercator(i.value().x(),
                                i.value().y(), /*isDeg=*/true);
    updateMinMax2(osmidToVertexLoc[i.key()], minBox, maxBox);
    //printf("1 M: %.2f %.2f MinBox %.2f %.2f MaxBox %.2f %.2f--> %.2f %.2f\n", osmidToVertexLoc[i.key()].x(), osmidToVertexLoc[i.key()].y(), minBox.x(), minBox.y(), maxBox.x(), maxBox.y(), maxBox.x() - minBox.x(), maxBox.y() - minBox.y());
    //exit(-1);
  }

  printf(">>Meters --> MinBox %f %f MaxBox %f %f--> %f %f\n", minBox.x(),
         minBox.y(), maxBox.x(), maxBox.y(), maxBox.x() - minBox.x(),
         maxBox.y() - minBox.y());

  // TERRAIN
  printf(">> Process terrain");
  float scale = 1.0f;
  float sqSideSz = std::max(maxBox.x() - minBox.x(),
                            maxBox.y() - minBox.y()) * scale * 0.5f; // half side
  QVector3D centerV(-minBox.x(), -minBox.y(), 0);
  QVector3D centerAfterSc(-sqSideSz, -sqSideSz, 0);
  G::boundingPolygon.clear();
  G::boundingPolygon.push_back(QVector3D(sqSideSz, -sqSideSz, 0.0f));
  glWidget3D->vboRenderManager.changeTerrainDimensions(sqSideSz * 2 + 400.0f,
      200);

  QString sfo_path("data/b2018.png");
  glWidget3D->vboRenderManager.vboTerrain.loadTerrain(sfo_path);

  ///////////////////////////////
  // ADD NODES
  printf(">> Process nodes");
  std::vector<RoadGraph::roadGraphVertexDesc> vertex;
  std::vector<RoadGraph::roadGraphVertexDesc> vertex_SIM;
  vertex.resize(osmidToVertexLoc.size());
  vertex_SIM.resize(osmidToVertexLoc.size());

  int index = 0;
  QHash<uint64, int> dynIndToInd;
  QHash<int, uint64> indToOsmid;

  for (i = osmidToVertexLoc.begin(); i != osmidToVertexLoc.end(); ++i) {
    uint64 ind = i.key();

    float x = osmidToVertexLoc[ind].x();
    float y = osmidToVertexLoc[ind].y();
    uchar bType = osmidToBType[ind];

    QVector3D pos(x, y, 0);
    pos += centerV;//center
    pos *= scale;
    pos += centerAfterSc;
    pos.setX(pos.x() * -1.0f); // seems vertically rotated

    vertex[index] = boost::add_vertex(inRoadGraph.myRoadGraph_BI);
    inRoadGraph.myRoadGraph_BI[vertex[index]].pt = pos;
    inRoadGraph.myRoadGraph_BI[vertex[index]].bType = bType;

    vertex_SIM[index] = boost::add_vertex(inRoadGraph.myRoadGraph);
    inRoadGraph.myRoadGraph[vertex_SIM[index]].pt = pos;
    inRoadGraph.myRoadGraph[vertex_SIM[index]].bType = bType;

    dynIndToInd[ind] = index;
    indToOsmid[index] = ind;
    index++;
  }

  printf("\n*** Readed in %d --> #Nod %d\n", timer.elapsed(), index);

  ///////////////////////////////
  // EDGES
  printf(">> Process edges\n");
  fileName = "berkeley_2018/edges_speed_capacity.csv";
  QFile linkFile(fileName); // Create a file handle for the file named

  if (!linkFile.open(QIODevice::ReadOnly | QIODevice::Text)) { // Open the file
    printf("Can't open file '%s'\n", fileName.toUtf8().constData());
    return;
  }

  QTextStream streamL(&linkFile); // Set the stream to read

  headers = (streamL.readLine()).split(",");
  const int indexId = headers.indexOf("uniqueid");
  const int indexU = headers.indexOf("u");
  const int indexV = headers.indexOf("v");
  const int indexLen = headers.indexOf("length");
  const int indexLanes = headers.indexOf("lanes");
  const int indexSpeedMH = headers.indexOf("speed_mph");

  printf("Link Index %d %d %d %d %d %d\n", indexId, indexU, indexV, indexLen,
         indexLanes, indexSpeedMH);
  qDebug() << sizeof(long);
  qDebug() << sizeof(uint64);
  QHash<int, std::pair<uint, uint>> dynEdgToEdge;
  std::pair<RoadGraph::roadGraphEdgeDesc_BI, bool> e0_pair;
  std::pair<RoadGraph::roadGraphEdgeDesc, bool> e0_pair_SIMP;
  float totalLeng = 0;
  int numEdges = 0;
  QSet<uint64> noAvailableNodes;
  const bool saveNoAvailableNodes = false;

  while (!streamL.atEnd()) {
    line = streamL.readLine();
    QStringList fields = line.split(',', QString::SkipEmptyParts);

    if (fields.size() < 3) {
      qDebug() << "ERROR line " << line << " --> SKIP";
      continue;
    }

    int ind = fields[indexId].toInt();

    uint64 start = fields[indexU].toLongLong();
    uint64 end = fields[indexV].toLongLong();

    if ((!dynIndToInd.contains(start)) || (!dynIndToInd.contains(end))) {
      if (saveNoAvailableNodes) {
        if (!dynIndToInd.contains(start)) {
          noAvailableNodes.insert(start);
        }

        if (!dynIndToInd.contains(end)) {
          noAvailableNodes.insert(end);
        }
      }

      qDebug() << "NO CONTAINS: start" << start << " end " << end;
      exit(-1);
      continue;
    }

    //qDebug() << "line" << line;
    //qDebug() << fields;
    //qDebug() << "start" << start << " end " << end;
    float length = fields[indexLen].toFloat();
    int numLanes = fields[indexLanes].toInt();
    float speedMS = fields[indexSpeedMH].toFloat() * 0.44704f; //m/h --> m/sec

    //printf("%d %d of %d): Leng %.2f #lanes %d speed %.2f\n", dynIndToInd[start], dynIndToInd[end], index, length, numLanes, speedMS);

    totalLeng += length;

    // add edge if not already there or update num lanes
    if (boost::edge(vertex_SIM[dynIndToInd[start]], vertex_SIM[dynIndToInd[end]],
                    inRoadGraph.myRoadGraph).second == false) {
      e0_pair_SIMP = boost::add_edge(vertex[dynIndToInd[start]],
                                     vertex[dynIndToInd[end]], inRoadGraph.myRoadGraph);
      inRoadGraph.myRoadGraph[e0_pair_SIMP.first].numberOfLanes = numLanes;
      inRoadGraph.myRoadGraph[e0_pair_SIMP.first].edgeLength = length;

    } else {
      inRoadGraph.myRoadGraph[boost::edge(vertex_SIM[dynIndToInd[start]],
                                          vertex_SIM[dynIndToInd[end]],
                                          inRoadGraph.myRoadGraph).first].numberOfLanes += numLanes;
    }

    e0_pair = boost::add_edge(vertex[dynIndToInd[start]], vertex[dynIndToInd[end]],
                              inRoadGraph.myRoadGraph_BI);
    inRoadGraph.myRoadGraph_BI[e0_pair.first].numberOfLanes = numLanes;
    inRoadGraph.myRoadGraph_BI[e0_pair.first].edgeLength = length;
    inRoadGraph.myRoadGraph_BI[e0_pair.first].maxSpeedMperSec = speedMS;

    // add to edge
    dynEdgToEdge[ind] = std::make_pair(dynIndToInd[start], dynIndToInd[end]);

    /*if (++numEdges > 1000) {
      break;
    }*/
  }

  // Save no available nodes to file.
  if (saveNoAvailableNodes) {
    QString filename = "noAvailableNodes.txt";
    saveSetToFile(noAvailableNodes, filename);
  }

  ///////////////////////////////
  // DEMAND
  printf(">> Process demand\n");
  fileName = "berkeley_2018/od.csv";
  QFile demandFile(fileName); // Create a file handle for the file named

  if (!demandFile.open(QIODevice::ReadOnly | QIODevice::Text)) { // Open the file
    printf("Can't open file '%s'\n", fileName.toUtf8().constData());
    return;
  }

  QTextStream streamD(&demandFile); // Set the stream to read
  headers = (streamD.readLine()).split(",");
  const int numPeopleIndex = headers.indexOf("PERNO");
  const int origIndex = headers.indexOf("orig");
  const int destIndex = headers.indexOf("dest");

  printf("Demand Index %d %d %d\n", numPeopleIndex, origIndex, destIndex);
  QSet<uint64> noAvailableNodesDemand;
  const bool saveNoAvailableNodesDemand = true;
  totalNumPeople = 0;

  while (!streamD.atEnd()) {
    line = streamD.readLine();
    QStringList fields = line.split(',', QString::SkipEmptyParts);

    if (fields.size() < 4) {
      qDebug() << "ERROR line " << line << " --> SKIP";
      continue;
    }

    int numPeople = fields[numPeopleIndex].toInt();
    totalNumPeople += numPeople;

    uint64 start = fields[origIndex].toLongLong();
    uint64 end = fields[destIndex].toLongLong();

    if ((!dynIndToInd.contains(start)) || (!dynIndToInd.contains(end))) {
      if (saveNoAvailableNodesDemand) {
        if (!dynIndToInd.contains(start)) {
          noAvailableNodesDemand.insert(start);
        }
        if (!dynIndToInd.contains(end)) {
          noAvailableNodesDemand.insert(end);
        }
      }
      //qDebug() << "NO CONTAINS: start" << start << " end " << end;
      exit(-1);
      continue;
    }
    demandB2018.push_back(DemandB2018(numPeople, dynIndToInd[start], dynIndToInd[end]));
  }

  // Save no available nodes to file.
  if (saveNoAvailableNodesDemand) {
    QString filename = "noAvailableNodesDemand.txt";
    saveSetToFile(noAvailableNodesDemand, filename);
  }

  ///////////////////////////////////////

  printf("\n*** Readed in %d\n", timer.elapsed());
  printf("       --> #Nod %d\n", index);
  printf("       --> #Edg %d -> Total length %2.2f\n", dynEdgToEdge.size(),
         totalLeng);
  printf("       --> #Peo %d\n", totalNumPeople);
  //printf("\nNodes readed in %d Nod %d Cen %d Link %d\n", timer.elapsed(), osmidToVertexLoc.size(), centroids.size(), links.size());
}

}  // namespace LC
