/************************************************************************************************
*		@desc Class to load the B2018 road graph
*		@author igarciad
************************************************************************************************/

#include <QHash>
#include <QVector2D>
#include <stdexcept>

#include "Geometry/client_geometry.h"
#include "bTraffic/bTrafficIntersection.h"
#include "global.h"
#include "roadGraphB2018Loader.h"
#include "OSMConstants.h"

namespace LC {

using namespace std::chrono;

std::vector<DemandB2018> RoadGraphB2018::demandB2018;
int RoadGraphB2018::totalNumPeople;
QHash<int, uint64_t> RoadGraphB2018::indToOsmid;

void updateMinMax2(const QVector2D &newPoint, QVector2D &minBox, QVector2D &maxBox) {
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

void saveSetToFile(QSet<uint64_t> &set, QString &filename) {
  QFile file(filename);

  if (file.open(QIODevice::ReadWrite)) {
    QTextStream stream(&file);
    QSetIterator<uint64_t> nA(set);

    while (nA.hasNext()) {
      stream << nA.next() << "\n";
    }
  }
}

//////////////////////////////////////////////////////////

void RoadGraphB2018::loadB2018RoadGraph(
    std::shared_ptr<RoadGraph> inRoadGraph,
    const QString & networkPath) {
  inRoadGraph->myRoadGraph.clear();
  inRoadGraph->myRoadGraph_BI.clear();

  QString nodesFileName = networkPath + "nodes.csv";
  std::cerr
    << "[Log] Loading nodes (using \"" << nodesFileName.toUtf8().constData()
    << "\"" << " as input)." << std::endl;

  QFile baseFile(nodesFileName); // Create a file handle for the file named
  QString line;
  if (!baseFile.open(QIODevice::ReadOnly | QIODevice::Text)) { // Open the file
    throw std::invalid_argument("RoadGraphB2018::loadB2018RoadGraph -> Can't open nodes files.");
  }

  QTextStream stream(&baseFile); // Set the stream to read from myFile
  QTime timer;
  timer.start();
  QVector2D minBox(FLT_MAX, FLT_MAX);
  QVector2D maxBox(-FLT_MAX, -FLT_MAX);
  QHash<uint64_t, QVector2D> osmidToVertexLoc;
  QHash<uint64_t, QVector2D> osmidToOriginalLoc;
  QHash<uint64_t, OSMConstant> OSMIdToOSMType; // node type

  QStringList headers = stream.readLine().split(",");
  const int indexOsmid = headers.indexOf("osmid");
  const int indexX = headers.indexOf("x");
  const int indexY = headers.indexOf("y");
  const int indexHigh = headers.indexOf("highway");

  while (!stream.atEnd()) {
    line = stream.readLine();
    QStringList fields = line.split(',');

    if (indexX >= fields.size() || indexY >= fields.size()) {
      qDebug() << "ERROR line " << line << " --> SKIP";
      continue;
    }

    float x = fields[indexX].toFloat();
    float y = fields[indexY].toFloat();
    uint64_t osmid = fields[indexOsmid].toLongLong();
    osmidToVertexLoc[osmid] = QVector2D(x, y);
    osmidToOriginalLoc[osmid] = QVector2D(x, y);
    updateMinMax2(QVector2D(x, y), minBox, maxBox);

    OSMIdToOSMType[osmid] = mapStringToOSMConstant(fields[indexHigh].toStdString());
  }

  // Update coordenades to East-North-Up coordinates;
  const float lat0 = (maxBox.x() + minBox.x()) * 0.5f;
  const float lon0 = (maxBox.y() + minBox.y()) * 0.5f;
  minBox = QVector2D(FLT_MAX, FLT_MAX);
  maxBox = QVector2D(-FLT_MAX, -FLT_MAX);
  QHash<uint64_t, QVector2D>::iterator i;

  for (i = osmidToVertexLoc.begin(); i != osmidToVertexLoc.end(); ++i) {
    osmidToVertexLoc[i.key()] = projLatLonToWorldMercator(i.value().x(),
                                i.value().y(), /*isDeg=*/true);
    updateMinMax2(osmidToVertexLoc[i.key()], minBox, maxBox);
  }

  // TERRAIN
  float scale = 1.0f;
  float sqSideSz = std::max<float>(maxBox.x() - minBox.x(),
                            maxBox.y() - minBox.y()) * scale * 0.5f; // half side
  QVector3D centerV(-minBox.x(), -minBox.y(), 0);
  QVector3D centerAfterSc(-sqSideSz, -sqSideSz, 0);
  G::boundingPolygon.clear();
  G::boundingPolygon.push_back(QVector3D(sqSideSz, -sqSideSz, 0.0f));

  ///////////////////////////////
  // ADD NODES
  std::vector<RoadGraph::roadGraphVertexDesc> vertex;
  std::vector<RoadGraph::roadGraphVertexDesc> vertex_SIM;

  int index = 0;
  QHash<uint64_t, int> dynIndToInd;

  // Create vertices
  for (i = osmidToVertexLoc.begin(); i != osmidToVertexLoc.end(); ++i) {
    const uint64_t ind = i.key();

    const auto new_bi_vertex_descriptor = boost::add_vertex(inRoadGraph->myRoadGraph_BI);
    vertex.push_back(new_bi_vertex_descriptor);
    const auto new_vertex_descriptor = boost::add_vertex(inRoadGraph->myRoadGraph);
    vertex_SIM.push_back(new_vertex_descriptor);

    dynIndToInd[ind] = index;
    indToOsmid[index] = ind;
    index++;
  }

  // Load data into vertices
  for (i = osmidToVertexLoc.begin(); i != osmidToVertexLoc.end(); ++i) {
    const uint64_t ind = i.key();
    const auto index = dynIndToInd[ind];

    const float x = osmidToVertexLoc[ind].x();
    const float y = osmidToVertexLoc[ind].y();

    QVector3D pos(x, y, 0);
    pos += centerV;//center
    pos *= scale;
    pos += centerAfterSc;
    pos.setX(pos.x() * -1.0f); // seems vertically rotated

    const auto bi_vertex_descriptor = vertex[index];
    inRoadGraph->myRoadGraph_BI[bi_vertex_descriptor].x = osmidToOriginalLoc[ind].x();
    inRoadGraph->myRoadGraph_BI[bi_vertex_descriptor].x = osmidToOriginalLoc[ind].x();
    inRoadGraph->myRoadGraph_BI[bi_vertex_descriptor].y = osmidToOriginalLoc[ind].y();
    inRoadGraph->myRoadGraph_BI[bi_vertex_descriptor].pt = pos;
    inRoadGraph->myRoadGraph[bi_vertex_descriptor].intersectionType = OSMIdToOSMType[ind];

    const auto vertex_descriptor = vertex_SIM[index];
    inRoadGraph->myRoadGraph[vertex_descriptor].pt = pos;
    inRoadGraph->myRoadGraph[vertex_descriptor].intersectionType = OSMIdToOSMType[ind];
  }

  QString edgeFileName = networkPath + "edges.csv";
  std::cerr
    << "[Log] Loading edges (using \"" << edgeFileName.toUtf8().constData()
    << "\"" << " as input)." << std::endl;

  QFile linkFile(edgeFileName); // Create a file handle for the file named
  if (!linkFile.open(QIODevice::ReadOnly | QIODevice::Text)) { // Open the file
    throw std::invalid_argument("RoadGraphB2018::loadB2018RoadGraph -> Can't open edges files.");
  }

  QTextStream streamL(&linkFile); // Set the stream to read

  headers = (streamL.readLine()).split(",");
  const int indexId = headers.indexOf("uniqueid");
  const int indexU = headers.indexOf("u");
  const int indexV = headers.indexOf("v");
  const int indexLen = headers.indexOf("length");
  const int indexLanes = headers.indexOf("lanes");
  const int indexSpeedMH = headers.indexOf("speed_mph");

  QHash<int, std::pair<uint, uint>> dynEdgToEdge;
  std::pair<RoadGraph::roadGraphEdgeDesc_BI, bool> e0_pair;
  std::pair<RoadGraph::roadGraphEdgeDesc, bool> e0_pair_SIMP;
  float totalLeng = 0;
  int numEdges = 0;
  QSet<uint64_t> noAvailableNodes;
  const bool saveNoAvailableNodes = false;

  while (!streamL.atEnd()) {
    line = streamL.readLine();
    QStringList fields = line.split(',', QString::SkipEmptyParts);

    if (fields.size() < 3) {
      qDebug() << "ERROR line " << line << " --> SKIP";
      continue;
    }

    uint ind = fields[indexId].toInt();

    uint64_t start = fields[indexU].toLongLong();
    uint64_t end = fields[indexV].toLongLong();

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

    float length = fields[indexLen].toFloat();
    int numLanes = std::max<int>(fields[indexLanes].toInt(),1); // at least one
    float speedMS = std::max<float>(0.01f, fields[indexSpeedMH].toFloat() * 0.44704f); //m/h --> m/sec // force to have a speed

    //printf("%d %d of %d): Leng %.2f #lanes %d speed %.2f\n", dynIndToInd[start], dynIndToInd[end], index, length, numLanes, speedMS);

    totalLeng += length;

    // add edge if not already there or update num lanes
    if (boost::edge(vertex_SIM[dynIndToInd[start]], vertex_SIM[dynIndToInd[end]],
                    inRoadGraph->myRoadGraph).second == false) {
      e0_pair_SIMP = boost::add_edge(vertex[dynIndToInd[start]],
                                     vertex[dynIndToInd[end]], inRoadGraph->myRoadGraph);
      inRoadGraph->myRoadGraph[e0_pair_SIMP.first].numberOfLanes = numLanes;
      inRoadGraph->myRoadGraph[e0_pair_SIMP.first].edgeLength = length;

    } else {
      inRoadGraph->myRoadGraph[boost::edge(vertex_SIM[dynIndToInd[start]],
                                          vertex_SIM[dynIndToInd[end]],
                                          inRoadGraph->myRoadGraph).first].numberOfLanes += numLanes;
    }

    e0_pair = boost::add_edge(vertex[dynIndToInd[start]], vertex[dynIndToInd[end]],
                              inRoadGraph->myRoadGraph_BI);
    inRoadGraph->myRoadGraph_BI[e0_pair.first].numberOfLanes = numLanes;
    inRoadGraph->myRoadGraph_BI[e0_pair.first].edgeLength = length;
    inRoadGraph->myRoadGraph_BI[e0_pair.first].maxSpeedMperSec = speedMS;
    inRoadGraph->myRoadGraph_BI[e0_pair.first].faci = ind;
    // add to edge
    dynEdgToEdge[ind] = std::make_pair(dynIndToInd[start], dynIndToInd[end]);
  }

  // Save no available nodes to file.
  if (saveNoAvailableNodes) {
    QString filename = "noAvailableNodes.txt";
    saveSetToFile(noAvailableNodes, filename);
  }

  QString odFileName = networkPath + "od_demand.csv";
  std::cerr
    << "[Log] Loading demand (using \"" << odFileName.toUtf8().constData()
    << "\"" << " as input)." << std::endl;

  QFile demandFile(odFileName); // Create a file handle for the file named

  if (!demandFile.open(QIODevice::ReadOnly | QIODevice::Text)) { // Open the file
    throw std::invalid_argument("RoadGraphB2018::loadB2018RoadGraph -> Can't open od demand files.");
  }

  QTextStream streamD(&demandFile); // Set the stream to read
  headers = (streamD.readLine()).split(",");
  const int numPeopleIndex = headers.indexOf("PERNO");
  const int origIndex = headers.indexOf("origin");
  const int destIndex = headers.indexOf("destination");

  if (numPeopleIndex < 0 || origIndex < 0 || destIndex < 0)
    throw std::runtime_error(
      "RoadGraphB2018::loadB2018RoadGraph -> Demand file has incorrect format. Required columns are 'PERNO', 'origin' and 'destination'.");

  QSet<uint64_t> noAvailableNodesDemand;
  const bool saveNoAvailableNodesDemand = false;
  totalNumPeople = 0;

  while (!streamD.atEnd()) {
    line = streamD.readLine();
    QStringList fields = line.split(',', QString::SkipEmptyParts);

    if (fields.size() < 4) {
      qDebug() << "ERROR line " << line << " --> SKIP";
      continue;
    }

    uint64_t start = fields[origIndex].toLongLong();
    uint64_t end = fields[destIndex].toLongLong();

    if ((!dynIndToInd.contains(start)) || (!dynIndToInd.contains(end))) {
      if (saveNoAvailableNodesDemand) {
        if (!dynIndToInd.contains(start)) {
          noAvailableNodesDemand.insert(start);
        }
        if (!dynIndToInd.contains(end)) {
          noAvailableNodesDemand.insert(end);
        }
      }
      continue;
    }
    int numPeople = fields[numPeopleIndex].toInt();
    totalNumPeople += numPeople;
    demandB2018.push_back(DemandB2018(numPeople, dynIndToInd[start], dynIndToInd[end]));
  }

  // Save no available nodes to file.
  if (saveNoAvailableNodesDemand) {
    QString filename = "noAvailableNodesDemand.txt";
    saveSetToFile(noAvailableNodesDemand, filename);
  }

  std::cerr
    << "[Log] Network loaded in " << timer.elapsed() << " milliseconds with "
    << num_vertices(inRoadGraph->myRoadGraph_BI) << " vertices, "
    << num_edges(inRoadGraph->myRoadGraph_BI) << " edges, "
    << demandB2018.size() <<  " pairs of demand and "
    << totalNumPeople << " people in total." << std::endl;
}

std::string RoadGraphB2018::loadABMGraph(const std::string& networkPath, const std::shared_ptr<abm::Graph>& graph_) {
  const std::string& edgeFileName = networkPath + "edges.csv";
  std::cout << edgeFileName << " as edges file\n";

  const std::string& nodeFileName = networkPath + "nodes.csv";
  std::cout << nodeFileName << " as nodes file\n";

  const std::string& odFileName = networkPath + "od_demand.csv";
  std::cout << odFileName << " as OD file\n";
  //const bool directed = true;
  //const auto graph = std::make_shared<abm::Graph>(directed);
  auto start = high_resolution_clock::now();
  //EDGES
  //Create the graph directly from the file (don't deal with the creation of the boost graph first or any associated weights calculations)
  graph_->read_graph_osm(edgeFileName);
  //printf("# of edges: %d\n", graph_->nedges());

  //NODES
  graph_->read_vertices(nodeFileName);

  assert(graph_->amount_of_vertices_ == graph_->vertex_osm_ids_to_lc_ids_.size());

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start);
  ///////////////////////////////
  return odFileName;
}




}  // Closing namespace LC

