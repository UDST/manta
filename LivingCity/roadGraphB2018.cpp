/************************************************************************************************
*		@desc Class to load the B2018 road graph
*		@author igarciad
************************************************************************************************/

#include "roadGraphB2018.h"

#include "Geometry/client_geometry.h"
#include "LC_UrbanMain.h"
#include "LC_GLWidget3D.h"

#include "global.h"
#include "bTraffic/bTrafficIntersection.h"

#include <QHash>
#include <QVector2D>

namespace LC {

  void updateMinMax2(QVector2D& newPoint, QVector2D& minBox, QVector2D& maxBox) {
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

  // WGS-84 geodetic constants
  const double a = 6378137.0;         // WGS-84 Earth semimajor axis (m)

  const double b = 6356752.314245;     // Derived Earth semiminor axis (m)
  const double f = (a - b) / a;           // Ellipsoid Flatness
  const double f_inv = 1.0 / f;       // Inverse flattening

  const double a_sq = a * a;
  const double b_sq = b * b;
  const double e_sq = f * (2 - f);    // Square of Eccentricity

# define M_PI           3.14159265358979323846  /* pi */

  static double DegreesToRadians(double degrees) {
    return M_PI / 180.0 * degrees;
  }

  static double RadiansToDegrees(double radians) {
    return 180.0 / M_PI * radians;
  }

  // Converts WGS-84 Geodetic point (lat, lon, h) to the 
  // Earth-Centered Earth-Fixed (ECEF) coordinates (x, y, z).
  static void GeodeticToEcef(double lat, double lon, double h,
    double& x, double& y, double& z) {
    // Convert to radians in notation consistent with the paper:
    double lambda = DegreesToRadians(lat);
    double phi = DegreesToRadians(lon);
    double s = sin(lambda);
    double N = a / sqrt(1 - e_sq * s * s);

    double sin_lambda = sin(lambda);
    double cos_lambda = cos(lambda);
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);

    x = (h + N) * cos_lambda * cos_phi;
    y = (h + N) * cos_lambda * sin_phi;
    z = (h + (1 - e_sq) * N) * sin_lambda;
  }

  // Converts the Earth-Centered Earth-Fixed (ECEF) coordinates (x, y, z) to 
  // East-North-Up coordinates in a Local Tangent Plane that is centered at the 
  // (WGS-84) Geodetic point (lat0, lon0, h0).
  static void EcefToEnu(double x, double y, double z,
    double lat0, double lon0, double h0,
    double& xEast, double& yNorth, double& zUp) {
    // Convert to radians in notation consistent with the paper:
    double lambda = DegreesToRadians(lat0);
    double phi = DegreesToRadians(lon0);
    double s = sin(lambda);
    double N = a / sqrt(1 - e_sq * s * s);

    double sin_lambda = sin(lambda);
    double cos_lambda = cos(lambda);
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);

    double x0 = (h0 + N) * cos_lambda * cos_phi;
    double y0 = (h0 + N) * cos_lambda * sin_phi;
    double z0 = (h0 + (1 - e_sq) * N) * sin_lambda;

    double xd, yd, zd;
    xd = x - x0;
    yd = y - y0;
    zd = z - z0;

    // This is the matrix multiplication
    xEast = -sin_phi * xd + cos_phi * yd;
    yNorth = -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd;
    zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd;
  }

  // Converts the geodetic WGS-84 coordinated (lat, lon, h) to 
  // East-North-Up coordinates in a Local Tangent Plane that is centered at the 
  // (WGS-84) Geodetic point (lat0, lon0, h0).
  static void GeodeticToEnu(double lat, double lon, double h,
    double lat0, double lon0, double h0,
    double& xEast, double& yNorth, double& zUp) {
    double x, y, z;
    GeodeticToEcef(lat, lon, h, x, y, z);
    EcefToEnu(x, y, z, lat0, lon0, h0, xEast, yNorth, zUp);
  }

  //////////////////////////////////////////////////////////

  void RoadGraphB2018::loadB2018RoadGraph(RoadGraph& inRoadGraph, LCGLWidget3D* glWidget3D) {

    printf(">>loadB2018RoadGraph\n");
    printf(">>Remove\n");
    inRoadGraph.myRoadGraph.clear();
    inRoadGraph.myRoadGraph_BI.clear();
    glWidget3D->cg.geoZone.blocks.clear();
    glWidget3D->vboRenderManager.removeAllStreetElementName("tree");
    glWidget3D->vboRenderManager.removeAllStreetElementName("streetLamp");
    printf("<<Remove\n");


    printf("1 loadB2018RoadGraph\n");
    QString fileName = "berkeley_2018/bay_area_full_strongly_nodes.csv";
    //QString fileName="data/Dynameq/smallTestNet_base.dqt";
    QFile baseFile(fileName); // Create a file handle for the file named

    QString line;
    if (!baseFile.open(QIODevice::ReadOnly | QIODevice::Text)) { // Open the file
      printf("Can't open file '%s'\n", fileName.toUtf8().constData());
      return;
    }
    QTextStream stream(&baseFile); // Set the stream to read from myFile
    /////////////////////////////////////////////////
    // READ FILE
    QTime timer;
    timer.start();
    QVector2D minBox(FLT_MAX, FLT_MAX);
    QVector2D maxBox(-FLT_MAX, -FLT_MAX);
    QHash<int, QVector2D> osmidToVertexLoc;

    stream.readLine();
    while (!stream.atEnd()) {
      line = stream.readLine();
      QStringList fields = line.split(',', QString::SkipEmptyParts);
      if (fields.size() < 3) {
        qDebug() << "ERROR line " << line << " --> SKIP";
        continue;
      }
      float x = fields[1].toFloat();
      float y = fields[2].toFloat();
      //qDebug() << "x " << x << " y " << y;
      osmidToVertexLoc[fields[0].toInt()] = QVector2D(x, y);
      updateMinMax2(QVector2D(x, y), minBox, maxBox);
      //printf("MinBox %f %f MaxBox %f %f--> %f %f\n", minBox.x(), minBox.y(), maxBox.x(), maxBox.y(), maxBox.x() - minBox.x(), maxBox.y() - minBox.y());
    }
    printf("Degrees --> MinBox %f %f MaxBox %f %f--> %f %f\n", minBox.x(), minBox.y(), maxBox.x(), maxBox.y(), maxBox.x() - minBox.x(), maxBox.y() - minBox.y());

    // Update coordenades to East-North-Up coordinates;

    double xEast, yNorth, zUp;
    const float h0 = a; //  6371000.0f;
    const float lat0 = (maxBox.x() + minBox.x())*0.5f;
    const float lon0 = (maxBox.y() + minBox.y())*0.5f;
    minBox = QVector2D(FLT_MAX, FLT_MAX);
    maxBox = QVector2D(-FLT_MAX, -FLT_MAX);
    QHash<int, QVector2D>::iterator i;
    for (i = osmidToVertexLoc.begin(); i != osmidToVertexLoc.end(); ++i) {

      GeodeticToEnu(i.value().x(), i.value().y(), h0,
        lat0, lon0, h0,
        xEast, yNorth, zUp);
      osmidToVertexLoc[i.key()] = QVector2D(xEast, yNorth);
      updateMinMax2(QVector2D(xEast, yNorth), minBox, maxBox);
    }
    printf("Meters? --> MinBox %f %f MaxBox %f %f--> %f %f\n", minBox.x(), minBox.y(), maxBox.x(), maxBox.y(), maxBox.x() - minBox.x(), maxBox.y() - minBox.y());

    printf("\nNodes readed in %d Nod %d\n", timer.elapsed(), osmidToVertexLoc.size());
    //printf("\nNodes readed in %d Nod %d Cen %d Link %d\n", timer.elapsed(), osmidToVertexLoc.size(), centroids.size(), links.size());
  }

}  // namespace LC