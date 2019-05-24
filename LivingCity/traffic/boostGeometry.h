#ifndef BOOST_GEOMETRY__H
#define BOOST_GEOMETRY__H

#include <limits>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

typedef boost::geometry::model::d2::point_xy<double> BoostPoint;
typedef boost::geometry::model::segment<BoostPoint> BoostSegment;

BoostPoint normalizeVector(const BoostPoint & source, const BoostPoint & target);

BoostPoint orthogonalize(const BoostPoint & source);

inline std::ostream& operator<<(std::ostream& stream, const BoostPoint & point) {
  stream << "(" << point.x() << ", " << point.y() << ")";
  return stream;
}

inline bool operator==(const BoostPoint & onePoint, const BoostPoint & anotherPoint) {
  const auto epsilon = std::numeric_limits<double>::epsilon();
  return
    std::fabs(onePoint.x() - anotherPoint.x()) < epsilon
    && std::fabs(onePoint.y() - anotherPoint.y()) < epsilon;
}

inline bool operator!=(const BoostPoint & onePoint, const BoostPoint & anotherPoint) {
  return !(onePoint == anotherPoint);
}


#endif  // BOOST_GEOMETRY__H

