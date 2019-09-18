#include <vector>
#include <unordered_map>
#include <functional>

#include "./RoadGraph/roadGraph.h"
#include "./types_definitions.h"
#include "./b18EdgeData.h"
#include "./boostGeometry.h"


namespace LC {


enum Direction {
  In,
  Out
};

struct EdgeInterface
{
  const BoostPoint center_;
  const BoostPoint side_;
  const BoostPoint direction_;

  EdgeInterface(
      const BoostPoint & center,
      const BoostPoint & side) :
    center_(center),
    side_(side),
    direction_(normalizeVector(side, center)) {}
};

using CoordinatesRetriever = std::function<std::pair<double, double>(const uint vertexIdx)>;

class LaneCoordinatesComputer
{

  public:
    LaneCoordinatesComputer(
      const CoordinatesRetriever & coordinatesRetriever,
      const std::vector<B18EdgeData> & edgesData,
      const std::vector<LC::Connection> & connections,
      const std::vector<LC::Intersection> & updatedIntersections
    );
    std::unordered_map<uint, BoostPoint> computeLanesCoordinatesFor(const uint vertexIdx);

  private:
    std::unordered_map<uint, EdgeInterface> edgesInterfaces_;
    std::unordered_map<uint, BoostPoint> lanesCoordinates_;
    uint maxLanesInEdge_;

    const CoordinatesRetriever & coordinatesRetriever_;
    const std::vector<B18EdgeData> & edgesData_;
    const std::vector<LC::Connection> & connections_;
    const std::vector<LC::Intersection> & updatedIntersections_;

    void computeMaxAmountOfLanes(const uint edgeIdx);
    void computeEdgeInterface(const uint edgeIdx, const Direction direction);
    void computeLaneCoordinates(const uint laneIdx, const uint edgeIdx);
};


}
