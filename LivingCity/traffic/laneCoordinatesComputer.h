#include <vector>
#include <unordered_map>

#include "./RoadGraph/roadGraph.h"
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

class LaneCoordinatesComputer
{

  public:
    LaneCoordinatesComputer(
      const RoadGraph & roadNetwork,
      const std::vector<B18EdgeData> & edgesData,
      const std::vector<LC::Connection> & connections,
      const std::vector<LC::Intersection> & updatedIntersections
    );
    std::unordered_map<uint, BoostPoint> computeLanesCoordinatesFor(const uint vertexIdx);

  private:
    std::unordered_map<uint, EdgeInterface> edgesInterfaces_;
    std::unordered_map<uint, BoostPoint> lanesCoordinates_;
    uint maxLanesInEdge_;

    const RoadGraph & roadNetwork_;
    const std::vector<B18EdgeData> & edgesData_;
    const std::vector<LC::Connection> & connections_;
    const std::vector<LC::Intersection> & updatedIntersections_;

    void computeMaxAmountOfLanes(const uint edgeIdx);
    void computeEdgeInterface(const uint edgeIdx, const Direction direction);
    void computeLaneCoordinates(const uint laneIdx, const uint edgeIdx);
};


}
