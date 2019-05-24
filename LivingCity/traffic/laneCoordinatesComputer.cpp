#include "./laneCoordinatesComputer.h"


namespace LC {

double D = 1.0;

void LaneCoordinatesComputer::computeMaxAmountOfLanes(const uint edgeIdx) {
  maxLanesInEdge_ = std::max(
    maxLanesInEdge_,
    static_cast<uint>(edgesData_.at(edgeIdx).numLines)
  );
}

void LaneCoordinatesComputer::computeEdgeInterface(
  const uint edgeIdx,
  const Direction direction)
{
  // If the edge's entry has already been computed return it
  if (edgesInterfaces_.find(edgeIdx) != edgesInterfaces_.end())
    return;

  const B18EdgeData edge = edgesData_.at(edgeIdx);
  const uint intersectionVertexIdx = direction == In
    ? edge.targetVertexIndex
    : edge.sourceVertexIndex;
  const uint otherVertexIdx = direction == In
    ? edge.sourceVertexIndex
    : edge.targetVertexIndex;

  // Otherwise compute it and store it
  const BoostPoint intersectionCoordinate(
    roadNetwork_.myRoadGraph_BI[intersectionVertexIdx].x,
    roadNetwork_.myRoadGraph_BI[intersectionVertexIdx].y
  );

  const BoostPoint otherExtremeCoordinate(
    roadNetwork_.myRoadGraph_BI[otherVertexIdx].x,
    roadNetwork_.myRoadGraph_BI[otherVertexIdx].y
  );

  const BoostPoint edgeOutDirection = normalizeVector(
    intersectionCoordinate,
    otherExtremeCoordinate
  );

  const BoostPoint center{
    intersectionCoordinate.x() + maxLanesInEdge_ * D * edgeOutDirection.x(),
    intersectionCoordinate.y() + maxLanesInEdge_ * D * edgeOutDirection.y(),
  };

  const BoostPoint orthogonal{
    edgeOutDirection.y(),
    -edgeOutDirection.x()
  };

  const double amountOfLanes = static_cast<double>(edgesData_.at(edgeIdx).numLines);
  const double factor = amountOfLanes * D * (direction == Out ? 1.0 : -1.0);
  const BoostPoint side{
    center.x() + factor * orthogonal.x(),
    center.y() + factor * orthogonal.y()
  };

  edgesInterfaces_.emplace(edgeIdx, EdgeInterface{center, side});
}

void LaneCoordinatesComputer::computeLaneCoordinates(
    const uint laneIdx,
    const uint edgeIdx)
{
  // If the lane coordinate has already been computed return
  if (lanesCoordinates_.find(laneIdx) != lanesCoordinates_.end())
    return;

  // Otherwise compute it and store it
  const EdgeInterface & edgeInterface = edgesInterfaces_.at(edgeIdx);
  const uint lanePositionFromLeft = laneIdx - edgeIdx;
  assert(lanePositionFromLeft < maxLanesInEdge_);
  const double amountOfLanes = static_cast<double>(edgesData_.at(edgeIdx).numLines);
  const double factor = 0.5 + D * static_cast<double>(lanePositionFromLeft);
  const BoostPoint laneCoordinate{
    edgeInterface.center_.x() + factor * edgeInterface.direction_.x(),
    edgeInterface.center_.y() + factor * edgeInterface.direction_.y()
  };

  lanesCoordinates_.insert({laneIdx, laneCoordinate});
}

LaneCoordinatesComputer::LaneCoordinatesComputer(
    const RoadGraph & roadNetwork,
    const std::vector<B18EdgeData> & edgesData,
    const std::vector<LC::Connection> & connections,
    const std::vector<LC::Intersection> & updatedIntersections) :
  roadNetwork_(roadNetwork),
  edgesData_(edgesData),
  connections_(connections),
  updatedIntersections_(updatedIntersections) {}

std::unordered_map<uint, BoostPoint> LaneCoordinatesComputer::computeLanesCoordinatesFor(
    const uint vertexIdx)
{
  const Intersection & intersection = updatedIntersections_.at(vertexIdx);

  // Reset local variables
  edgesInterfaces_.clear();
  lanesCoordinates_.clear();
  maxLanesInEdge_ = 0;

  // Compute the maximum amount of lanes in the arriving and departing edges
  for (
      uint connectionIdx = intersection.connectionGraphStart;
      connectionIdx < intersection.connectionGraphEnd;
      connectionIdx++) {
    const Connection & connection = connections_.at(connectionIdx);
    computeMaxAmountOfLanes(connection.outEdgeNumber);
    computeMaxAmountOfLanes(connection.inEdgeNumber);
  }

  // Compute edges interface
  for (
      uint connectionIdx = intersection.connectionGraphStart;
      connectionIdx < intersection.connectionGraphEnd;
      connectionIdx++) {
    // TODO(ffigari): Implement this
    const Connection & connection = connections_.at(connectionIdx);
    computeEdgeInterface(
      connection.inEdgeNumber,
      In);
    computeEdgeInterface(
      connection.outEdgeNumber,
      Out);
  }

  // Compute lanes coordinates
  for (
      uint connectionIdx = intersection.connectionGraphStart;
      connectionIdx < intersection.connectionGraphEnd;
      connectionIdx++) {
    const Connection & connection = connections_.at(connectionIdx);
    computeLaneCoordinates(
      connection.inLaneNumber,
      connection.inEdgeNumber);
    computeLaneCoordinates(
      connection.outLaneNumber,
      connection.outEdgeNumber);
  }

  return lanesCoordinates_;
}


}

