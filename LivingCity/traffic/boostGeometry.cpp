#include "./boostGeometry.h"

BoostPoint normalizeVector(const BoostPoint & source, const BoostPoint & target) {
  const double deltaX = target.x() - source.x();
  const double deltaY = target.y() - source.y();
  const double directionNorm = std::pow(std::pow(deltaX, 2) + std::pow(deltaY, 2), 0.5);
  return {
    deltaX / directionNorm,
    deltaY / directionNorm
  };
}

BoostPoint orthogonalize(const BoostPoint & source) {
  return {
    source.y(),
    -source.x()
  };
}

