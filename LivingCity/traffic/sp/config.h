#ifndef _ABM_CONFIG_H_
#define _ABM_CONFIG_H_

#include "../types_definitions.h"

namespace abm {
namespace graph {
//! Vertex id type
using vertex_t = osm_id_type;
//! Weight type, that can be added with +
using weight_t = double;
}  // namespace graph
}  // namespace abm

#endif  // _ABM_CONFIG_H_
