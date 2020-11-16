#ifndef _ABM_CONFIG_H_
#define _ABM_CONFIG_H_

namespace abm {
namespace graph {
//! Vertex id type
using vertex_t = long long;
using edge_id_t = long long;
//! Weight type, that can be added with +
using weight_t = double;
}  // namespace graph
}  // namespace abm

struct parameters{
    float a;
    float b;
    float T;
    float s_0;
};

#endif  // _ABM_CONFIG_H_


