#include "graph.h"

// Add edge
inline void abm::Graph::add_edge(
    abm::graph::vertex_t source_osm_id,
    abm::graph::vertex_t target_osm_id,
    std::vector<float> edge_values,
    abm::graph::vertex_t edgeid = std::numeric_limits<abm::graph::vertex_t>::max()) {
	abm::graph::weight_t weight = edge_values[0];
  // Create a map of vertices
  if (vertices_.find(source_osm_id) == vertices_.end())
    vertices_[source_osm_id] = vertices_.size();
  if (vertices_.find(target_osm_id) == vertices_.end())
    vertices_[target_osm_id] = vertices_.size();

  if (!this->directed_)
    if (source_osm_id > target_osm_id) std::swap(source_osm_id, target_osm_id);

  // Create an edge
  auto edge = std::make_shared<Graph::Edge>(
      std::make_pair(std::make_pair(source_osm_id, target_osm_id), edge_values));
  edges_[std::make_tuple(source_osm_id, target_osm_id)] = edge;

  // Add edge id
  if (edgeid == std::numeric_limits<abm::graph::vertex_t>::max()) {
    edge_ids_[std::make_tuple(source_osm_id, target_osm_id)] = this->edgeid_;
    // Add edge cost
    edge_costs_[this->edgeid_] = weight;
    this->edgeid_ += 1;
  } else {
    edge_ids_[std::make_tuple(source_osm_id, target_osm_id)] = edgeid;
    edge_ids_to_vertices[edgeid] = std::make_tuple(source_osm_id, target_osm_id);
    // Add edge cost
    edge_costs_[edgeid] = weight;
  }

  // Vertex 1
  auto vertex1_edges = vertex_edges_[source_osm_id];
  vertex1_edges.emplace_back(edge);
  vertex_edges_[source_osm_id] =
      std::vector<std::shared_ptr<Graph::Edge>>(vertex1_edges);
  //out edges vertex 1
  auto vertex1_out_edges = vertex_out_edges_[source_osm_id];
  vertex1_out_edges.emplace_back(edge);
  vertex_out_edges_[source_osm_id] =
      std::vector<std::shared_ptr<Graph::Edge>>(vertex1_out_edges);
  //in edges vertex 2
  auto vertex2_in_edges = vertex_in_edges_[target_osm_id];
  vertex2_in_edges.emplace_back(edge);
  vertex_in_edges_[target_osm_id] =
      std::vector<std::shared_ptr<Graph::Edge>>(vertex2_in_edges);
  if (!this->directed_) {
    // Vertex 2
    auto vertex2_edges = vertex_edges_[target_osm_id];
    vertex2_edges.emplace_back(edge);
    vertex_edges_[target_osm_id] =
        std::vector<std::shared_ptr<Graph::Edge>>(vertex2_edges);

  //out edges vertex 2
  auto vertex2_out_edges = vertex_out_edges_[target_osm_id];
  vertex2_out_edges.emplace_back(edge);
  vertex_out_edges_[target_osm_id] =
      std::vector<std::shared_ptr<Graph::Edge>>(vertex2_out_edges);
  //in edges vertex 1
  auto vertex1_in_edges = vertex_in_edges_[source_osm_id];
  vertex1_in_edges.emplace_back(edge);
  vertex_in_edges_[source_osm_id] =
      std::vector<std::shared_ptr<Graph::Edge>>(vertex1_in_edges);
  }
}

// Update edge
void abm::Graph::update_edge(abm::graph::vertex_t vertex1,
                             abm::graph::vertex_t vertex2,
                             abm::graph::weight_t weight) {
  // Get pointer to specified edge connecting vertex 1 and 2
  auto edge = edges_.at(std::make_tuple(vertex1, vertex2));
  // Update edge weight
  edge->second[0] = weight;
  //std::cout << "weight = " << std::get<1>(x) << "\n";
}

// Remove edge
void abm::Graph::remove_edge(abm::graph::vertex_t vertex1,
                             abm::graph::vertex_t vertex2) {
  auto edge = edges_[std::make_tuple(vertex1, vertex2)];
  edges_.erase(edges_.find(std::make_tuple(vertex1, vertex2)));

  if (edge_ids_.size() > 0 &&
      edge_ids_.find(std::make_tuple(vertex1, vertex2)) != edge_ids_.end())
    edge_ids_.erase(edge_ids_.find(std::make_tuple(vertex1, vertex2)));

  auto v1edge = vertex_edges_.at(vertex1);
  auto v2edge = vertex_edges_.at(vertex2);

  v1edge.erase(std::remove(v1edge.begin(), v1edge.end(), edge));
  v2edge.erase(std::remove(v2edge.begin(), v2edge.end(), edge));

  vertex_edges_[vertex1] = v1edge;
  vertex_edges_[vertex2] = v2edge;
}

bool abm::Graph::read_graph_osm(const std::string& filename) {
  bool status = true;
  try {
    csvio::CSVReader<6> in(filename);
    abm::graph::vertex_t edge_osm_id, source_osm_id, target_osm_id;
    std::vector<float> edge_values(3);
    float length, lanes, speed_mph;
    abm::graph::vertex_t edge_lc_id = 0;
    in.read_header(csvio::ignore_no_column, "uniqueid", "u", "v", "length", "lanes", "speed_mph");
    while (in.read_row(edge_osm_id, source_osm_id, target_osm_id, length, lanes, speed_mph)) {
      edge_values[0] = length;
      edge_values[1] = lanes;
      edge_values[2] = speed_mph;

      //Don't add if there is already an edge with the same vertices
      if (edges_.find(std::make_pair(source_osm_id, target_osm_id)) == edges_.end()) {
        this->add_edge(source_osm_id, target_osm_id, edge_values, edge_osm_id);
      }

      //map edge vertex ids to smaller values
      edge_vertex_map_[source_osm_id] = edge_lc_id;
      ++edge_lc_id;
    }
    
    std::cout << "# of edges: " << this->edges_.size() << "\n";

  } catch (std::exception& exception) {
    std::cout << "Read OSM file: " << exception.what() << "\n";
    status = false;
  }

  return status;
}

void abm::Graph::read_vertices(const std::string& filename) {
  vertex_osm_ids_to_lc_ids_.clear();
  vertex_lc_ids_to_osm_ids_.clear();
  vertices_data_.clear();
  amount_of_vertices_ = 0;

  QVector2D minBox(FLT_MAX, FLT_MAX);
  QVector2D maxBox(-FLT_MAX, -FLT_MAX);
  float scale = 1.0f;
  float sqSideSz = std::max<float>(maxBox.x() - minBox.x(), maxBox.y() - minBox.y()) * scale * 0.5f; // half side
  QVector3D centerV(-minBox.x(), -minBox.y(), 0);
  QVector3D centerAfterSc(-sqSideSz, -sqSideSz, 0);
  try {
    abm::graph::vertex_t osm_id;
    float lat, lon;
    std::string osm_string_id;

    csvio::CSVReader<4> in(filename);
    abm::graph::vertex_t living_city_id = 0;
    in.read_header(csvio::ignore_extra_column, "osmid", "x", "y", "highway");
    while (in.read_row(osm_id, lat, lon, osm_string_id)) {
      vertex_osm_ids_to_lc_ids_[osm_id] = living_city_id;
      vertex_lc_ids_to_osm_ids_[living_city_id] = osm_id;
      ++living_city_id;

      QVector3D position(lat, lon, 0);
      position += centerV;
      position *= scale;
      position += centerAfterSc;
      position.setX(position.x() * -1.0f);
      vertices_data_[osm_id] = position;

      vertex_OSM_type_[osm_id] = mapStringToOSMConstant(osm_string_id);
    }
    amount_of_vertices_ = living_city_id;
  } catch (std::exception& exception) {
    std::cout
      << "abm::Graph::read_vertices -> Nodes loading failed with: "
      << exception.what()
      << std::endl;
    throw;
  }
}

// Dijktra shortest paths from src to a vertex
std::vector<abm::graph::vertex_t> abm::Graph::dijkstra(
    abm::graph::vertex_t source, abm::graph::vertex_t destination) {

  // Using lambda to compare elements.
  auto compare =
      [](std::pair<abm::graph::weight_t, abm::graph::vertex_t> left,
         std::pair<abm::graph::weight_t, abm::graph::vertex_t> right) {
        return left.first > right.first;
      };


  // Create a priority queue to store weights and vertices
  std::priority_queue<
      std::pair<abm::graph::weight_t, abm::graph::vertex_t>,
      std::vector<std::pair<abm::graph::weight_t, abm::graph::vertex_t>>,
      decltype(compare)>
      priority_queue(compare);

  // Create a vector for distances and initialize all to max
  std::vector<graph::weight_t> distances;
  distances.resize(this->vertices_.size(),
                   std::numeric_limits<abm::graph::weight_t>::max());
  // Parent array to store shortest path tree
  std::vector<graph::vertex_t> parent;
  parent.resize(this->vertices_.size(), -1);

  std::vector<abm::graph::vertex_t> path;
  if (vertices_.find(source) == vertices_.end() ||
      vertices_.find(destination) == vertices_.end())
    return path;

  // Insert source itself in priority queue & initialize its distance as 0.
  priority_queue.push(std::make_pair(0., source));
  distances[vertices_.at(source)] = 0.;

  // Looping till priority queue becomes empty (or all
  // distances are not finalized)
  while (!priority_queue.empty()) {
    // {min_weight, vertex} sorted based on weights (distance)
    abm::graph::vertex_t u = priority_queue.top().second;
    priority_queue.pop();

    // Break if destination is reached
    if (u == destination) break;

    // Get all adjacent vertices of a vertex
    for (const auto& edge : vertex_edges_[u]) {
      // Get vertex label and weight of neighbours of u.
      const abm::graph::vertex_t neighbour = edge->first.second;
      const abm::graph::weight_t weight = edge->second[0];

      // Distance from source to neighbour
      // distance_u = distance to current node + weight of edge u to
      // neighbour
      const abm::graph::vertex_t uidx = vertices_.at(u);
      const abm::graph::vertex_t nidx = vertices_.at(neighbour);

      const abm::graph::weight_t distance_u = distances.at(uidx) + weight;
      // If there is shorted path to neighbour vertex through u.
      if (distances.at(nidx) > distance_u) {
        parent[nidx] = u;
        // Update distance of the vertex
        distances.at(nidx) = distance_u;
        priority_queue.push(std::make_pair(distance_u, neighbour));
      }
    }
  }

  path.emplace_back(destination);
  // Iterate until source has been reached
  while (destination != source && destination != -1) {
    destination = parent.at(vertices_.at(destination));
    if (destination != source && destination != -1)
      path.emplace_back(destination);
  }
  path.emplace_back(source);
  // Reverse to arrange path from source to destination
  std::reverse(std::begin(path), std::end(path));

  return path;
}

// Dijktra shortest paths from src to a vertex return vertices
std::vector<std::array<abm::graph::vertex_t, 2>> abm::Graph::dijkstra_vertices(
    abm::graph::vertex_t source, abm::graph::vertex_t destination) {

  const auto path = this->dijkstra(source, destination);

  std::vector<std::array<abm::graph::vertex_t, 2>> route_vertices;
  if (path.size() > 0) {
    for (auto itr = path.begin(); itr != path.end() - 1; ++itr) {
      auto nitr = itr + 1;
      if (itr != path.end()) {
        std::array<abm::graph::vertex_t, 2> edges = {
            static_cast<abm::graph::vertex_t>(*itr),
            static_cast<abm::graph::vertex_t>(*nitr)};
        route_vertices.emplace_back(edges);
      }
    }
  }
  return route_vertices;
}

// Dijktra shortest paths from src to a vertex return vertices
std::vector<abm::graph::vertex_t> abm::Graph::dijkstra_vertices_ual(
    abm::graph::vertex_t source, abm::graph::vertex_t destination) {

  const auto path = this->dijkstra(source, destination);

  std::vector<abm::graph::vertex_t> route_vertices;
  if (path.size() > 0) {
    for (auto itr = path.begin(); itr != path.end(); ++itr) {
      if (itr != path.end())
        route_vertices.emplace_back(static_cast<abm::graph::vertex_t>(*itr));
    }
    route_vertices.emplace_back(-1);
  }
  return route_vertices;
}

// Dijktra shortest paths from src to a vertex return edges
std::vector<abm::graph::vertex_t> abm::Graph::dijkstra_edges(
    abm::graph::vertex_t source, abm::graph::vertex_t destination) {

  const auto path = this->dijkstra(source, destination);
  //printf("path size = %d\n", path.size());

  std::vector<abm::graph::vertex_t> route_edges;
  if (path.size() > 0) {
    // Reverse to arrange from source to destination
    for (auto itr = path.begin(); itr != path.end() - 1; ++itr) {
      auto nitr = itr + 1;
      if (itr != path.end()) {
        auto map_itr = edge_ids_.find(std::make_tuple(*itr, *nitr));
        if (map_itr != edge_ids_.end())
          route_edges.emplace_back((*map_itr).second);
      }
    }
  }
  route_edges.emplace_back(-1);
  return route_edges;
}

// Determine cost of path
abm::graph::weight_t abm::Graph::path_cost(
    const std::vector<std::array<abm::graph::vertex_t, 2>>& path) {
  abm::graph::weight_t cost = 0.;
  for (const auto& vertices : path)
    cost += (edges_.at(std::make_tuple(vertices[0], vertices[1])))->second[0];
  return cost;
}

// Determine cost of path
abm::graph::weight_t abm::Graph::path_cost(
    const std::vector<abm::graph::vertex_t>& path) {
  abm::graph::weight_t cost = 0.;
  for (const auto& edge : path) cost += edge_costs_.at(edge);
  return cost;
}
