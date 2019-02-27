#include "graph.h"

// Add edge
inline void abm::Graph::add_edge(
    abm::graph::vertex_t vertex1, abm::graph::vertex_t vertex2,
    abm::graph::weight_t weight = 1,
    abm::graph::vertex_t edgeid =
        std::numeric_limits<abm::graph::vertex_t>::max()) {
  // Create a map of vertices
  if (vertices_.find(vertex1) == vertices_.end())
    vertices_[vertex1] = vertices_.size();
  if (vertices_.find(vertex2) == vertices_.end())
    vertices_[vertex2] = vertices_.size();

  if (!this->directed_)
    if (vertex1 > vertex2) std::swap(vertex1, vertex2);

  // Create an edge
  auto edge = std::make_shared<Graph::Edge>(
      std::make_pair(std::make_pair(vertex1, vertex2), weight));
  edges_[std::make_tuple(vertex1, vertex2)] = edge;

  // Add edge id
  if (edgeid == std::numeric_limits<abm::graph::vertex_t>::max()) {
    edge_ids_[std::make_tuple(vertex1, vertex2)] = this->edgeid_;
    // Add edge cost
    edge_costs_[this->edgeid_] = weight;
    this->edgeid_ += 1;
  } else {
    edge_ids_[std::make_tuple(vertex1, vertex2)] = edgeid;
    // Add edge cost
    edge_costs_[edgeid] = weight;
  }

  // Vertex 1
  auto vertex1_edges = vertex_edges_[vertex1];
  vertex1_edges.emplace_back(edge);
  vertex_edges_[vertex1] =
      std::vector<std::shared_ptr<Graph::Edge>>(vertex1_edges);

  if (!this->directed_) {
    // Vertex 2
    auto vertex2_edges = vertex_edges_[vertex2];
    vertex2_edges.emplace_back(edge);
    vertex_edges_[vertex2] =
        std::vector<std::shared_ptr<Graph::Edge>>(vertex2_edges);
  }
}

// Update edge
void abm::Graph::update_edge(abm::graph::vertex_t vertex1,
                             abm::graph::vertex_t vertex2,
                             abm::graph::weight_t weight) {
  // Get pointer to specified edge connecting vertex 1 and 2
  auto edge = edges_.at(std::make_tuple(vertex1, vertex2));
  // Update edge weight
  edge->second = weight;
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

// Read MatrixMarket graph file format
bool abm::Graph::read_graph_matrix_market(const std::string& filename) {
  bool status = true;
  try {
    std::fstream file;
    file.open(filename.c_str(), std::ios::in);
    if (file.is_open() && file.good()) {
      // Line
      std::string line;
      bool header = true;
      double ignore;
      while (std::getline(file, line)) {
        std::istringstream istream(line);
        int v1, v2;
        double weight;
        unsigned nvertices;
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('%') == std::string::npos) && (line != "")) {
          if (header) {
            // Ignore header
            istream >> nvertices;
            while (istream.good()) istream >> ignore;
            header = false;
            this->assign_nvertices(nvertices + 1);
          }
          while (istream.good()) {
            // Read vertices edges and weights
            istream >> v1 >> v2 >> weight;
            this->add_edge(v1, v2, weight);
          }
        }
      }
      std::cout << "Graph summary #edges: " << this->edges_.size()
                << " #vertices: " << this->nvertices_ << "\n";
    } else {
      throw std::runtime_error("Input file not found");
    }
  } catch (std::exception& exception) {
    std::cout << "Read matrix market file: " << exception.what() << "\n";
    status = false;
  }
  return status;
}

// Read MatrixMarket graph file format
bool abm::Graph::read_graph_osm(const std::string& filename) {
  bool status = true;
  try {
    io::CSVReader<4> in(filename);
    in.read_header(io::ignore_extra_column, "uniqueid", "u", "v", "length");
    abm::graph::vertex_t edgeid, v1, v2;
    abm::graph::weight_t weight;
    abm::graph::vertex_t nvertices = 0;
    while (in.read_row(edgeid, v1, v2, weight)) {
      this->add_edge(v1, v2, weight, edgeid);
      ++nvertices;
    }
    this->assign_nvertices(nvertices);
    std::cout << "Graph summary #edges: " << this->edges_.size()
              << " #vertices: " << this->nvertices_ << "\n";

  } catch (std::exception& exception) {
    std::cout << "Read OSM file: " << exception.what() << "\n";
    status = false;
  }

  return status;
}

void abm::Graph::generate_simple_graph() {
  this->assign_nvertices(7);
  // set up a simple graph
  this->add_edge(1, 2, 1.5);
  this->add_edge(1, 3, 9.1);
  this->add_edge(1, 6, 14.3);
  this->add_edge(2, 3, 15.9);
  this->add_edge(2, 4, 5.5);
  this->add_edge(3, 1, 5.6);
  this->add_edge(3, 4, 11.6);
  this->add_edge(3, 6, 2.4);
  this->add_edge(4, 3, 0.2);
  this->add_edge(4, 5, 6.2);
  this->add_edge(5, 6, 9.7);
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
      const abm::graph::weight_t weight = edge->second;

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
    route_edges.emplace_back(-1);
  }
  return route_edges;
}

// Determine cost of path
abm::graph::weight_t abm::Graph::path_cost(
    const std::vector<std::array<abm::graph::vertex_t, 2>>& path) {
  abm::graph::weight_t cost = 0.;
  for (const auto& vertices : path)
    cost += (edges_.at(std::make_tuple(vertices[0], vertices[1])))->second;
  return cost;
}

// Determine cost of path
abm::graph::weight_t abm::Graph::path_cost(
    const std::vector<abm::graph::vertex_t>& path) {
  abm::graph::weight_t cost = 0.;
  for (const auto& edge : path) cost += edge_costs_.at(edge);
  return cost;
}
