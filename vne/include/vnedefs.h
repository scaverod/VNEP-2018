#pragma once

#define NDEBUG
#include <vector>
#include <limits>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/lookup_edge.hpp>
#include "../include/util.h"

typedef int Cpu;
typedef int Bandwidth;

namespace vne {
class TimeLimitExceeded : public std::exception {};
class Infeasible : public std::exception {};
}

struct Parameters {
  int timelimit_in_s;
  float ro;
  float phi;
  float tmin;
  float tmax;
  int maxiter;
  int ants;
  int hops;
  int d;
  bool onlythebestreinforces;
};

struct VertexResources {
  Cpu cpu;
  inline bool operator<(const VertexResources& res) const {
    return cpu < res.cpu;
  }
  inline bool hasEnoughResources(const VertexResources& vres) const {
    return cpu >= vres.cpu;
  }
  inline void allocate(const VertexResources& vres) {
    cpu -= vres.cpu;
  }
  inline void deallocate(const VertexResources& vres) {
    cpu += vres.cpu;
  }
};

struct EdgeResources {
  Bandwidth bandwidth;
  size_t id;

  bool operator<(const EdgeResources& other) const {
    return bandwidth < other.bandwidth;
  }
};

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS,
  VertexResources, EdgeResources> Graph;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS,
  VertexResources, EdgeResources> DirGraph;
typedef Graph::vertex_descriptor  Vertex;
typedef Graph::vertex_iterator    VertexIterator;
typedef Graph::out_edge_iterator  EdgeIterator;
typedef DirGraph::out_edge_iterator DirEdgeIterator;
typedef boost::graph_traits<Graph>::edge_iterator EdgesIterator;
typedef boost::graph_traits<DirGraph>::edge_iterator DirEdgesIterator;
typedef Graph::edge_descriptor    Edge;
typedef DirGraph::edge_descriptor DirEdge;
enum VTypes { NONE = std::numeric_limits<Vertex>::max()};

struct Instance {
  Graph& substrate;
  const Graph& virtual1;
};

class Mapping {
public:
  // virtual vertex -> substrate vertex
  std::vector<Vertex> vertex;
  // whether a substrate node was mapped
  std::vector<bool> substrate_mapped;
  std::vector<bool> virtual_mapped;
  Bandwidth cost;

  Mapping(const Mapping& mapping);

  explicit Mapping(Bandwidth _cost = 0) : cost(_cost) {}
  
  Mapping(const Instance& instance) :
      vertex(num_vertices(instance.virtual1), NONE),
      substrate_mapped(num_vertices(instance.substrate), 0),
      virtual_mapped(num_vertices(instance.virtual1), 0), cost(0),
      edge(num_vertices(instance.virtual1)) { 
    for (Vertex v = 1; v < num_vertices(instance.virtual1); ++v) {
      edge[v] = std::vector<std::vector<Edge>*>(v+1, nullptr);
      EdgeIterator e, eend;
      for (tie(e, eend) = out_edges(v, instance.virtual1); e != eend; ++e) {
        const Vertex w = target(*e, instance.virtual1);
        if(v > w)
          edge[v][w] = new std::vector<Edge>();
      }
    }
  }

  Mapping(const Graph& substrate, const Graph& virtual1) :
      vertex(num_vertices(virtual1), NONE),
      substrate_mapped(num_vertices(substrate), 0),
      virtual_mapped(num_vertices(virtual1), 0), cost(0),
      edge(num_vertices(virtual1)) { 
    for (Vertex v = 1; v < num_vertices(virtual1); ++v) {
      edge[v] = std::vector<std::vector<Edge>*>(v+1, nullptr);
      EdgeIterator e, eend;
      for (tie(e, eend) = out_edges(v, virtual1); e != eend; ++e) {
        const Vertex w = target(*e, virtual1);
        if(v > w)
          edge[v][w] = new std::vector<Edge>();
      }
    }
  }

  ~Mapping() {
    for(Vertex v = 1; v < edge.size(); ++v) {
      for(Vertex u = 0; u < edge[v].size(); u++)
        if(edge[v][u] != nullptr)
          delete edge[v][u];
    }
  }
  
  void swap(Mapping& mapping) {
    using std::swap;
    swap(vertex, mapping.vertex);
    swap(substrate_mapped, mapping.substrate_mapped);
    swap(cost, mapping.cost);
    edge.swap(mapping.edge);
  }

  Mapping& operator=(Mapping mapping) {
    swap(mapping);
    return *this;
  }

  std::vector<std::vector<Edge>*>& edge_map_ref(Vertex u) {
    return edge[u];
  }

  std::vector<Edge>* edge_map_ref(Vertex u, Vertex v) {
    if (u < v)
      std::swap(u, v);
    return edge[u][v];
  }

  std::vector<Edge>* edge_map_ref(Vertex u, Vertex v) const {
    if (u < v)
      std::swap(u, v);
    return edge[u][v];
  }

  void map_vertex(const Vertex v, const Vertex s) {
    vertex[v] = s;
    substrate_mapped[s] = virtual_mapped[v] = true;
  }
  void unmap_vertex(const Vertex v) {
    substrate_mapped[vertex[v]] = virtual_mapped[v] = false;
    vertex[v] = NONE;
  }
  void map_path(Instance& instance, const Edge vedge, std::vector<Edge>& path) {
    Bandwidth vedge_cost = instance.virtual1[vedge].bandwidth;
    cost += vedge_cost * path.size();
    for (Edge e : path) {
      instance.substrate[e].bandwidth -= vedge_cost;
      assert(instance.substrate[e].bandwidth >= 0);
    }
    *edge_map_ref(source(vedge, instance.virtual1),
        target(vedge, instance.virtual1)) = path;
  }
  void map_path_using_bandwidth(Graph& substrate, const Graph& virtual1,
      const Edge v_edge, std::vector<Edge>&& path) {
    const Bandwidth demand = virtual1[v_edge].bandwidth;
    cost += demand * path.size();
    for (Edge e : path) {
      substrate[e].bandwidth -= demand;
      assert(substrate[e].bandwidth >= 0);
    }
    *edge_map_ref(source(v_edge, virtual1),
        target(v_edge, virtual1)) = path;
  }
  void map_path(const Graph& virtual1, const Edge k, const std::vector<Edge>& path) {
    const Bandwidth demand = virtual1[k].bandwidth;
    cost += demand * path.size();
    *edge_map_ref(source(k, virtual1),
        target(k, virtual1)) = path;
  }

  void unmap_path(Instance& instance, const Edge vedge) {
    Bandwidth vedge_cost = instance.virtual1[vedge].bandwidth;
    auto* path = edge_map_ref(source(vedge, instance.virtual1),
        target(vedge, instance.virtual1));
    cost -= vedge_cost * path->size();
    for (const Edge e : *path) 
      instance.substrate[e].bandwidth += vedge_cost;
    path->clear();
  }
  void unmap_all_paths(Instance& instance) {
    for (Vertex v = 1; v < edge.size(); ++v)
      for (Vertex u = 0; u < edge[v].size(); u++)
        if (edge[v][u] != nullptr)
          unmap_path(instance, boost::lookup_edge(v, u, instance.virtual1).first);
  }
private:  
  // virtual link -> substrate links
  std::vector<std::vector<std::vector<Edge>*>> edge;
};

void read_file(const char* filename, Graph* g);
void print(const Graph& g, const char prefix);
void print(const Graph& s, Graph& v);
void print(const Mapping& mapping);
int sumCpu(const Graph& g);
Bandwidth sum_bandwidth(const Graph& g);
Bandwidth out_bandwidth(const Graph& g, const Vertex v);
Bandwidth maxBandwidth(const Graph& g);
bool verify(const Mapping& mapping, Graph substrate, const Graph& virtual1);
