#ifndef BACKTRACK_H_
#define BACKTRACK_H_

#include <iostream>
#include <queue>
#include <functional>

#include "../include/util.h"
#include "../include/vnedefs.h"


namespace vne{
using std::cout;
using std::endl;
using std::numeric_limits;
using std::list;
using std::chrono::seconds;
using std::vector;


struct Bundle {
  Mapping& current;
  Mapping& best;
  // lb[v] = lower bound on the solution for v to n
  Bandwidth upper_bound;
  int timelimit_in_s;
  Timer<seconds> timer;

  vector<Vertex> substrate_nodes;
  vector<Vertex> virtual_nodes;
  vector<Bandwidth> lower_bound;
  vector<vector<Bandwidth>> min_hop;
};

class Backtrack {
public:
  static const bool sort_vertices = false;
  static const bool calc_initial_solution = false;
  static const bool per_link_lb_calc = false; 
  static void precalc(const Instance&, Bundle&);
  static void initial_solution(const Instance&, Bundle&);
  static void sort(const Instance&, Bundle&);
};

class HopBacktrack : public Backtrack {
public:
  static const bool per_link_lb_calc = true;
  static void precalc(const Instance&, Bundle&);
};

class GreedyInitialBacktrack : public HopBacktrack {
public:
  static const bool calc_initial_solution = true;
  static void initial_solution(Instance&, Bundle&) throw (Infeasible);
};

// Sort
class NoSort {
public:
  static const bool sort_vertices = false;
  static void sort(const Instance&, Bundle&);
};

template <class VertexComparator>
class Sort {
public:
  static const bool sort_vertices = true;
  static void sort(const Instance&, Bundle&);
};

class MoreCpu {
public:
  MoreCpu(const Graph& graph) : graph_(graph) { }
  bool operator()(const Vertex v, const Vertex w) const {
    return graph_[v].cpu > graph_[w].cpu;
  }
private:
  const Graph& graph_;
};

class GreaterDegree {
public:
  GreaterDegree(const Graph& graph) : graph_(graph) { }
  bool operator()(const Vertex v, const Vertex w) const {
    return degree(v, graph_) > degree(w, graph_);
  }
private:
  const Graph& graph_;
};

template <class VertexComparator>
void Sort<VertexComparator>::sort(const Instance& instance, Bundle& bundle) {
  std::sort(bundle.virtual_nodes.begin(), bundle.virtual_nodes.end(),
      VertexComparator(instance.virtual1));
  std::sort(bundle.substrate_nodes.begin(), bundle.substrate_nodes.end(),
      VertexComparator(instance.substrate));
}

struct VWPath {
  Vertex w;
  std::vector<Edge> path;
  std::vector<bool> visited;
  VWPath(Vertex _w, size_t subsize) : w(_w),
      visited(subsize) { 
    visited[w] = true;
  }
};

template <class BACKTRACK>
bool backtrack(const size_t, Instance&, Bundle&) throw(TimeLimitExceeded);

// map all virtual edges from v to the vertices in to_map
// for each possible map, call backtrack
template <class BACKTRACK> 
bool backtrack_paths(const size_t vindex, const list<Edge>::iterator w,
    const list<Edge>& to_map, Instance& instance,
    Bundle& bundle) {
  // base case
  if (w == to_map.end()) {
    return backtrack<BACKTRACK>(vindex+1, instance, bundle);
  }
  const Vertex v {bundle.virtual_nodes[vindex]};
  const Vertex source {bundle.current.vertex[v]},
               vdest {target(*w, instance.virtual1)},
               dest {bundle.current.vertex[vdest]};
  // lower bound
  if (BACKTRACK::per_link_lb_calc) {
    Bandwidth lb {bundle.current.cost};
    for (auto e = w; e != to_map.end(); ++e) {
      const Vertex edest {bundle.current.vertex[target(*e, instance.substrate)]};
      if (source > edest)
        lb += bundle.min_hop[source][edest] *
            instance.virtual1[*e].bandwidth;
      else
        lb += bundle.min_hop[edest][source] *
            instance.virtual1[*e].bandwidth;
    }
    if (lb >= bundle.upper_bound) { 
      return false;
    }
  }
  Bandwidth virtual_edge_cost {instance.virtual1[*w].bandwidth};
  bool valid_mapping_found = false;
  std::queue<VWPath> Q;
  Q.push(VWPath(source, num_vertices(instance.substrate)));
  while (!Q.empty()) {
    if (Q.front().w == dest) {
      // maps virtual link (v,w) to substrate path
      bundle.current.map_path(instance, *w, Q.front().path);
      if (bundle.current.cost < bundle.upper_bound &&
          backtrack_paths<BACKTRACK>(vindex, std::next(w), to_map, instance,
            bundle))
        valid_mapping_found = true;
      bundle.current.unmap_path(instance, *w);
    } else {
      Q.front().visited[Q.front().w] = true;
      EdgeIterator e, eend;
      for (tie(e, eend) = out_edges(Q.front().w, instance.substrate);
          e != eend; ++e) {
        if (instance.substrate[*e].bandwidth >= virtual_edge_cost) {
          const Vertex r = target(*e, instance.substrate);
          if (!Q.front().visited[r]) {
            VWPath rpath = Q.front();
            rpath.w = r;
            rpath.path.push_back(*e);
            Q.push(rpath);
          }
        }
      }
    }
    Q.pop();
  }
  return valid_mapping_found;
}

// backtrack over all solutions, returns whether a valid mapping was found
template <class BACKTRACK>
bool backtrack(const size_t vindex, Instance& instance,
               Bundle& bundle) throw(TimeLimitExceeded) {
  // time limit exceeded
  if (bundle.timer.reachedTimeLimit())
    throw TimeLimitExceeded();

  // base case
  if (vindex >= num_vertices(instance.virtual1)) {
    if (bundle.current.cost < bundle.best.cost) {
      bundle.upper_bound = bundle.current.cost;
      bundle.best = bundle.current;
    }
    return true;
  }
  const Vertex v {bundle.virtual_nodes[vindex]};
  // test simple lower bound
  if (bundle.lower_bound[vindex] + bundle.current.cost >=
      bundle.upper_bound)
    return false; 
  // map v to every substrate node s that is not used
  bool valid_mapping_found = false;
  for (Vertex s : bundle.substrate_nodes) {
    // bool s_can_be_mapped_to_v
    if (!bundle.current.substrate_mapped[s] &&
        instance.substrate[s].hasEnoughResources(instance.virtual1[v])) {
      // cout << "map " << v << " to " << s << endl;
      bundle.current.map_vertex(v, s);
      // make a list of all unmapped edges
      EdgeIterator e, eend;
      list<Edge> to_map;
      for (tie(e, eend) = out_edges(v, instance.virtual1); e != eend; ++e) {
        if (bundle.current.virtual_mapped[target(*e, instance.virtual1)]) {
          to_map.push_back(*e);
        }
      }
      if (backtrack_paths<BACKTRACK>(vindex, to_map.begin(), to_map, instance, bundle))
        valid_mapping_found = true;
      bundle.current.unmap_vertex(v);
    }
  }
  return valid_mapping_found;
}

template <class BACKTRACK, class SORT = NoSort>
bool backtrack(Graph& substrate, const Graph& virtual1, const Parameters param,
    Mapping& out_mapping) {
  Instance instance {substrate, virtual1};
  Mapping current(instance);
  out_mapping.cost = numeric_limits<Bandwidth>::max();
  Timer<seconds> timer;
  timer.start();
  Bundle bundle{
    current,
    out_mapping,
    /* upper_bound */ std::numeric_limits<Bandwidth>::max(),
    param.timelimit_in_s,
    timer};
  bundle.virtual_nodes = vector<Vertex>(num_vertices(instance.virtual1));
  std::iota(bundle.virtual_nodes.begin(), bundle.virtual_nodes.end(), 0);
  bundle.substrate_nodes = vector<Vertex>(num_vertices(instance.substrate));
  std::iota(bundle.substrate_nodes.begin(), bundle.substrate_nodes.end(), 0);
  try {
    if (SORT::sort_vertices)
      SORT::sort(instance, bundle);
    if (BACKTRACK::calc_initial_solution)
      BACKTRACK::initial_solution(instance, bundle);
    BACKTRACK::precalc(instance, bundle);
    return backtrack<BACKTRACK>(0, instance, bundle);
  } catch (Infeasible e) {
    cout << "infeasible" << endl;
    return false;
  } catch (TimeLimitExceeded e) {
    bool aSolutionWasFound = out_mapping.cost < numeric_limits<Bandwidth>::max();
    if (aSolutionWasFound)
      cout << "BestCost:" << out_mapping.cost << endl;
    return false;
  }
}
}
#endif
