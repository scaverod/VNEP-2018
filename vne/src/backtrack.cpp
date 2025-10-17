// Copyright 2013 Leonardo Moura
#include <cstdio>
#include <sys/times.h>
#include <cassert>
#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <limits>
#include <algorithm>
#include <exception>
#include <functional>
#include "../include/backtrack.h"
#include "../include/vnedefs.h"

namespace vne {
using std::pair;
using std::make_pair;
using std::queue;
using std::swap;

// PRECALC functions
void Backtrack::precalc(const Instance& instance, Bundle& bundle) {
  bundle.lower_bound = vector<Bandwidth>(num_vertices(instance.virtual1));

  for (size_t i = num_vertices(instance.virtual1); i-- > 0;) {
    if (i < num_vertices(instance.virtual1)-1)
      bundle.lower_bound[i] = bundle.lower_bound[i+1];
    else
      bundle.lower_bound[i] = 0;
    const Vertex v {bundle.virtual_nodes[i]};
    EdgeIterator e, eend;
    for (tie(e, eend) = out_edges(v, instance.virtual1); e != eend; ++e) {
      if (v > target(*e, instance.virtual1))
        bundle.lower_bound[i] += instance.virtual1[*e].bandwidth;
    }
  }
}

void HopBacktrack::precalc(const Instance& instance, Bundle& bundle) {
  Backtrack::precalc(instance, bundle);

  bundle.min_hop = vector<vector<Bandwidth>>(num_vertices(instance.substrate));
  std::vector<bool> visited(num_vertices(instance.substrate));
  for (int s = 1; s < num_vertices(instance.substrate); ++s) {
    bundle.min_hop[s] = vector<Bandwidth>(s); 
    std::vector<int>  distance(num_vertices(instance.substrate));
    visited.clear();
    std::queue<Vertex> Q;
    Q.push(s);
    while (!Q.empty()) {
      Vertex r = Q.front();
      Q.pop();
      EdgeIterator e, eend;
      for (tie(e, eend) = out_edges(r, instance.substrate); e != eend; ++e) {
        Vertex t = target(*e, instance.substrate);
        if (!visited[t]) {
          visited[t] = true;
          if(s > t)
              bundle.min_hop[s][t] = distance[t] = distance[r] + 1;
          Q.push(t);
        }
      }
    }
  }
}

// finds the path with least hops in instance.substrate from s to t
// saves it in out_path, if there's no path returns false
bool find_min_path(const Instance& instance, const Bandwidth demand,
    const Vertex s, const Vertex t, vector<Edge>& out_path) {
  std::vector<bool> visited(num_vertices(instance.substrate));
  std::vector<Edge> pred(num_vertices(instance.substrate));
  std::queue<Vertex> Q;
  Q.push(s);
  visited[s] = true;
  while (!Q.empty()) {
    Vertex r = Q.front();
    if (r == t) {
      while (r != s) {
        out_path.push_back(pred[r]);
        r = source(pred[r], instance.substrate);
      }
      return true;
    }
    Q.pop(); 
    EdgeIterator e, eend;
    for (tie(e, eend) = out_edges(r, instance.substrate); e != eend; ++e) {
      if (instance.substrate[*e].bandwidth < demand)
        continue;
      const Vertex y = target(*e, instance.substrate);
      if (!visited[y]) {
        visited[y] = true;
        pred[y] = *e;
        Q.push(y);
      }
    }
  }
  return false;
}

// TODO use the fact that bundle.virtual_nodes are sorted
void GreedyInitialBacktrack::initial_solution(Instance& instance,
    Bundle& bundle) throw(Infeasible) {
  typedef pair<int, Vertex> cpu_index;
  vector<cpu_index> s_nodes;
  s_nodes.reserve(num_vertices(instance.substrate));
  for (Vertex s {0}; s < num_vertices(instance.substrate); ++s) {
    s_nodes.push_back(make_pair(instance.substrate[s].cpu, s));
  }  
  std::sort(s_nodes.begin(), s_nodes.end(), std::greater<cpu_index>());
  vector<cpu_index> v_nodes;
  v_nodes.reserve(num_vertices(instance.virtual1));
  for (Vertex v {0}; v < num_vertices(instance.virtual1); ++v) {
    v_nodes.push_back(make_pair(instance.substrate[v].cpu, v));
  }
  std::sort(v_nodes.begin(), v_nodes.end(), std::greater<cpu_index>());

  // map vertices
  Mapping mapping(instance);
  for (Vertex v {0}; v < num_vertices(instance.virtual1); ++v) {
    if (v_nodes[v].first <= s_nodes[v].first)
      mapping.map_vertex(v_nodes[v].second, s_nodes[v].second);
    else
      throw Infeasible();
  }
  // map edges
  for (Vertex v {0}; v < num_vertices(instance.virtual1); ++v) {
    const Vertex s {mapping.vertex[v]};
    EdgeIterator e, eend;
    for (tie(e, eend) = out_edges(v, instance.virtual1); e != eend; ++e) {
      const Vertex t {mapping.vertex[target(*e, instance.virtual1)]};
      vector<Edge> minpath;
      if (!find_min_path(instance, instance.virtual1[*e].bandwidth,
          s, t, minpath)) {
        mapping.unmap_all_paths(instance);
        cout << "initial solution construction failed" << endl;
        return;
      }
      mapping.map_path(instance, *e, minpath);
    }
  }
  swap(mapping, bundle.best);
  bundle.upper_bound = bundle.best.cost;
}

}

/*
  struct State{
    Vertex v, s;

    bool operator<(State s2) const {
      return v < s2.v;
    }
  };
  queue<State> Q;

  const State initial {0,0};
  Q.push(initial);

  vector<bool> mapped(num_vertices(substrate));

  while(!Q.empty()) {
    State s = Q.front();
    Q.pop();

    if(s.v == num_vertices(virtual1)) {
      //evaluate cost
      //undo last mapping
      s.v-1;
    } else {
      do{
        if(s.s >= num_vertices(substrate)) {
          //unmap v
          break;
        }
      }while(mapped[s.s]);
      //map v on s
    }

    Q.push(s);
  }
  
  return true;

*/
/*
bool backtrack(const Graph& substrate, const Graph& virtual1, Mapping& outMapping) {
  vector<Vertex> m(num_vertices(virtual1));
  vector<bool> mapped(num_vertices(substrate) + 1);
  Vertex v = 0;
  while(true) {
    if(v == num_vertices(virtual1)){
      cout << endl;
      --v;
      mapped[m[v]] = false;
      ++m[v];
      continue;
    } 
    while(mapped[m[v]]) {
      if(m[v] >= num_vertices(substrate)){
        break;
      }
      m[v]++;
    }
    if(m[v] >= num_vertices(substrate)){
      if(v == 0)
        break;
      m[v] = 0;
      --v;
      mapped[m[v]] = false;
      ++m[v];
      continue;
    }
    // map v on m[v]
    mapped[m[v]] = true; 
    ++v;
  }
  return true;
}*/

