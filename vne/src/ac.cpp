// Copyright 2014 Leonardo Moura
#include <cstdio>
#include <sys/times.h>
#include <cassert>
#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <limits>
#include <algorithm>
#include <numeric>
#include <exception>
#include <functional>
#include "../include/ac.h"
#include "../include/vnedefs.h"
#include "PriorityQueue.h"

using std::endl;
using std::min;
using std::max;
using std::cout;
using std::pair;
using std::make_pair;
using std::queue;
using std::swap;
using std::tie;
using std::vector;
using std::set;
using std::numeric_limits;

double truncateBetween(double x, double xmin, double xmax) {
  return max(min(x, xmax), xmin);
}

double best_path(Vertex sSource, Vertex sTarget, int bandwidth,
    Graph& substrate, vector<Edge>& outPred) {
  double* d;
  PriorityQueue Q(num_vertices(substrate), &d);
  vector<int> hops(num_vertices(substrate)),
    bn(num_vertices(substrate));
  vector<bool> visited(num_vertices(substrate)+1);
  vector<bool> discovered(num_vertices(substrate)+1);

  d[sTarget] = numeric_limits<int>::max();
  d[sSource] = 0;
  bn[sSource] = numeric_limits<int>::max();
  hops[sSource] = 0;
  discovered[sSource] = true;
  Q.insert(sSource);

  while (!Q.empty()) {
    Vertex u = Q.key_top();
    Q.pop();
    if ( u == sTarget) {
      return d[sTarget];
    }
    visited[u] =  true;
    EdgeIterator e, eend;
    for (tie(e,eend) = out_edges(u, substrate); e != eend; ++e) {
      Vertex v = target(*e, substrate);
      if (visited[v]) 
        continue;
      int bottleneck = min(bn[u], substrate[*e].bandwidth);
      if (bottleneck < bandwidth)
        continue;
      double newdist = (hops[u] + 1) / bottleneck;
      if (!discovered[v]) {
        d[v] = newdist;
        bn[v] = bottleneck;
        hops[v] = hops[u]+1;
        outPred[v] = *e;
        Q.insert(v);
        discovered[v] = true;
      } else if (d[v] > newdist) {
        outPred[v] = *e;
        bn[v] = bottleneck;
        hops[v] = hops[u]+1;
        d[v] = newdist;
        Q.decrease(v);
      }
    }
  }
  return -1;
}

/**
 * Map a path between sSource and sTarget with a minimum of bandwidth
 * return the path cost (the number or arcs * bandwidth)
 * return -1 if no paths with the required bandwidth was found
 */
Bandwidth edge_map(Vertex sSource, Vertex sTarget, int bandwidth, Graph& substrate, vector<Edge>* outPath) {
  // fill the vector pred with a marker indicating no edge
  /*vector<Edge> pred(num_vertices(substrate), numeric_limits<Edge>::max());
  if (best_path(sSource, sTarget,bandwidth, substrate, pred) < 0) {
    return -1;
  }
  Bandwidth pathCost {0};
  // compute the best path from the vector pred
  do {
    pathCost += bandwidth;
    substrate[pred[sTarget]].bandwidth -= bandwidth;
    assert(substrate[pred[sTarget]].bandwidth >= 0);
    outPath->push_back(pred[sTarget]);
    sTarget = source(pred[sTarget], substrate);
  } while (sTarget != sSource);
  return pathCost;
  */
  return 0;
}

/**
 * Maps virtual vVertex into substrate sVertex, mapping every path
 * between the vVertex already mapped neighbors
 */
bool vertex_map(Vertex vVertex, Vertex sVertex, Graph& substrate,
    const Graph& virtual1, Mapping& mapping) {
  assert(substrate[sVertex].hasEnoughResources(virtual1[vVertex]));
  mapping.vertex[vVertex] = sVertex;
  mapping.substrate_mapped[sVertex] = true;
  substrate[sVertex].allocate(virtual1[vVertex]);
  // map edges between mapped nodes
  // for every edge e in the neighborhood of vVertex
  EdgeIterator e, eend;
  for (tie(e, eend) = out_edges(vVertex, virtual1); e != eend; ++e) {
    const Vertex vt = target(*e, virtual1);
    
    // access was probably already mapped?
    //const Vertex st = (virtual1[vt].access) ? baseMapping.vMap_[vt] : mapping.vMap_[vt];
    const Vertex st = mapping.vertex[vt];
    // if target(e) is already mapped
    // find the minimum cost path between the mapped nodes
    if (st != NONE) {
      cout << "map edge between " << vVertex << " mapped to " << sVertex << " and " << vt << " mapped to " << st << endl;
      Bandwidth pathCost = edge_map(sVertex, st, virtual1[*e].bandwidth, substrate, 
        mapping.edge_map_ref(vVertex, vt));
      if (pathCost < 0)
        return false;
      assert(mapping.edge_map_ref(vVertex, vt)->size() > 0);
      mapping.cost += pathCost; 
    }
  }
  return true;
}

void getSubstrateCandidates(Vertex u, int hops, Graph& substrate, int cpu,
    Mapping& mapping, set<Vertex>& outSubstrateCandidates) {
  if(hops > 0) {
    EdgeIterator e, eend;
    for(tie(e, eend) = out_edges(u, substrate); e != eend; ++e) {
      Vertex v = target(*e, substrate);
      getSubstrateCandidates(v, hops-1, substrate, cpu, mapping, 
        outSubstrateCandidates); 
    }
  }
  if(!mapping.substrate_mapped[u] 
      && substrate[u].cpu >= cpu)
    outSubstrateCandidates.insert(u);
}

void free_vmap(Graph& substrate, const Graph& virtual1, Mapping& mapping) {
  for(size_t u = 0; u < num_vertices(virtual1); ++u) {
    // free paths
    EdgeIterator e, eend;
    for(tie(e, eend) = out_edges(u, virtual1); e != eend; ++e) {
      if (u < target(*e, virtual1))
        continue;
      vector<Edge>* path = mapping.edge_map_ref(u, target(*e, virtual1));
      for(const Edge& ep: *path) {
          assert(substrate[ep].bandwidth >= 0);
          substrate[ep].bandwidth += virtual1[*e].bandwidth;
          assert(substrate[ep].bandwidth >= 0);
          cout << "free " << ep << " b = " << virtual1[ep].bandwidth << endl;
      }
    }
    const Vertex s = mapping.vertex[u];
    if(s == NONE)
      continue;
    mapping.substrate_mapped[s] = false;
    substrate[s].deallocate(virtual1[u]);
  }
}


bool vne::ac(Graph& substrate, const Graph& virtual1,
    const Parameters param, Mapping& outMapping) {
  cout << "starting" << endl;
  // This is the mapping of the virtual access nodes, every solution will
  // share this mapping
  outMapping = Mapping(substrate, virtual1);

  // Build solution components -------------------------------------------
  // * map access nodes
  /* this used to map access nodes
  size_t mapped {0};
  VertexIterator v, vend;

  vector<pair<int, Vertex>> arr;
  
  for(Vertex v = 0; v < num_vertices(virtual1); ++v)
    if(virtual1[v].access)
      arr.push_back(make_pair(-virtual1[v].cpu, v));

  sort(arr.begin(), arr.end());

  for(auto vv : arr) {
    Vertex v = vv.second;

    cout << v << endl;
    if(virtual1[v].access) {
      Vertex smax {NONE};
      VertexIterator s, send;
      for(tie(s, send) = vertices(substrate); s != send; ++s) {
        if(substrate[*s].access && !outMapping.substrateMapped[*s]){
        }
        if(outMapping.substrateMapped[*s] || !substrate[*s].access
            || !substrate[*s].hasEnoughResources(virtual1[v]))
          continue;
        if((smax == NONE) || (substrate[smax] < substrate[*s]))
          smax = *s;
      }
      // no possible substrate nodes
      // TODO: if vertex cannot be mapped to this substrate vertex, try others
      if(smax == NONE 
          || !vertex_map(v, smax, substrate, virtual1, outMapping, outMapping)) {
        LOG("rejected") << " no substrate access nodes " << endl;
        free_vmap(substrate, virtual1, outMapping);
        return false;
      }
      assert(outMapping.vMap_[v] == smax);
      ++mapped;
    }
  }
  */
  // * order nodes by hanging links
  // * * count number of hanging links
  /*
  cout << "counting hanging links" << endl;
  vector<pair<int, Vertex>> componentCandidates;
  VertexIterator v, vend;
  EdgeIterator e, eend;
  for (tie(v, vend) = vertices(virtual1); v != vend; ++v) {
    // the components are those virtual nodes that were not mapped
    if (outMapping.vertex[*v] != NONE)
      continue;
    componentCandidates.push_back(make_pair(0, *v));
    for (tie(e, eend) = out_edges(*v, virtual1); e != eend; ++e) {
      int alreadyMapped = outMapping.vertex[target(*e, virtual1)] != NONE;
      componentCandidates.back().first -= alreadyMapped;
    }
  }*/
  
  // * * sort components by number of hanging links
  //sort(componentCandidates.begin(), componentCandidates.end());
  vector<Vertex> components(num_vertices(virtual1));
  std::iota(components.begin(), components.end(), 0);
  //transform(componentCandidates.begin(), componentCandidates.end(),
  //    components.begin(), [](pair<int, Vertex> cand) { return cand.second; });
  sort(components.begin(), components.end(), [&](Vertex u, Vertex v) { return out_degree(u, virtual1) < out_degree(v, virtual1); });

  // * initialize a priori information (\eta) and a posteriori information
  // * find substrate candidates for components
  cout << "initializing heuristic information" << endl;
  vector<vector<double>> eta(components.size()),
    tau(components.size());
  vector<vector<Vertex>> candidates(components.size());
  const Edge noedge = *out_edges(0, substrate).second;
  vector<Edge> pred(num_vertices(substrate), noedge);
  for(size_t i = 0; i < components.size(); ++i) {
    const Vertex v {components[i]}; 

    for(Vertex s = 0; s < num_vertices(substrate); ++s) {
      if (substrate[s].cpu < virtual1[v].cpu)
        continue;

      // sum the best paths between the candidate and
      // all the mapped neighbors of the compenent
      double bestPathSum = 0;
      EdgeIterator e, eend;
      for(tie(e, eend) = out_edges(components[i], virtual1); e != eend; ++e) {
        const Vertex sk = outMapping.vertex[target(*e, virtual1)];
    
        if(sk != NONE) {
          const double bp { best_path(s, sk, 
          numeric_limits<int>::max(), substrate, pred)};
          assert(bp < num_vertices(substrate));
          // TODO: there must be a path!
          bestPathSum += (bp > 0) ? bp : num_vertices(substrate);
        }
      }
      candidates[i].push_back(s);
      eta[i].push_back(pow(static_cast<double>(
        substrate[s].cpu + bestPathSum),2));
      tau[i].push_back(param.tmax);
    }

    if (candidates[i].size() == 0)
      return false;
  }

  for (int i = 0; i < components.size(); i++) {
    const Vertex v {components[i]};
    cout << v << " has " << candidates[i].size() << " candidates " << endl;

    for (int j = 0; j < candidates[i].size(); j++)
      cout << "   * " << candidates[i][j] << " t = " << tau[i][j] << ", eta = " << eta[i][j] << endl; 
  }

  // * init 
  //
  // Heuristic Search ----------------------------------------------------
  Mapping bestSolution(numeric_limits<Bandwidth>::max());
  for (int iter = 0; iter < param.maxiter; ++iter) {
    cout << "::::iteration " << iter << endl;
    vector<Mapping> solutions(param.ants, 
        Mapping(substrate, virtual1));
    for (auto& solution : solutions) {
      assert(solution.vertex.size() > 0);
      vector<double> p(num_vertices(substrate));
      for (size_t c = 0; c < components.size(); c++) {
        const vector<Vertex>& cCandidates {candidates[c]};
        //  * assign probabilities
        double sum = 0;
        for (size_t t = 0; t < cCandidates.size(); ++t) {
          if (!solution.substrate_mapped[cCandidates[t]]) {
            sum += p[t] = tau[c][t] * eta[c][t];
          }
        }
        // no candidates for this component
        if (sum == 0) {
          solution.cost = numeric_limits<Bandwidth>::max();
          break;
        }
        //  * select one candidate randomly with probability p_ab
        double res = sum * drand48();
        auto selected = -1;
        do{
          ++selected;
          if (solution.substrate_mapped[cCandidates[selected]])
            continue;
          res -= p[selected];
        } while(res >= 0);
        assert(!solution.substrate_mapped[cCandidates[selected]]);
        if (!vertex_map(components[c], cCandidates[selected], substrate, 
            virtual1, solution)) {
          cout << "failed mapping" << endl;
          solution.cost = numeric_limits<Bandwidth>::max();
          break;
        }
      }
      assert(solution.cost > 0);
      cout << "free " << endl;
      free_vmap(substrate, virtual1, solution);
      assert(verify(solution, substrate, virtual1));
    }
    // select best solution
    auto bestIterationSolution =
      min_element(solutions.begin(), solutions.end(), 
      [](Mapping& sol1, Mapping& sol2){ return sol1.cost < sol2.cost; });
    // store the best overall solution
    bool noFeasibleSolutions = 
      bestIterationSolution->cost == numeric_limits<Bandwidth>::max();
    if (noFeasibleSolutions) {
      cout << "no feasible solution" << endl;
      continue;
    }
    cout << "feasible solution found" << endl;
    if (bestIterationSolution->cost < bestSolution.cost)
      bestSolution = *bestIterationSolution;
    assert(verify(bestSolution, substrate, virtual1));
    // update pheromone trail
    for (size_t i = 0; i < components.size(); ++i) {
      for (size_t s = 0; s < candidates.size(); ++s) {
        tau[i][s] *= param.ro;
      }
      Vertex v = bestIterationSolution->vertex[components[i]];
      int c = distance(candidates[i].begin(), 
          lower_bound(candidates[i].begin(), candidates[i].end(), v)); 
      assert(candidates[i][c] == v);
      if (param.onlythebestreinforces){
        tau[i][c] = truncateBetween(tau[i][c] + 
          (param.phi / bestIterationSolution->cost), param.tmin, param.tmax);
      } else {
        for (auto& solution: solutions) {
          const Vertex s = solution.vertex[i];
          if (s != NONE)
            tau[i][s] = truncateBetween(tau[i][s]
              + (param.phi / solution.cost), param.tmin, param.tmax);
        }
      }
    }
  }
  // apply mapping
  for (size_t v = 0; v != bestSolution.vertex.size(); ++v) {
    Vertex s = bestSolution.vertex[v];
    if (s != NONE) {
      outMapping.vertex[v] = s;
#ifndef NDEBUG
      int cput = substrate[s].cpu;
#endif
      substrate[s].allocate(virtual1[v]);
      assert(cput > substrate[s].cpu);
      assert(substrate[s].cpu >= 0);
    }
    //occupy paths 
    for (auto emap_v : bestSolution.edge_map_ref(v)) {
      if (emap_v == nullptr)
        continue;
      for (Vertex u = 0; u < emap_v->size(); u++) {

        EdgeIterator e, eend;
        tie(e, eend) = out_edges(v, virtual1);
        Bandwidth bandwidth = virtual1[*find_if(e, eend, 
            [u, &virtual1](Edge e) 
            { return target(e, virtual1) == u; })].bandwidth;

        for(Edge& e: emap_v[u])
          substrate[e].bandwidth -= bandwidth;
        
        //assert(outemap.size() == 0);
        
        // TODO this assignment
        //std::swap(emap_v[u], outMapping.edge_map_ref(u, v));
        //assert(outemap.size() > 0);
      }
    }
  }
  outMapping.cost += bestSolution.cost;
  assert(verify(outMapping, substrate, virtual1));
  print(outMapping);
  return true;
}
