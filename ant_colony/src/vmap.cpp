//#define NDEBUG
#include <cassert>
#include <memory>
#include <iostream>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

#include "vmap.h"
#include "util.h"
#include "PriorityQueue.h"

using namespace std;

double best_path(Vertex sSource, Vertex sTarget, int bandwidth,
    Graph& substrate, vector<Edge>& outPred) {
  LOG("VMAP") << "best_path " << sSource << " -> " << sTarget << endl;
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

  while(!Q.empty()) {
    Vertex u = Q.key_top();
    Q.pop();
    if( u == sTarget) {
      return d[sTarget];
    }
    visited[u] =  true;
    EdgeIterator e, eend;
    for(tie(e,eend) = out_edges(u, substrate); e != eend; ++e) {
      Vertex v = target(*e, substrate);
      if(visited[v]) 
        continue;
      int bottleneck = min(bn[u], substrate[*e].bandwidth);
      if(bottleneck < bandwidth)
        continue;
      double newdist = (hops[u] + 1) / bottleneck;
      if(!discovered[v]) {
        d[v] = newdist;
        bn[v] = bottleneck;
        hops[v] = hops[u]+1;
        outPred[v] = *e;
        Q.insert(v);
        discovered[v] = true;
      } else if(d[v] > newdist) {
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
 * Maps a path between sSource and sTarget with a minimum of bandwidth
 * returns the path cost (the number or arcs * bandwidth)
 * retruns -1 if no paths with the required bandwidth were found
 */
Bandwidth edge_map(Vertex sSource, Vertex sTarget, int bandwidth, Graph& substrate, vector<Edge>& outPath) {
  const Edge noedge = *out_edges(sTarget, substrate).second;
  vector<Edge> pred(num_vertices(substrate), noedge);
  if(best_path(sSource, sTarget,bandwidth, substrate, pred) < 0) {
    return -1;
  }
  Bandwidth pathCost {0};
  while(true) {
    pathCost += bandwidth;
    substrate[pred[sTarget]].bandwidth -= bandwidth;
    outPath.push_back(pred[sTarget]);
    sTarget = source(pred[sTarget], substrate);
    if(sTarget == sSource)
      return pathCost;
  }
}

/**
 * Maps virtual vVertex into substrate sVertex, mapping every path
 * between the vVertex already mapped neighbors
 */
bool vertex_map(Vertex vVertex, Vertex sVertex, Graph& substrate,
    Graph& virtual1, Mapping& mapping, Mapping& baseMapping) {
  assert(substrate[sVertex].hasEnoughResources(virtual1[vVertex]));
  mapping.vMap_[vVertex] = sVertex;
  mapping.substrateMapped[sVertex] = true;
  substrate[sVertex].allocate(virtual1[vVertex]);
  // map edges between mapped nodes
  // for every edge e in the neighborhood of vVertex
  EdgeIterator e, eend;
  for(tie(e, eend) = out_edges(vVertex, virtual1); e != eend; ++e) {
    const Vertex vt = target(*e, virtual1),
      st = (virtual1[vt].access) ? baseMapping.vMap_[vt] : mapping.vMap_[vt];
    // if target(e) is already mapped
    // find the minimum cost path between the mapped nodes
    if(st != NONE) {
      Bandwidth pathCost = edge_map(sVertex, st, virtual1[*e].bandwidth, substrate, 
        mapping.edge_map_ref(vVertex, vt));
      if(pathCost < 0)
        return false;
      mapping.cost += pathCost; 
    }
  }
  return true;
}

void getSubstrateCandidates(Vertex u, int hops, Graph& substrate, int cpu,
    int memory, Mapping& mapping, set<Vertex>& outSubstrateCandidates) {
  if(hops > 0) {
    EdgeIterator e, eend;
    for(tie(e, eend) = out_edges(u, substrate); e != eend; ++e) {
      Vertex v = target(*e, substrate);
      getSubstrateCandidates(v, hops-1, substrate, cpu, memory, mapping, 
        outSubstrateCandidates); 
    }
  }
  if(!mapping.substrateMapped[u] && !substrate[u].access 
      && substrate[u].cpu >= cpu && substrate[u].memory >= memory)
    outSubstrateCandidates.insert(u);
}

bool VNEAC::vmap(Graph& substrate, Graph& virtual1, const Parameters& param,
    Mapping& outMapping) {
  // This is the mapping of the virtual access nodes, every solution will
  // share this mapping
  outMapping = Mapping(num_vertices(substrate), num_vertices(virtual1));

  // Build solution components -------------------------------------------
  // * map access nodes
  size_t mapped {0};
  VertexIterator v, vend;

  vector<pair<int, Vertex>> arr;
  for(Vertex v = 0; v < num_vertices(virtual1); ++v)
    if(virtual1[v].access)
      arr.push_back(make_pair(-virtual1[v].cpu, v));

  sort(arr.begin(), arr.end());

  for(auto vv : arr) {
    Vertex v = vv.second;

  //for(tie(v, vend) = vertices(virtual1); v != vend; ++v) {
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
  if(mapped==num_vertices(virtual1)){
    LOG("accept") << "only access nodes " << endl;
    verify(outMapping, substrate, virtual1);
    return true;
  }
  // * order nodes by hanging links
  // * * count number of hanging links
  vector<pair<int, Vertex>> componentCandidates;
  EdgeIterator e, eend;
  for(tie(v, vend) = vertices(virtual1); v != vend; ++v) {
    // the components are those virtual nodes that were not mapped
    if(outMapping.vMap_[*v] != NONE)
      continue;
    componentCandidates.push_back(make_pair(0, *v));
    for(tie(e, eend) = out_edges(*v, virtual1); e != eend; ++e) {
      int alreadyMapped = outMapping.vMap_[target(*e, virtual1)] != NONE;
      componentCandidates.back().first -= alreadyMapped;
    }
  }
  
  // * * sort components by number of hanging links
  sort(componentCandidates.begin(), componentCandidates.end());
  vector<Vertex> components(componentCandidates.size());
  transform(componentCandidates.begin(), componentCandidates.end(),
      components.begin(), [](pair<int, Vertex> cand) { return cand.second; });

  // * initialize a priori information (\eta) and a posteriori information
  // * find substrate candidates for components
  vector<vector<double>> eta(components.size()),
    tau(components.size());
  vector<vector<Vertex>> candidates(components.size());
  const Edge noedge = *out_edges(0, substrate).second;
  vector<Edge> pred(num_vertices(substrate), noedge);
  for(size_t i = 0; i < components.size(); ++i) {
    // find barycenter
    double x = 0, y = 0;
    int mappedNeighbors = 0;
    for(tie(e, eend) = out_edges(components[i], virtual1); e != eend; ++e) {
      Vertex s = outMapping.vMap_[target(*e, virtual1)];
      if(s != NONE) {
        ++mappedNeighbors;
        x += substrate[s].x;
        y += substrate[s].y;
      }
    }
    x /= mappedNeighbors;
    y /= mappedNeighbors;

    tie(v, vend) = vertices(substrate);
    Vertex p = *min_element(v, vend, [&x, &y, &substrate](Vertex u, Vertex v) 
      { return pow(substrate[u].x - x, 2) + pow(substrate[u].y - y, 2) < 
               pow(substrate[v].x - x, 2) + pow(substrate[v].y - y, 2); });
    set<Vertex> substrateCandidates;
    getSubstrateCandidates(p, param.Hops, substrate, virtual1[components[i]].cpu,
        virtual1[components[i]].memory, outMapping, substrateCandidates);

    if(substrateCandidates.empty()) {
      LOG("rejected") << " no component candidates" << endl;
      free_vmap(substrate, virtual1, outMapping);
      return false;
    }
    for(Vertex v: substrateCandidates) {
      // sum the best paths between the candidate and
      // all the mapped neighbors of the compenent
      double bestPathSum = 0;
      EdgeIterator e, eend;
      for(tie(e, eend) = out_edges(components[i], virtual1); e != eend; ++e) {
        const Vertex sk = outMapping.vMap_[target(*e, virtual1)];
    
        if(sk != NONE) {
          const double bp { best_path(v, sk, 
          numeric_limits<int>::max(), substrate, pred)};
          assert(bp < num_vertices(substrate));
          // TODO: there must be a path!
          bestPathSum += (bp > 0) ? bp : num_vertices(substrate);
        }
      }
      candidates[i].push_back(v);
      eta[i].push_back(pow(static_cast<double>(
        substrate[v].cpu + substrate[v].memory + bestPathSum),2));
      tau[i].push_back(param.Tmax);
    }
  }
  // * init 
  //
  // Heuristic Search ----------------------------------------------------
  Mapping bestSolution(numeric_limits<Bandwidth>::max());
  for(int iter = 0; iter < param.IterMax; ++iter) {
    vector<Mapping> solutions(param.NAnts, 
        Mapping(num_vertices(substrate), num_vertices(virtual1)));
    assert(solutions[0].vMap_.size() > 0);
    for(auto& solution : solutions) {
      assert(solution.vMap_.size() > 0);
      vector<double> p(num_vertices(substrate));
      for(size_t c = 0; c < components.size(); c++) {
        const vector<Vertex>& cCandidates {candidates[c]};
        //  * assign probabilities
        double sum = 0;
        for(size_t t = 0; t < cCandidates.size(); ++t) {
          if(!solution.substrateMapped[cCandidates[t]]) {
            sum += p[t] = tau[c][t] * eta[c][t];
          }
        }
        // no candidates for this component
        if(sum == 0) {
          solution.cost = numeric_limits<Bandwidth>::max();
          break;
        }
        //  * select one candidate randomly with probability p_ab
        double res = sum * drand48();
        auto selected = -1;
        do{
          ++selected;
          if(solution.substrateMapped[cCandidates[selected]])
            continue;
          res -= p[selected];
        } while(res >= 0);
        assert(!solution.substrateMapped[cCandidates[selected]]);
        if(!vertex_map(components[c], cCandidates[selected], substrate, 
            virtual1, solution, outMapping)) {
          solution.cost = numeric_limits<Bandwidth>::max();
          break;
        }
      }
      assert(solution.cost > 0);
      free_vmap(substrate, virtual1, solution);
    }
    // select best solution
    auto bestIterationSolution =
      min_element(solutions.begin(), solutions.end(), 
      [](Mapping& sol1, Mapping& sol2){ return sol1.cost < sol2.cost; });
    // store the best overall solution
    bool noFeasibleSolutions = 
      bestIterationSolution->cost == numeric_limits<Bandwidth>::max();
    // TODO: if no feasible solutions, go to another iteration
    if(noFeasibleSolutions) {
      LOG("rejected") << " no feasible solutions" << endl;
      free_vmap(substrate, virtual1, outMapping);
      return false;
    }
    if(bestIterationSolution->cost < bestSolution.cost)
      bestSolution = *bestIterationSolution;
    // update pheromone trail
    for(size_t i = 0; i < components.size(); ++i) {
      for(size_t s = 0; s < candidates.size(); ++s) {
        tau[i][s] *= param.Ro;
      }
      Vertex v = bestIterationSolution->vMap_[components[i]];
      int c = distance(candidates[i].begin(), 
          lower_bound(candidates[i].begin(), candidates[i].end(), v)); 
      assert(candidates[i][c] == v);
      if(param.onlyTheBestReinforces){
        tau[i][c] = truncateBetween(tau[i][c] + 
          (param.Phi / bestIterationSolution->cost), param.Tmin, param.Tmax);
      } else {
        for(auto& solution: solutions) {
          const Vertex s = solution.vMap_[i];
          if(s != NONE)
            tau[i][s] = truncateBetween(tau[i][s]
              + (param.Phi / solution.cost), param.Tmin, param.Tmax);
        }
      }
    }
  }
  // join mapping
  for(size_t v = 0; v != bestSolution.vMap_.size(); ++v) {
    Vertex s = bestSolution.vMap_[v];
    if(s != NONE) {
      outMapping.vMap_[v] = s;
#ifndef NDEBUG
      int cput = substrate[s].cpu;
#endif
      substrate[s].allocate(virtual1[v]);
      assert(cput > substrate[s].cpu);
      assert(substrate[s].cpu >= 0 && substrate[s].memory >= 0);
    }
    //occupy paths 
    for(auto& emap: bestSolution.eMap_[v]) {
      if(emap.second.size() <= 0)
        continue;
      Vertex u = emap.first;
      EdgeIterator e, eend;
      tie(e, eend) = out_edges(v, virtual1);
      Bandwidth bandwidth = virtual1[*find_if(e, eend, 
          [u, &virtual1](Edge e) 
          { return target(e, virtual1) == u; })].bandwidth;

      for(Edge& e: emap.second)
        substrate[e].bandwidth -= bandwidth;
      
      auto& outemap = outMapping.edge_map_ref(u, v);
      //assert(outemap.size() == 0);
      std::swap(emap.second, outemap);
      //assert(outemap.size() > 0);
    }
  }
  outMapping.cost += bestSolution.cost;
  assert(verify(outMapping, substrate, virtual1));
  return true;
}

void VNEAC::free_vmap(Graph& substrate, Graph& virtual1, Mapping& mapping) {
  for(size_t u = 0; u < num_vertices(virtual1); ++u) {
    // free paths
    EdgeIterator e, eend;
    for(tie(e, eend) = out_edges(u, virtual1); e != eend; ++e) {
      if(u > target(*e, virtual1))
        continue;
      vector<Edge>& path = mapping.edge_map_ref(u, target(*e, virtual1));
      for(const Edge& ep: path) {
          substrate[ep].bandwidth += virtual1[*e].bandwidth;
      }
    }
    const Vertex s = mapping.vMap_[u];
    if(s == NONE)
      continue;
    mapping.substrateMapped[s] = false;
    substrate[s].deallocate(virtual1[u]);
  }
}
