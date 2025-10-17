#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <queue>

#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include <ilcplex/ilocplex.h>
#pragma GCC diagnositc pop

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/lookup_edge.hpp>

#include "../include/vnedefs.h"
#include "../include/cg.h"
#include "../include/bprice.h"

using std::cout;
using std::priority_queue;
using std::cerr;
using std::endl;
using std::numeric_limits;
using std::vector;
using std::pair;
using vne::CGModel;
using vne::Path;
using vne::PathVar;
using vne::CplexValue;

namespace{
struct BPBundle {
  IloEnv& env;
  Graph& aux;
  Graph& substrate;
  const Graph& virtual1;
  CGModel& model;
  Timer<std::chrono::milliseconds>& timer;
  vector<Path>& paths;

  CplexValue ub;
  CplexValue lb;
  //TODO store best integer x and z
  Mapping& incumbent;
  
  // debug variables
  int depth;
  int nodes_explored;  

  BPBundle(IloEnv& _env, Graph& _aux, Graph& _substrate, const Graph& _virtual1,
      CGModel& _model, Timer<std::chrono::milliseconds>& _timer, 
      vector<Path>& _paths, Mapping& out_mapping)
      : env(_env), aux(_aux), substrate(_substrate), virtual1(_virtual1),
        model(_model), timer(_timer), 
        paths(_paths), ub(numeric_limits<CplexValue>::max()), lb(0),
        incumbent(out_mapping) {
  incumbent.cost = numeric_limits<Bandwidth>::max();
  depth = 0;
  nodes_explored = 0;
}
};

// generate necessary paths, if a infeasibility is detected, return false
bool genPaths(const Vertex v, const Vertex s, const Graph& substrate,
    const Graph& virtual1, const Graph& aux, CGModel& model, 
    vector<Path>& paths) {
  EdgeIterator e, eend;
  for (tie(e, eend) = out_edges(v, virtual1); e != eend; ++e) {
    const size_t id = virtual1[*e].id;
    cout << "** " << v << " -- " << id << endl;
    if (!model.doesPathExistUsingAuxEdge(id, v, s)) {
      Path path {id, source(*e, virtual1), s,
        target(*e, virtual1), 0, virtual1[*e].bandwidth};

      if (model.getPathThatCoversAuxEdge(aux,
          substrate, model, path)) {
        cout << "adding path";
        print(path, substrate);
        paths.push_back(path);
      } else {
        return false;
      }
    }
  }
  return true;
}

// return true if optimal integer solution is found
bool bp(BPBundle& bundle) {
  if (bundle.timer.reachedTimeLimit())
    return false;
  bundle.nodes_explored++;

  cout << "start solving LP" << endl;
  //solve LP
  CplexValue current_obj_value = vne::cgunsplittable_iteration(bundle.env,
      bundle.aux, bundle.substrate, bundle.virtual1, bundle.timer, bundle.model,
      bundle.paths);

#ifdef EXPORT_MODEL
  std::stringstream filename;
  filename << "bp " << bundle.nodes_explored << ".lp";
  cout << "print " << filename.str() << endl;
  bundle.model.cplex.exportModel(filename.str().c_str());
#endif

  cout << "solved LP " << endl;

  if (current_obj_value >= bundle.ub)
    return false;
  
  bool found_optimal_solution = false;
  //find fractional x[v,s]
  /* this is going to be used if the decision variables with a higher
    value are to be branched first
  struct DecValue {
    Vertex s;
    double value;

    bool operator<(vsValue& other) {
      return value >= other.value;
    }
  };
  priority_queue<DecValue> decValues;*/
  //find fractional 
  for (Vertex v = 0; v < num_vertices(bundle.virtual1); v++) {
    bool fractional {false};
    for (Vertex s = 0; s < num_vertices(bundle.substrate); ++s) {
      const double vsval = bundle.model.cplex.getValue(bundle.model.x[v][s]);
      cout << " " << vsval;
      if (DBL_EQL(1.0, vsval)) {
        break;
      }
      if (!DBL_EQL(0.0, vsval)) {
        fractional = true;
        break;
      }
    }
    cout << endl;
    if (fractional) {
      for (Vertex s = 0; s < num_vertices(bundle.substrate); ++s) {
        bool cannot_host = bundle.model.isNodeFixed(s)
            || bundle.substrate[s].cpu < bundle.virtual1[v].cpu;
        if (cannot_host || bundle.timer.reachedTimeLimit())
          continue;

        cout << "branch on x[" << v << ","<< s << endl;
        // TODO this do not work if v nodes are not mapped in order
        bundle.model.fixNode(v, s);
        if (genPaths(v, s, bundle.substrate, bundle.virtual1, bundle.aux,
            bundle.model, bundle.paths)) {
          bundle.depth++;
          if (bp(bundle))
            return true;
          bundle.depth--;
        }
        bundle.model.freeNode(v, s);
      }

      return false;
    }
  }

  cout << "found all x integer solution" << endl;

  // find two fractional z[p] (p in P^k) with the highest value
  EdgesIterator e, eend;
  for (tie(e, eend) = edges(bundle.virtual1); e != eend; ++e) {
    const size_t k_id = bundle.virtual1[*e].id;
    // strategy: branch on all edges that appear in a path that is partially
    // used in the relaxed solution (z_p > 0.0)
    typedef pair<CplexValue, size_t> ValueIndex;
    priority_queue<ValueIndex, vector<ValueIndex>,
        std::greater<ValueIndex>> pathVarsValuesQ;
    for (int i = 0; i < bundle.model.paths[k_id].size(); ++i) {
      const PathVar& pathVar = bundle.model.paths[k_id][i];
      const double vsval = bundle.model.cplex.getValue(pathVar.var);
      if (DBL_EQL(1.0, vsval))
        break;
      if (!DBL_EQL(0.0, vsval)) {
        pathVarsValuesQ.emplace(vsval, i);
        cout << "insert " << i << " value " << vsval << endl;
      }
    }
    // if there are non integer path variables in P^k
    if (!pathVarsValuesQ.empty()) {
      assert(pathVarsValuesQ.size() > 1);

      if (bundle.timer.reachedTimeLimit())
        return false;
      
      cout << "fractional paths" << endl;
      
      size_t p1 = pathVarsValuesQ.top().second;
      pathVarsValuesQ.pop();
      size_t p2 = pathVarsValuesQ.top().second;
      
      assert(p1 != p2);

      // p1 is the smallest
      if (bundle.model.paths[k_id][p1].path.size() > 
          bundle.model.paths[k_id][p2].path.size())
        std::swap(p1, p2);

      const vector<Edge>& path1 = bundle.model.paths[k_id][p1].path.path_;
      const vector<Edge>& path2 = bundle.model.paths[k_id][p2].path.path_;

      cout << "sizes: " << path1.size() << " and " << path2.size() << endl;
      print(bundle.model.paths[k_id][p1].path, bundle.substrate);
      print(bundle.model.paths[k_id][p2].path, bundle.substrate);
      // branch on first edge that is present in p1 but not in p2
      for (int i = 0; i < path1.size(); i++) {
        size_t id1 = bundle.substrate[path1[i]].id,
               id2 = bundle.substrate[path2[i]].id;
        assert(id2 < num_edges(bundle.substrate));
        assert(id1 < num_edges(bundle.substrate));
        if (id1 != id2) {
          cout << "branch on edge " << id1 << endl;
          bundle.model.edges_hosted[k_id][id1].setUB(0.0);
          bool found_optimal_solution = bp(bundle);
          bundle.model.edges_hosted[k_id][id1].setUB(1.0);
          if (found_optimal_solution)
            return true;
          bundle.model.edges_hosted[k_id][id2].setUB(0.0);
          found_optimal_solution = bp(bundle);
          bundle.model.edges_hosted[k_id][id2].setUB(1.0);

          return found_optimal_solution;
        }
      }
      // should have found a edge by this point?
      assert(false);
    }
  } 
  cout << "total int solution: " << current_obj_value << endl;

  // no fractional decision variables
  if (current_obj_value < bundle.incumbent.cost) {
    bundle.incumbent = bundle.model.getCurrentMapping(bundle.substrate,
        bundle.virtual1);
    bundle.ub = bundle.incumbent.cost;
    if (DBL_EQL(bundle.ub, bundle.lb))
      return true;
  }

  return false;
}
}

bool vne::bprice_rec(Graph& substrate, const Graph& virtual1,
    const Parameters param, Mapping& outMapping) {
  Timer<std::chrono::milliseconds> timer(param.timelimit_in_s * 1000);
  timer.start();
  IloEnv env;
  const CplexValue M {static_cast<CplexValue>(num_edges(virtual1))};
  const CplexValue K {
    static_cast<CplexValue>(num_edges(virtual1) * sum_bandwidth(substrate))};

  Bandwidth lb {0}, ub {numeric_limits<Bandwidth>::max()};
  float cplex_time {0},  // time spent on the LP
         path_time {0};   // time spent on Dijkstra
  vector<Path> paths;
  Graph aux;
  buildAuxiliaryGraph(substrate, virtual1, aux);
  
  if (!vne::get_initial_paths_from_aux(substrate, virtual1, aux, paths)) {
    cout << "infeasible problem" << endl;
    return false;
  }
  cout << "@initialColumns:"<< paths.size() << endl;

  CGModel model(env, M, K);
  model.init(env, substrate, virtual1);
  // Presolve
  model.cplex.setParam(IloCplex::PreInd, false);
  model.cplex.setParam(IloCplex::PreDual, -1);
  //model.cplex.setParam(IloCplex::SimDisplay, 2);
  //model.cplex.setParam(IloCplex::AggInd, false);
  // Perturbation switch
  //model.cplex.setParam(IloCplex::PerInd, true);
  model.cplex.setParam(IloCplex::EpOpt, 1e-7);
  //model.cplex.setParam(IloCplex::PreLinear, 0);
  //model.cplex.setParam(IloCplex::RootAlg, CPX_ALG_PRIMAL);

  BPBundle bundle(env, aux, substrate, virtual1, model, timer, paths, outMapping);
  bundle.lb = sum_bandwidth(virtual1);

  bp(bundle);

  bool maxCost =
      outMapping.cost >= numeric_limits<decltype(outMapping.cost)>::max();
  cout << "@bpTimeLimit:" << timer.reachedTimeLimit() << endl;
  cout << "@bpFeasible:"  << !maxCost << endl;
  cout << "@bpInfeasible:"
      << (!timer.reachedTimeLimit() && maxCost) << endl;
  cout << "@bpOptimal:" << (!timer.reachedTimeLimit() && !maxCost) << endl;
  cout << "@bpCost: " << outMapping.cost << endl;
  cout << "@bpnodes_explored: " << bundle.nodes_explored << endl;
  cout << "@bpTime:" << timer.total() / 1000.0 << endl;
  env.end();
  if (outMapping.cost < numeric_limits<Bandwidth>::max())
    return true;
  else
    return false;
}
