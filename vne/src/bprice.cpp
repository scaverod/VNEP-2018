#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <queue>
#include <numeric>

#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include <ilcplex/ilocplex.h>
#pragma GCC diagnositc pop

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/lookup_edge.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>

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

using vne::BPNode;

extern float primal_time, pricing_time, addc_time, aux_time; 

BPNode::BPNode(size_t _depth, int _father,
    size_t _c1, size_t _c2, CplexValue _dual, float _infeasibilities) :  depth(_depth), father(_father),
    c1(_c1), c2(_c2), dual(_dual), infeasibilities(_infeasibilities) { }

namespace {

/*
struct bfs_node {
};

struct dfs_node {
  const float dual;
};

struct projection_node {
  const float infeasibility;
};
*/

struct BPRootNode : public BPNode {
  BPRootNode() : BPNode(0,0,0,0,0,0) { }
  void apply(CGModel& model) const { };
  void undo(CGModel& model) const { };
};

struct BPForbidNode : public BPNode {
  BPForbidNode(size_t _depth, int _father,
      size_t _c1, size_t _c2, CplexValue _dual, float _infeasibilities) : BPNode(_depth,
      _father,
      _c1, _c2, _dual, _infeasibilities) { }

  void apply(CGModel& model) const {
    model.x[c1][c2].setUB(0);
  }
  
  void undo(CGModel& model) const {
    model.x[c1][c2].setUB(1);
  }
};

struct BPFixNode : public BPNode {
  BPFixNode(size_t _depth, int _father,
      size_t _c1, size_t _c2, CplexValue _dual, float _infeasibilities) : BPNode(_depth,
      _father,
      _c1, _c2, _dual, _infeasibilities) { }

  void apply(CGModel& model) const {
    model.fixNode(c1, c2);
  }
  
  void undo(CGModel& model) const {
    model.freeNode(c1, c2);
  }
};

// set x[c1, [c2, end_s]] to 0
struct BPDichotomicForbidNode : public BPNode {
  size_t end_s;

  BPDichotomicForbidNode(size_t _depth, int _father,
      size_t _c1, size_t _c2, size_t _end_s, CplexValue _dual,
      float _infeasibilities) : BPNode(_depth,
      _father, _c1, _c2, _dual, _infeasibilities), end_s(_end_s) { }

  void apply(CGModel& model) const {
    for (size_t i = c2; i <= end_s; i++)
      model.x[c1][i].setUB(0);
  }
  
  void undo(CGModel& model) const {
    for (size_t i = c2; i <= end_s; i++)
      model.x[c1][i].setUB(1);
  }
};

struct BPEdgeNode : public BPNode {
  BPEdgeNode(size_t _depth, int _father,
      size_t _c1, size_t _c2, CplexValue _dual, float _infeasibilities) : BPNode(_depth,
      _father, _c1, _c2, _dual, _infeasibilities) { }

  void apply(CGModel& model) const {
    model.edges_hosted[c1][c2].setUB(0);
  }

  void undo(CGModel& model) const {
    model.edges_hosted[c1][c2].setUB(1);
  }
};

void print_node(const BPNode*  node) {
  if (const BPNode* n = dynamic_cast<const BPRootNode*>(node)) {
    cout << "*** root" << endl;
  } else if (const BPNode* n = dynamic_cast<const BPFixNode*>(node)) {
    cout << "*** fix node " << node->c1 << ", " << node->c2 << endl; 
  } else if (const BPNode* n = dynamic_cast<const BPForbidNode*>(node)) {
    cout << "*** forbid node " << node->c1 << ", " << node->c2 << endl; 
  } else {
    cout << "*** edge node?" << endl;
  }
}

bool heuristic(const CGModel& model, const Graph& substrate,
    const Graph& virtual1, Mapping& out_mapping) {
  out_mapping = Mapping(substrate, virtual1);
  Graph temp_substrate(substrate);
  // map nodes according to current solution
  for (Vertex v = 0; v < num_vertices(virtual1); v++) {
    CplexValue max_value = 0;
    Vertex best_s = numeric_limits<Vertex>::max();
    for (Vertex s = 0; s < num_vertices(substrate); ++s) {
      const CplexValue val = model.cplex.getValue(model.x[v][s]);
      if (val > max_value && !out_mapping.substrate_mapped[s]) {
        max_value = val;
        best_s = s;
        if (val == 1.0)
          break;
      }
    }
    if (best_s == numeric_limits<Vertex>::max())
      return false;
    out_mapping.vertex[v] = best_s;
    out_mapping.substrate_mapped[best_s] = true;
  }

  // map edges
  EdgesIterator e, eend;
  for (tie(e, eend) = edges(virtual1); e != eend; ++e) {
    const Vertex v = source(*e, virtual1);
    const Vertex w = target(*e, virtual1);

    vector<Edge> min_path;
    // cout << "search " << v << " to " << w << endl;
    if (!vne::find_min_path(temp_substrate, virtual1[*e].bandwidth,
          out_mapping.vertex[v], out_mapping.vertex[w], min_path))
      return false;
  
    out_mapping.map_path_using_bandwidth(temp_substrate, virtual1, *e,
        std::move(min_path));
  }
  return true;
}

// min hop between any two substrate nodes
// only half of the matrix is filled
void precalc(const Graph& substrate, vector<vector<Bandwidth>>& min_hop) {
  min_hop = vector<vector<Bandwidth>>(num_vertices(substrate));
  vector<bool> visited(num_vertices(substrate));
  for (int s = 1; s < num_vertices(substrate); ++s) {
    min_hop[s] = vector<Bandwidth>(s); 
    vector<int> distance(num_vertices(substrate));
    visited.clear();
    std::queue<Vertex> Q;
    Q.push(s);
    while (!Q.empty()) {
      Vertex curr = Q.front();
      Q.pop();
      EdgeIterator e, eend;
      for (tie(e, eend) = out_edges(curr, substrate); e != eend; ++e) {
        Vertex t = target(*e, substrate);
        if (!visited[t]) {
          visited[t] = true;
          if(s > t)
            min_hop[s][t] = distance[t] = distance[curr] + 1;
          Q.push(t);
        }
      }
    }
  }
}

bool select_variable_first_best(const size_t sub_size,
    const vector<Vertex>& v_nodes, const CGModel& model,
    Vertex& sel_v, Vertex& sel_s) {
  //CplexValue max_infeasibility = 0;
  CplexValue min_infeasibility = .5;
  bool fractional = false;
  for (Vertex v : v_nodes) { 
    for (Vertex s = 0; s < sub_size; ++s) {
      const double vsval = model.cplex.getValue(model.x[v][s]);
      if (DBL_EQL(1.0, vsval)) {
        break;
      }
      if (DBL_EQL(.0, vsval)) {
        continue;
      }
      const double infeasibility =
          std::min(vsval, 1.0 - vsval);
      //if (infeasibility > max_infeasibility) {
      if (infeasibility < min_infeasibility) {
        //max_infeasibility = infeasibility;
        min_infeasibility = infeasibility;
        sel_v = v;
        sel_s = s;
        fractional = true;
      }
    }
    if (fractional) {
      return true;
    }
  }
  return false;
}
}

void vne::branch_on_nodes(const CGModel&, const Vertex sel_v, const Vertex sel_s,
    const size_t current_i, BPTree& tree) {
  const BPNode* current_node = tree.nodes[current_i];
  tree.push(new BPFixNode(current_node->depth + 1,
      current_i, sel_v, sel_s, current_node->dual, current_node->infeasibilities));
  tree.push(new BPForbidNode(current_node->depth + 1,
      current_i, sel_v, sel_s, current_node->dual, current_node->infeasibilities));
}

void vne::BPTree::push(BPNode* node) {
  nodes.emplace_back(node);
  Q.push(nodes.size()-1);
}
size_t vne::BPTree::pop() {
  const size_t current_i = Q.top();
  Q.pop();
  return current_i;
}

void vne::branch_on_nodes_dichotomic(const CGModel& model, const Vertex sel_v,
    const Vertex sel_s,
    const size_t current_i, BPTree& tree) { 
  // find interval begin
  int begin = sel_s;
  for (;begin > 0 && model.x[sel_v][begin].getUB() > 0; begin--) {}

  // find interval end
  int end = sel_s;
  for (;end < model.x[sel_v].getSize()-1 && model.x[sel_v][end].getUB() > 0; end++) {}

  //cout << begin << " -- " << end << endl;
  /*if (end - begin == 1) {
    cout << "base case" << endl;
    nodes.emplace_back(new BPFixNode(nodes[current_i]->depth + 1,
        current_i, sel_v, begin, nodes[current_i]->dual));
    queue.push(nodes.size()-1);
    nodes.emplace_back(new BPFixNode(nodes[current_i]->depth + 1,
        current_i, sel_v, end, nodes[current_i]->dual));
    queue.push(nodes.size()-1);
  } else {*/
    //cout << " * "  << (end-begin)/2 << endl;
    const BPNode* current_node = tree.nodes[current_i];
    tree.push(new BPDichotomicForbidNode(current_node->depth + 1,
        current_i, sel_v, begin, (end-begin)/2, current_node->dual,
        current_node->infeasibilities));

    tree.push(new BPDichotomicForbidNode(tree.nodes[current_i]->depth + 1,
        current_i, sel_v, ((end-begin)/2)+1, end, current_node->dual,
        current_node->infeasibilities));
  //}
}

template <vne::BPStrategy Strategy,
    vne::BranchOnNodesFun branch_on_nodes_fun = &vne::branch_on_nodes>
bool vne::bprice(Graph& substrate, const Graph& virtual1,
    const Parameters param, Mapping& out_mapping) {
  Timer<std::chrono::milliseconds> timer(param.timelimit_in_s * 1000);
  Timer<std::chrono::milliseconds> profiler;
  timer.start();
  IloEnv env;

  const Bandwidth available_band = sum_bandwidth(substrate);
  
  const CplexValue M {static_cast<CplexValue>(num_edges(virtual1))};
  const CplexValue K {
    static_cast<CplexValue>(available_band + 1)};

  CplexValue lb {static_cast<CplexValue>(sum_bandwidth(virtual1))},
             ub {static_cast<CplexValue>(available_band)};
  float cplex_time {0},  // time spent on the LP
         path_time {0};   // time spent on Dijkstra
  float select_time {0}; // time selecting variables
  float firstint_time {-1}; // time to obtain the first int sol
  vector<Path> paths;
  Graph aux;
  vne::buildAuxiliaryGraphRestrict(substrate, virtual1, aux);
  
  /*if (!vne::get_initial_paths_from_aux(substrate, virtual1, aux, paths)) {
    cout << "infeasible problem" << endl;
  } else {
    cout << "@initialColumns:"<< paths.size() << endl;
  */
    out_mapping.cost = numeric_limits<decltype(out_mapping.cost)>::max();

    CGModel model(env, M, K);
    model.init(env, substrate, virtual1);
    // Presolve
    model.cplex.setParam(IloCplex::IloCplex::Threads, 1);
    model.cplex.setParam(IloCplex::PreInd, false);
    model.cplex.setParam(IloCplex::PreDual, -1);
    model.cplex.setParam(IloCplex::SimDisplay, 0);
    // the invoking object will do basis pivots in order to maintain a 
    // valid basis when variables or constraints are removed
    model.cplex.setDeleteMode(IloCplex::FixBasis);
    //model.cplex.setParam(IloCplex::AggInd, false);
    // Perturbation switch
    //model.cplex.setParam(IloCplex::PerInd, true);
    model.cplex.setParam(IloCplex::EpOpt, 1e-7);
    //model.cplex.setParam(IloCplex::PreLinear, 0);
    //model.cplex.setParam(IloCplex::RootAlg, CPX_ALG_PRIMAL);
    
    //vector<vector<Bandwidth>> min_hop;
    //precalc(substrate, min_hop);

    vector<Vertex> v_nodes(num_vertices(virtual1));
    std::iota(v_nodes.begin(), v_nodes.end(), 0);
    std::sort(v_nodes.begin(), v_nodes.end(),
        [&](Vertex u, Vertex v) { return degree(u, virtual1) > degree(v, virtual1); });

    //vector<BPNode*> nodes;
    //vne::QueueType Q;
    BPTree tree;

    if (Strategy == vne::DFS) {
      tree.Q = priority_queue<size_t, vector<size_t>,
              std::function<bool(size_t, size_t)>>(
              std::greater<size_t>());
    } else if (Strategy == vne::BFS) {
      tree.Q = priority_queue<size_t, vector<size_t>,
              std::function<bool(size_t, size_t)>>(
            [&](size_t i, size_t j)
            { return tree.nodes[i]->dual > tree.nodes[j]->dual; } );
    } else {
      tree.Q = priority_queue<size_t, vector<size_t>,
              std::function<bool(size_t, size_t)>>(
            [&](size_t i, size_t j)
            { return (tree.nodes[i]->dual +  (ub - tree.nodes[0]->dual)
                  / tree.nodes[0]->infeasibilities
                  * tree.nodes[i]->infeasibilities)
                > (tree.nodes[j]->dual + (ub - tree.nodes[0]->dual)
                  / tree.nodes[0]->infeasibilities
                  * tree.nodes[j]->infeasibilities); });
    }
    
    // add root to the queue
    tree.push(new BPRootNode());
    
    int nodes_explored = 0, // including root
        //vmap_explored = 0,  // node-type cuts explored
        integer_solutions_found = 0;

    size_t previous_i = 0;
    while (!tree.Q.empty()) {
      if (timer.reachedTimeLimit())
        break;
      const size_t current_i = tree.pop();
      nodes_explored++;

#ifdef VERBOSE
      cout << "Q size = " << tree.Q.size() << " time: " << timer.total()/1000.0 << endl;
      cout << "(" << lb << ", "<< ub << ")" << endl;
      cout << "dual:" << tree.nodes[current_i]->dual << endl; 
#endif

      //vmap_explored += (nodes[current_i].type == NODE);

      // find common ancestor
      size_t a_prev = previous_i, a_curr = current_i;
      // put them at the same depth
      while (tree.nodes[a_prev]->depth > tree.nodes[a_curr]->depth) {
        //print_node(tree.nodes[a_prev]);
        a_prev = tree.nodes[a_prev]->father;
      }
      while (tree.nodes[a_curr]->depth > tree.nodes[a_prev]->depth) {
        //print_node(tree.nodes[a_curr]);
        a_curr = tree.nodes[a_curr]->father;
      }
      while (a_curr != a_prev) {
        //print_node(tree.nodes[a_prev]);
        //print_node(tree.nodes[a_curr]);
        a_curr = tree.nodes[a_curr]->father;
        a_prev = tree.nodes[a_prev]->father;
      }
      
      //cout << "undo starting with " << a_prev << endl;
      for (int i = previous_i; i != a_prev; i = tree.nodes[i]->father) {
        //cout << "undo " << i << " "; print_node(tree.nodes[i]);
        tree.nodes[i]->undo(model);
      }

      /*for (Vertex v = 0; v < num_vertices(virtual1); v++) {
        for (Vertex s = 0; s < num_vertices(substrate); s++) {
          assert(model.x[v][s].getLB() == 0);
          assert(model.x[v][s].getUB() == 1);
        }
      }*/

      previous_i = current_i;

      //cout << "apply starting with " << current_i << endl;
      for (int i = current_i; i != a_curr; i = tree.nodes[i]->father) {
        //cout << "apply " << i << " "; print_node(tree.nodes[i]);
        tree.nodes[i]->apply(model);
      }
      profiler.start();     
      BPNode& current = *tree.nodes[current_i];
      current.dual = vne::cgunsplittable_iteration(env, aux,
          substrate, virtual1, timer, model, paths);
      cplex_time += profiler.total();
#ifdef VERBOSE
      cout << "current node:" << current_i << " "; // TODO print(nodes[current_i]);
      cout << "lb = " << current.dual << endl;
#endif

      // if current is infeasible, this value is Infinity
      if (current.dual >= ub)
        continue;

      Mapping heuristic_mapping;
      if (heuristic(model, substrate, virtual1, heuristic_mapping)) {
        if (heuristic_mapping.cost < ub) {
          if (firstint_time < 0)
            firstint_time = timer.total();
          integer_solutions_found++;
          out_mapping = heuristic_mapping;
          ub = out_mapping.cost;
#ifdef VERBOSE
          cout << "ub:"<< ub << " > " << (timer.total() / 1000.0) <<endl;
#endif
          if (ub == lb)
            break;
        }
      }

      if (Strategy == BEST_PROJECTION) {
        current.infeasibilities = 0;
        for (Vertex v = 0; v < num_vertices(virtual1); v++) {
          for (Vertex s = 0; s < num_vertices(substrate); s++) {
            const CplexValue val = model.cplex.getValue(model.x[v][s]);
            if (val == 1.0)
              break;
            current.infeasibilities += std::min(val - floor(val),
                ceil(val) - val);
          }
        }
      }

      profiler.start();
      // branch on nodes
      Vertex sel_v, sel_s;
      bool fractional = select_variable_first_best(num_vertices(substrate),
          v_nodes, model, sel_v, sel_s);

      if (fractional) {
        
        (*branch_on_nodes_fun)(model, sel_v, sel_s, current_i, tree);
#ifdef VERBOSE
        cout << "branch on x[" << sel_v << "," << sel_s << "]" << endl;
        cout << "insert " << tree.nodes.size() << " on the queue" << endl;
        // TODO print(nodes[nodes.size()-2]);
        // print(nodes[nodes.size()-1]);
#endif
        continue;
      }
      // x is integral, search for fractional z
      // find two fractional z[p] (p in P^k) with the highest value
      EdgesIterator e, eend;
      for (tie(e, eend) = edges(virtual1); e != eend; ++e) {
        const size_t k_id = virtual1[*e].id;
        // strategy: branch on all edges that appear in a path that is partially
        // used in the relaxed solution (0 < z_p < 1)
        typedef pair<CplexValue, size_t> ValueIndex;
        priority_queue<ValueIndex, vector<ValueIndex>,
            std::greater<ValueIndex>> pathVarsValuesQ;
        for (int i = 0; i < model.paths[k_id].size(); ++i) {
          const PathVar& pathVar = model.paths[k_id][i];
          const double vsval = model.cplex.getValue(pathVar.var);
          if (DBL_EQL(1.0, vsval))
            break;
          if (!DBL_EQL(0.0, vsval)) {
            fractional = true;
            pathVarsValuesQ.emplace(std::min(vsval, 1 - vsval), i);
          }
        }
        // if there are non integer path variables in P^k
        if (!pathVarsValuesQ.empty()) {
          assert(pathVarsValuesQ.size() > 1);
          size_t p1 = pathVarsValuesQ.top().second;
          pathVarsValuesQ.pop();
          size_t p2 = pathVarsValuesQ.top().second;
          
          assert(p1 != p2);

          // p1 is the smallest
          if (model.paths[k_id][p1].path.size() > 
              model.paths[k_id][p2].path.size())
            std::swap(p1, p2);

          const vector<Edge>& path1 = model.paths[k_id][p1].path.path_;
          const vector<Edge>& path2 = model.paths[k_id][p2].path.path_;

#ifdef VERBOSE
          cout << "branch on z[i] in P^" << k_id << endl;
          cout << "sizes: " << path1.size() << " and " << path2.size() << endl;
          //print(model.paths[k_id][p1].path, substrate);
          //print(model.paths[k_id][p2].path, substrate);
#endif
          // branch on first edge that is present in p1 but not in p2
          for (int i = 0; i < path1.size(); i++) {
            size_t id1 = substrate[path1[i]].id,
                   id2 = substrate[path2[i]].id;
            assert(id2 < num_edges(substrate));
            assert(id1 < num_edges(substrate));
            if (id1 != id2) {
              tree.push(new BPEdgeNode(
                  current.depth + 1, current_i, k_id, id1, current.dual,
                  current.infeasibilities));
              tree.push(new BPEdgeNode(
                  current.depth + 1, current_i, k_id, id2, current.dual,
                  current.infeasibilities));
              break;
            }
          }
        }
        if (fractional)
          break;
      }

      if (!fractional) {
        if (firstint_time < 0)
          firstint_time = timer.total();
        integer_solutions_found++;
        out_mapping = model.getCurrentMapping(substrate, virtual1);
#ifdef VERBOSE
        cout << " --- found integer solution with cost = "
            << out_mapping.cost << endl;
        print(out_mapping);
#endif
        ub = out_mapping.cost;
        cout << "ub:"<< ub << " > " << (timer.total() / 1000.0) <<endl;
        // cout << "new ub " << ub << " current " << current.dual << endl;
        assert(current.dual == ub);
        if (DBL_EQL(ub, lb))
          break;
      }
      select_time += profiler.total();
    }
  //}
  timer.stop();
  cout << "@sel_time:" << (select_time / 1000.0) << endl;
  cout << "@nodes_explored:" << nodes_explored << endl;
  bool maxCost = out_mapping.cost >= numeric_limits<Bandwidth>::max();
  cout << "@infeasible:"
    << (!timer.reachedTimeLimit() && maxCost) << endl;
  cout << "@optimal:" << (!timer.reachedTimeLimit() && !maxCost) << endl;
  if (!maxCost) {
    cout << "@cost:" << out_mapping.cost << endl;
    cout << "@firstint:" << firstint_time / 1000.0 << endl;
  }
  cout << "@time:" << (timer.total() / 1000.0) << endl;
  /*cout << "@v_explored:" << vmap_explored << endl;
  cout << "@e_explored:"
      << (nodes_explored - vmap_explored - 1) << endl;*/
  cout << "@ifound:" << integer_solutions_found << endl;
  const double cg_time {cplex_time / 1000.0};
  cout << std::setprecision(3) << "@cg_time:"  << cg_time  << endl;
  cout << "@addc_time:" << (addc_time / 1000.0) / cg_time << endl;
  cout << "@primal_time:" << (primal_time / 1000.0) / cg_time << endl;
  cout << "@pricing_time:" << (pricing_time / 1000.0) / cg_time << endl;
  cout << "@aux_time:" << (aux_time / 1000.0) / cg_time << endl;
  env.end();
  if (out_mapping.cost < numeric_limits<Bandwidth>::max())
    return true;
  else
    return false;
}

template bool vne::bprice<vne::BPStrategy::DFS>(Graph&,
    const Graph&, const Parameters, Mapping&);
template bool vne::bprice<vne::BPStrategy::BFS>(Graph&,
    const Graph&, const Parameters, Mapping&);
template bool vne::bprice<vne::BPStrategy::BEST_PROJECTION>(Graph&,
    const Graph&, const Parameters, Mapping&);
template bool vne::bprice<vne::BPStrategy::DFS,
         &vne::branch_on_nodes_dichotomic>(Graph&,
    const Graph&, const Parameters, Mapping&);
template bool vne::bprice<vne::BPStrategy::BFS,
         &vne::branch_on_nodes_dichotomic>(Graph&,
    const Graph&, const Parameters, Mapping&);
template bool vne::bprice<vne::BPStrategy::BEST_PROJECTION,
         &vne::branch_on_nodes_dichotomic>(Graph&,
    const Graph&, const Parameters, Mapping&);

