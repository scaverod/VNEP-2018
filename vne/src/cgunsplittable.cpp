#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <unordered_map>

#include <boost/heap/binomial_heap.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/lookup_edge.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>

#include "../include/vnedefs.h"
#include "../include/cg.h"
#include "../include/PriorityQueue.h"

ILOSTLBEGIN

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::map;
using std::queue;
using std::stringstream;
using std::unordered_map;

using vne::Path;
using vne::CplexValue;

float primal_time {0}, pricing_time {0}, addc_time {0}, aux_time {0}; 

// TODO extract Path to path.h
size_t Path::size() const {
  return path_.size();
}

void Path::push_back(Edge e) {
  path_.push_back(e);
}

Bandwidth Path::cost(const Graph& g, const vector<vne::CplexValue>& cvec)
    const {
  Bandwidth pcost = 0;
  for (Edge e : path_) {
    const int id = g[e].id;
    pcost += vcost_ * (cvec[id] + 1);
  }
  return pcost;
}

bool Path::operator==(const Path& p) const {
  return v_ == p.v_ && w_ == p.w_ 
      && vs_ == p.vs_ && ws_ == p.ws_
      && path_ == p.path_;
}

void vne::print(const Path& path, const Graph& g) {
  cout << "[" << path.v_ << " to " << path.w_ << "] "
      << "<" << path.vs_ << ", " << path.ws_ << ">";
  for (Edge e : path.path_) {
    cout << " - (" << source(e, g) << ", "
    << target(e, g) << ") #" << g[e].id;
  assert(path.vcost_ <= g[e].bandwidth);
  }
  cout << endl;
}

bool vne::loopless_path(const Path& path, const Graph& g) {
  std::set<Vertex> vertices_in_path;
  bool loopless = true;

  for (Edge e : path.path_) {
    vertices_in_path.insert(target(e, g));
    if (vertices_in_path.find(source(e, g)) != vertices_in_path.end()) {
      loopless = false;
      break;
    }
  }
  
  if (!loopless)
    print(path, g);

  return loopless;
}

CplexValue vne::cgunsplittable_iteration(IloEnv& env, Graph& aux,
    Graph& substrate, const Graph& virtual1,
    const Timer<std::chrono::milliseconds>& global_timer,
    CGModel& cg_model, vector<Path>& paths) {
  int iterations {0};
  float cplex_time {0},  // time spent on the LP
        path_time {0};  // time spent on Dijkstra


  try {
    // init dual variables
    IloNumArray ye(env, num_edges(substrate));
    IloNumArray lambda(env, num_edges(virtual1));
    IloArray<IloNumArray> pi(env, num_vertices(virtual1));
    IloArray<IloNumArray> psi(env, num_edges(virtual1));

    for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
      pi[v] = IloNumArray(env, num_vertices(substrate));
    }
    for (size_t eid = 0; eid < num_edges(virtual1); eid++)
      psi[eid] = IloNumArray(env, num_edges(substrate));

    vector<CplexValue> cost_vec(num_edges(aux),
        numeric_limits<CplexValue>::max());
    Timer<std::chrono::milliseconds> timer_ms;

    // cplex fresh starts
    //cg_model.cplex.setParam(IloCplex::AdvInd, 0);

    do {
      if (global_timer.reachedTimeLimit())
        return numeric_limits<CplexValue>::max();
      iterations++;
      // add paths to model
      timer_ms.start();
      for (const Path& path : paths) {
#ifdef VERBOSE
        cout << path.v_ << " to " << path.vs_ << endl;
        cout << path.w_ << " to " << path.ws_ << endl;
#endif
        // print(path, aux);
        IloNumColumn col = cg_model.obj(path.size() * path.vcost_)
            + cg_model.edge_demands[path.id_](1) 
            + cg_model.flow_fx[path.v_][path.vs_](1)
            + cg_model.flow_fx[path.w_][path.ws_](1);
        // add each edge separately in edge_capacities
        for (Edge e : path.path_) {
          const size_t e_id = substrate[e].id;
          assert(e_id < num_edges(substrate));
          col += cg_model.edge_capacities[e_id](path.vcost_);
          col += cg_model.edges_hosted[path.id_][e_id](1);
        }
#ifdef CHECK_REP_PATHS
        for (const PathVar& pathVar : cg_model.paths[path.id_]) {
          //cout << cg_model.cplex.getValue(pathVar.var) << " " << endl;
          //print(pathVar.path, substrate);
          if (pathVar.path == path) {
            print(pathVar.path, substrate);
            print(path, substrate);
            cout << "error: REPEATED PATHS" << endl;
            assert(false);
          }
        }
#endif
        PathVar pv {path, IloNumVar(col)};
        pv.var.setLB(0.0);
        cg_model.model.add(pv.var);

        assert(pv.path.id_ < num_edges(virtual1)); 
        cg_model.paths[pv.path.id_].push_back(pv);
      }
      assert(cg_model.checkPaths(substrate));
      paths.clear();
#ifdef EXPORT_MODEL
      stringstream filename;
      filename << "cg" << iterations << ".lp";
      cg_model.cplex.exportModel(filename.str().c_str());
#endif
      addc_time += timer_ms.total();
      // solve lp'
      timer_ms.start();
      if (!cg_model.cplex.solve()) {
        //cout << "Failed to optimize LP" << endl;
        // cout << "Status: " << cg_model.cplex.getStatus() << endl;
        break;
      }
      cplex_time += timer_ms.total();
      primal_time += timer_ms.total();
#ifdef VERVOSE
      cout << "status: " << cg_model.cplex.getStatus() << endl;
      cout << "z:" << cg_model.cplex.getObjValue() << endl;
#endif
      timer_ms.start();
      cg_model.cplex.getDuals(ye, cg_model.edge_capacities);
      cg_model.cplex.getDuals(lambda, cg_model.edge_demands);
      for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
        cg_model.cplex.getDuals(pi[v], cg_model.flow_fx[v]);
      }
      for (size_t e_id = 0; e_id < num_edges(virtual1); e_id++)
        cg_model.cplex.getDuals(psi[e_id], cg_model.edges_hosted[e_id]);
      aux_time += timer_ms.total();
#ifdef PRINT_DUALS
      IloNumArray yv(env, num_vertices(virtual1));
      IloNumArray ys(env, num_vertices(substrate));
      /*for (size_t i = 0; i < cg_model.pathVars.size(); ++i) {
        cout << " " << cg_model.cplex.getValue(cg_model.pathVars[i]);
      }
      cout << endl;*/

      cg_model.cplex.getDuals(yv, cg_model.vmap);
      for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
          cout << "yv["<<v<<"] = " << yv[v] << endl; 
      }

      cg_model.cplex.getDuals(ys, cg_model.smap);
      for (Vertex s = 0; s < num_vertices(substrate); ++s) {
          cout << "ys["<<s<<"] = " << ys[s] << endl; 
      }

      for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
        for (Vertex s = 0; s < num_vertices(substrate); ++s) {
          cout << "-pi["<<v<<", "<<s<<"] = "
              << -static_cast<CplexValue>(pi[v][s]) << endl;
          cout << "x = " << cg_model.cplex.getValue(cg_model.x[v][s]) << endl;
        }
      }
      
      for (size_t k_id = 0; k_id < num_edges(virtual1); k_id++) {
        for (size_t e_id =0; e_id < num_edges(substrate); e_id++) {
          cout << "psi["<<k_id<<", "<<e_id<<"] =" << psi[k_id][e_id] << endl;
        }
      }

      IloNumArray slacks(env, num_edges(substrate));
      cg_model.cplex.getSlacks(slacks, cg_model.edge_demands);

      for (Vertex s = 0; s < num_edges(substrate); ++s) {
        //assert(ye[s] == cg_model.cplex.getDual(cg_model.edge_capacities[s])); 
        cout << "-ye[" << s << "] = " << -static_cast<CplexValue>(ye[s]) 
            << endl;
        cout << " cv = "
            << cg_model.cplex.getValue(cg_model.cap_violation[s]) << endl;
      }
      for (Vertex s = 0; s < num_edges(virtual1); ++s) {
        cout << "lambda[" << s << "] = " << lambda[s] << endl;
        cout << " slack = " << slacks[s] << endl;
      }
#endif
      Timer<std::chrono::milliseconds> dtimer;
      float djk_time {0};
      // change aux edge's weights
      timer_ms.start();
      // get paths in auxiliary graph

      EdgesIterator e, eend;  
      for (tie(e, eend) = edges(virtual1); e != eend; ++e) {
        const size_t id = virtual1[*e].id;
        Path path {id, source(*e, virtual1), 0,
          target(*e, virtual1), 0, virtual1[*e].bandwidth};
        // set auxiliary edges
        for (Vertex s = 0; s < num_vertices(substrate); ++s) {
          Edge f;
          bool edgeExists;
          tie(f, edgeExists) = lookup_edge(num_vertices(substrate)+path.v_, s,
              aux);
          if (edgeExists) {
            if (cg_model.isNodeFixed(s) && cg_model.node_fixed[s] != path.v_) {
              cost_vec[aux[f].id] = numeric_limits<CplexValue>::max();
              //cout << path.v_ << " - " << s << " blocked " << endl;
            } else {
              const CplexValue pivs = -static_cast<CplexValue>(pi[path.v_][s]);
              assert(DBL_GE(pivs, 0.0));
              cost_vec[aux[f].id] = max(0.0, pivs);
            }
          }
          tie(f, edgeExists) = lookup_edge(num_vertices(substrate)+path.w_, s,
              aux);
          if (edgeExists) {
            if (cg_model.isNodeFixed(s) && cg_model.node_fixed[s] != path.w_) {
              cost_vec[aux[f].id] = numeric_limits<CplexValue>::max();
              //cout << path.w_ << " - " << s << " blocked " << endl;
            } else {
              const CplexValue piws = -static_cast<CplexValue>(pi[path.w_][s]);
              assert(DBL_GE(piws, 0.0));
              cost_vec[aux[f].id] = max(0.0, piws);
            }
          }
        }
        // set substrate graph weights
        EdgesIterator f, fend;
        for (tie(f, fend) = edges(substrate); f != fend; ++f) {
            const size_t e_id = substrate[*f].id;
            const CplexValue yeneg = -static_cast<CplexValue>(ye[e_id]);
            const CplexValue psikeneg = -static_cast<CplexValue>(psi[id][e_id]);
            assert(DBL_GE(yeneg, 0.0) && DBL_GE(psikeneg, 0.0)); // check in the precision
            cost_vec[e_id] = max(0.0, yeneg) + max(0.0, psikeneg);
            assert(cost_vec[id] >= 0.0);
        }

        dtimer.start();
        const CplexValue max_cost = lambda[id];
        CplexValue cost = vne::min_path_aux_graph(substrate, aux, cost_vec,
            num_vertices(substrate) + path.v_,
            num_vertices(substrate) + path.w_, path.vs_, path.ws_,
            virtual1[*e].bandwidth, max_cost, path.path_);
        djk_time += dtimer.total();

        // remove auxiliary edges
        for (Vertex s = 0; s < num_vertices(substrate); ++s) {
          Edge f;
          bool edgeExists;
          tie(f, edgeExists) = lookup_edge(num_vertices(substrate)+path.v_, s,
              aux);
          if (edgeExists)
            cost_vec[aux[f].id] = numeric_limits<CplexValue>::max();
          tie(f, edgeExists) = lookup_edge(num_vertices(substrate)+path.w_, s,
              aux);
          if (edgeExists)
            cost_vec[aux[f].id] = numeric_limits<CplexValue>::max();
        }

        if (DBL_L(cost, max_cost)) {
          assert(loopless_path(path, substrate));
          //print(path, substrate);
          //assert(cost > 0);
          assert(path.cost(substrate, cost_vec) < lambda[id]);
          paths.push_back(path);
        }
      }
      path_time += timer_ms.total();
      pricing_time += timer_ms.total() - djk_time;
      aux_time += djk_time;
    } while (!paths.empty());

    const bool feasible = cg_model.cplex.isPrimalFeasible()
        && cg_model.isFeasible(substrate, virtual1);
#ifdef PER_NODE_FULL_OUTPUT
    cout << "@cgiterations:" << iterations << endl;
    const float cgtime = global_timer.total();
    cout << "@cglptime:" << (cplex_time / 1000.0)  / static_cast<float>(cgtime)
        << endl;
    cout << "@cgpathtime:" << (path_time / 1000.0) / static_cast<float>(cgtime)
        << endl;
    cout << "@cgmodelfeasible:" << feasible << endl;
    if (feasible) {
      for(int id = 0; id < cg_model.paths.size(); id++) {
        for (const PathVar& pathVar : cg_model.paths[id]) {
          const CplexValue val =
              static_cast<CplexValue>(cg_model.cplex.getValue(pathVar.var));
          cout << id << ") " << val << " ";
          print(pathVar.path, substrate);
          assert(val >= 0.0);
        }
      }

      int xfrac = 0;
      for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
        for (Vertex s = 0; s < num_vertices(substrate); ++s) { 
          const double vsval = cg_model.cplex.getValue(cg_model.x[v][s]);
          xfrac += !( DBL_EQL(0.0, vsval) || DBL_EQL(1.0, vsval));
        }
      }
      cout << "@xfrac:"<<xfrac<<endl;
      int ffrac = 0;
      for (const auto& pathVars_k : cg_model.paths) {
        for (const auto& pathVar : pathVars_k) {
          const double vsval = cg_model.cplex.getValue(pathVar.var);
          ffrac += !( DBL_EQL(0.0, vsval) || DBL_EQL(1.0, vsval));
        }
      }
      cout << "@ffrac:"<<ffrac<<endl;

      if (global_timer.reachedTimeLimit()) {
        cout << "@cguFeasible:1" << endl;
        cout << "@cguCost:";
      } else {
        cout << "@cgunoptimal:1" << endl;
        cout << "@cgubestCost:";
      }
      cout << cg_model.cplex.getObjValue() << endl;
      return cg_model.cplex.getObjValue();
    }
#endif
    if (feasible) {
      return cg_model.cplex.getObjValue();
    }
  } catch (IloException& ex) {
    cerr << "Error " << ex << endl;
    throw(ex);
  } 
  return numeric_limits<CplexValue>::max();
}

bool vne::cgunsplittable(Graph& substrate, const Graph& virtual1,
    const Parameters param,
    Mapping& outMapping) {
  Timer<std::chrono::milliseconds> timer(param.timelimit_in_s * 1000);
  timer.start();
  IloEnv env;
  const CplexValue M {static_cast<CplexValue>(num_edges(virtual1))};
  const CplexValue K {static_cast<CplexValue>(
      num_edges(virtual1) * sum_bandwidth(substrate))};

  float cplex_time {0},  // time spent on the LP
         path_time {0};   // time spent on Dijkstra
  vector<Path> paths;
  if (!vne::get_initial_paths(substrate, virtual1, paths)) {
    cout << "infeasible" << endl;
    return false;
  }
  cout << "@initialColumns:"<< paths.size() << endl;

  CGModel cg_model(env, M, K);
  cg_model.init(env, substrate, virtual1);

  Graph aux;
  buildAuxiliaryGraph(substrate, virtual1, aux);

  vne::cgunsplittable_iteration(env, aux, substrate, virtual1, timer,
      cg_model, paths);
  env.end();

  // the solution is not necessarily integral
  return false;
}

bool vne::find_min_path(const Graph& g, const Bandwidth demand, const Vertex s,
    const Vertex t, vector<Edge>& out_path) {
  std::vector<bool> visited(num_vertices(g));
  std::vector<Edge> pred(num_vertices(g));
  queue<Vertex> Q;
  Q.push(s);
  visited[s] = true;
  while (!Q.empty()) {
    Vertex r = Q.front();
    if (r == t) {
      while (r != s) {
        out_path.push_back(pred[r]);
        r = source(pred[r], g);
      }
      return true;
    }
    Q.pop(); 
    EdgeIterator e, eend;
    for (tie(e, eend) = out_edges(r, g); e != eend; ++e) {
      if (g[*e].bandwidth < demand)
        continue;
      const Vertex y = target(*e, g);
      if (!visited[y]) {
        visited[y] = true;
        pred[y] = *e;
        Q.push(y);
      }
    }
  }
  return false;
}

bool vne::get_initial_paths(const Graph& substrate, const Graph& virtual1,
    vector<Path>& out_paths) {
  typedef pair<int, Vertex> cpu_index;
  vector<cpu_index> s_nodes;
  s_nodes.reserve(num_vertices(substrate));
  for (Vertex s {0}; s < num_vertices(substrate); ++s) {
    s_nodes.push_back(make_pair(substrate[s].cpu, s));
  }  
  std::sort(s_nodes.begin(), s_nodes.end(), std::greater<cpu_index>());
  vector<cpu_index> v_nodes;
  v_nodes.reserve(num_vertices(virtual1));
  for (Vertex v {0}; v < num_vertices(virtual1); ++v) {
    v_nodes.push_back(make_pair(substrate[v].cpu, v));
  }
  std::sort(v_nodes.begin(), v_nodes.end(), std::greater<cpu_index>());

  map<Vertex, Vertex> mapping;
  for (Vertex v {0}; v < num_vertices(virtual1); ++v) {
    mapping[v_nodes[v].second] = s_nodes[v].second;
  }

  EdgesIterator e, eend;
  for (tie(e, eend) = edges(virtual1); e != eend; ++e) {
    const Vertex v = source(*e, virtual1);
    const Vertex w = target(*e, virtual1);
    Vertex vs = mapping[v];
    Vertex ws = mapping[w];
    Path path {virtual1[*e].id, v, mapping[v], w, mapping[w],
        virtual1[*e].bandwidth};
    find_min_path(substrate, virtual1[*e].bandwidth, path.vs_, path.ws_,
        path.path_);
    if (path.path_.size() > 0) {
      out_paths.push_back(path);
    } else {
      // greedy does not work try all vertices
      for (Vertex vsi = 0; vsi < num_vertices(substrate); vsi++) {  
        path.vs_ = s_nodes[vsi].second;
        if (substrate[path.vs_].cpu < virtual1[v].cpu) {
          break;
        }
        for (Vertex wsi = 0; wsi < vsi; wsi++) {
          path.ws_ = s_nodes[wsi].second;
          if (substrate[path.ws_].cpu < virtual1[w].cpu) {
            break;
          }
          find_min_path(substrate, virtual1[*e].bandwidth, path.vs_, path.ws_,
              path.path_);
          if (path.path_.size() > 0) {
            out_paths.push_back(path);
            break;
          }
        }
        if (path.path_.size() > 0)
          break;
      }
      if (path.path_.size() == 0)
        return false; 
    }
  }
  return true;
}

bool vne::get_initial_paths_from_aux(const Graph& substrate, const Graph& virtual1,
    const Graph& aux, std::vector<Path>& paths) {
  vector<CplexValue> cost_vec(num_edges(aux), 1);
  std::fill(cost_vec.begin() + num_edges(substrate), cost_vec.end(),
      numeric_limits<CplexValue>::max());

  EdgesIterator e, eend;
  for (tie(e, eend) = edges(virtual1); e != eend; ++e) {
    Path path {virtual1[*e].id, source(*e, virtual1), 0,
        target(*e, virtual1), 0, virtual1[*e].bandwidth};
    for (Vertex s = 0; s < num_vertices(substrate); ++s) {
      Edge f;
      bool edgeExists;
      tie(f, edgeExists) = lookup_edge(num_vertices(substrate)+path.v_, s,
          aux);
      if (edgeExists) {
          cost_vec[aux[f].id] = 0;
      }
      tie(f, edgeExists) = lookup_edge(num_vertices(substrate)+path.w_, s,
          aux);
      if (edgeExists) {
          cost_vec[aux[f].id] = 0;
      }
    }
    CplexValue cost = vne::min_path_aux_graph(substrate, aux, cost_vec,
        num_vertices(substrate) + path.v_,
        num_vertices(substrate) + path.w_, path.vs_, path.ws_,
        path.vcost_, numeric_limits<CplexValue>::max(), path.path_);

    std::fill(cost_vec.begin() + num_edges(substrate), cost_vec.end(),
        numeric_limits<CplexValue>::max());
    if (path.size() == 0)
      return false;

    paths.push_back(path);
  }
  return true;
}

// creates num_vertices(virtual1) vertices in an auxiliary graph
// vertex n+v corresponds to the vertex v
void vne::buildAuxiliaryGraph(const Graph& substrate, const Graph& virtual1,
  Graph& aux_out) {
  aux_out = Graph(substrate);
  const size_t S = num_vertices(substrate);
  size_t edges = num_edges(substrate);
  for (Vertex v = 0; v < num_vertices(virtual1); v++) {
    add_vertex(aux_out);
    for (Vertex s = 0; s < S; ++s) {
      if(aux_out[s].cpu >= virtual1[v].cpu) {
        EdgeResources eres {numeric_limits<Bandwidth>::max(), edges++};
        Edge f;
        bool success;
        tie(f, success) = add_edge(S+v, s, eres, aux_out);
      }
    }
  }
}

void vne::buildAuxiliaryGraphRestrict(const Graph& substrate,
    const Graph& virtual1, Graph& aux_out) {
  aux_out = Graph(substrate);
  const size_t S = num_vertices(substrate);
  size_t edges = num_edges(substrate);

  vector<Bandwidth> sflow(num_vertices(substrate));
  for (Vertex s = 0; s < S; ++s)
    sflow[s] = out_bandwidth(substrate, s);

  for (Vertex v = 0; v < num_vertices(virtual1); v++) {
    const Bandwidth vflow = out_bandwidth(virtual1, v); 
    add_vertex(aux_out);
    for (Vertex s = 0; s < S; ++s) {
      if(aux_out[s].cpu >= virtual1[v].cpu && sflow[s] >= vflow) {
        EdgeResources eres {numeric_limits<Bandwidth>::max(), edges++};
        add_edge(S+v, s, eres, aux_out);
      }
    }
  }
}

CplexValue vne::min_path_aux_graph(const Graph& substrate,
    const Graph& aux, const vector<CplexValue>& cost_vec,
    const Vertex s, const Vertex t, Vertex& ss,
    Vertex& ts, const Bandwidth demand, const CplexValue max_cost,
    vector<Edge>& out_path) {
  typedef tuple<CplexValue, Vertex, Vertex> VertexNode;
  // min heap
  typedef boost::heap::binomial_heap<VertexNode,
      boost::heap::compare<std::greater<VertexNode>>> MinHeap;
  typedef MinHeap::handle_type Handle;

  MinHeap Q;
  vector<map<Vertex, Handle>> handles(num_vertices(substrate));
  vector<map<Vertex, pair<bool, Edge>>> visited_pred(num_vertices(substrate));
  
  // add candidate substrate nodes to the queue
  EdgeIterator e, eend;
  for (tie(e,eend) = out_edges(s, aux); e != eend; ++e) {
    const size_t id = aux[*e].id;
    const Vertex v = target(*e, aux);
    visited_pred[v][s].first = true;
    if (!DBL_EQL(cost_vec[id], numeric_limits<CplexValue>::max()))
      Q.push(make_tuple(cost_vec[id], v, v));
  }
  // Dijkstra
  while (!Q.empty()) {
    const Vertex u = get<1>(Q.top());
    const Vertex origin = get<2>(Q.top());
    const CplexValue cost = get<0>(Q.top());
    //cout << "visitou " << u << " com origem " << origin << endl;
    if (u == t) {
      // recover path
      Edge e = visited_pred[origin][t].second;
      ts = source(e, aux);
      do {
        e = visited_pred[origin][source(e, aux)].second;
        out_path.push_back(e);
        //cout << source(e, aux) << " -- " << target(e, aux) << endl;
      } while (source(e, aux) != origin);
      ss = source(e, aux);
      //cout << s - num_vertices(substrate) << " mapped to " << ss << " cost = "
      //    << endl;
      //cout << t - num_vertices(substrate) << " mapped to " << ts << endl;
      return cost;
    }
    if (cost > max_cost)
      return numeric_limits<CplexValue>::max();
    visited_pred[origin][u].first = true;
    Q.pop(); 
    // for every neighbor v
    for (tie(e,eend) = out_edges(u, aux); e != eend; ++e) {
      const size_t id =aux[*e].id;
      const Vertex v = target(*e, aux);
      if (DBL_G(cost_vec[id], max_cost)
          || visited_pred[origin][v].first
          || aux[*e].bandwidth < demand 
          || (v == t && u == origin))
        continue;
      //cout << " * " << v << " (" << cost_vec[id] << ")" << endl;
      const CplexValue new_cost = cost + ( (v == t) ? cost_vec[id]
          : (demand * (1.0 + cost_vec[id])));
      // a better cost can only be obtained while the node was not yet visited
      auto it = handles[origin].find(v);
      if (it == handles[origin].end()) {
        handles[origin][v] = Q.push(make_tuple(new_cost, v, origin));
        visited_pred[origin][v].second = *e;
        //cout << "init: " << source(*e, aux) << " -- " << target(*e, aux) << endl;
      } else {
        const CplexValue& cost = get<0>(*(it->second));
        if (new_cost < cost) {
          //cout << new_cost << " < " << cost << endl;
          //cout << "<decrease key>: ";
          //cout << source(*e, aux) << " -- " << target(*e, aux) << endl;
          Q.decrease(it->second, make_tuple(new_cost, v, origin));
          visited_pred[origin][v].second = *e;
        }
      }
    }
  }

  return numeric_limits<CplexValue>::max();
}

// TODO delete function ant its dependencies
CplexValue vne::min_path_aux_graph2(const Graph& substrate,
    const Graph& aux, const vector<CplexValue>& cost_vec,
    const Vertex s, const Vertex t, Vertex& ss,
    Vertex& ts, const Bandwidth demand, const CplexValue maxCost,
    vector<Edge>& out_path) {
  CplexValue* d;
  // 2i is the node 2i+1 is the node coming from an auxiliary node
  PriorityQueue Q(num_vertices(aux) * 2, &d);
  vector<bool> visited(num_vertices(aux)*2, false);
  vector<bool> discovered(num_vertices(aux)*2, false);
  vector<Vertex> pred(num_vertices(aux)*2);
  vector<Vertex> origin(num_vertices(aux)*2);

  d[t*2] = numeric_limits<CplexValue>::max();
  d[s*2] = 0;
  discovered[s*2] = true;
  // forbid going back
  origin[s*2] = numeric_limits<Vertex>::max();
  Q.insert(s*2);

  //cout << " map " << s << " to " << t << endl;
  //cout << "demand = " << demand << endl;

  while (!Q.empty()) {
    Vertex u = Q.key_top()/2;
    if (u == t) {
      break;
    }
    CplexValue baseCost = d[Q.key_top()];
    if (baseCost > maxCost)
      break;
    cout << "visitou " << u << endl;
    cout << "k: " << Q.key_top() << endl;
    cout << "d: " << d[Q.key_top()] << endl;
    bool fromAux = Q.key_top() % 2;
    visited[Q.key_top()] =  true;
    EdgeIterator e, eend;
    for (tie(e,eend) = out_edges(u, aux); e != eend; ++e) {
      const size_t id =aux[*e].id;
      if (DBL_EQL(cost_vec[id], numeric_limits<CplexValue>::max()) ||
          aux[*e].bandwidth < demand)
        continue;

      Vertex v = target(*e, aux);
      cout << " - v = " << v << endl;
      cout << " - cost_vec[" << id << "] = " << cost_vec[id] << endl;
      cout << " - band = " << aux[*e].bandwidth << endl;
      Vertex vindex = v*2 + (u==s);

      if (visited[vindex]) 
        continue;

      if (fromAux && v == t)
        continue;
      cout << "check origin: " << origin[Q.key_top()] << endl; 
      if (origin[Q.key_top()] == v)
        continue;

      CplexValue newCost;
      if (u == s || v == t)
        newCost = baseCost + cost_vec[id];
      else
        newCost = baseCost + demand * (1.0 + cost_vec[id]);

      cout << "dist " << u << ", " << v << " - " << newCost << endl;
      //assert(DBL_GE(newdist, 0.0));
      //cout << "  - " << vindex << " - " << newdist << endl;
      //cout << "  - origin vindex = " << origin[vindex] << endl;
      //cout << "  - origin u      = " << origin[u] << endl;
      if (!discovered[vindex]) {
        d[vindex] = newCost;
        Q.insert(vindex);
        discovered[vindex] = true;
        pred[vindex] = Q.key_top();
        if (u == s)
          origin[vindex] = v;
        else
          origin[vindex] = origin[Q.key_top()];
      } else if (d[vindex] > newCost) {
        d[vindex] = newCost;
        Q.decrease(vindex);
        pred[vindex] = Q.key_top();
        origin[vindex] = origin[Q.key_top()];
      }
    }
    Q.pop();
  }

  if (d[t*2] >= numeric_limits<CplexValue>::max())
    return numeric_limits<CplexValue>::max();

  Vertex v = pred[t*2];
  ts = v/2;
  while ((pred[v]/2) != s) {
    //cout << " v = " << v/2 << " - > " << pred[v]/2 << endl;
    out_path.push_back(lookup_edge(v/2, pred[v]/2, substrate).first);
    v = pred[v];
  }
  ss = v/2;
  cout << s - num_vertices(substrate) << " mapped to " << ss << " cost = "
      << d[t*2] << endl;
  cout << t - num_vertices(substrate) << " mapped to " << ts << endl;
  return d[t*2];
} 

void vne::CGModel::init(IloEnv& env, const Graph& substrate,
      const Graph& virtual1) {
    paths = vector<vector<PathVar>>(num_edges(virtual1));
    node_fixed = vector<Vertex>(num_vertices(substrate),
          numeric_limits<Vertex>::max());

    // optimization
    env.setNormalizer(IloFalse);
    
    // init x decision variables and constraints (they do not change)
    x = IloArray<DecisionVariableArray>(env, num_vertices(virtual1));

    vmap = IloRangeArray(env, num_vertices(virtual1), 1, 1);

    for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
      x[v] = DecisionVariableArray(env, num_vertices(substrate), 0, 1);
    //  model.add(IloSum(x[v]) == 1);
      vmap[v].setExpr(IloSum(x[v]));
    }
    model.add(vmap);

    smap = IloRangeArray(env, num_vertices(substrate), -IloInfinity, 1);
    for (Vertex s = 0; s < num_vertices(substrate); ++s) {
      IloExpr sum_s(env);
      sum_s += 0;
      for(Vertex v = 0; v < num_vertices(virtual1); ++v)
        sum_s += x[v][s];
      //model.add(sum_s <= 1);
      smap[s].setExpr(sum_s);
    }
    model.add(smap);

    // bandwidth demands
    dem_violation = IloNumVarArray(env, num_edges(virtual1));
    edge_demands = IloRangeArray(env, num_edges(virtual1));
    EdgesIterator e, eend;
    for (tie(e, eend) = edges(virtual1); e != eend; ++e) {
      const size_t id = substrate[*e].id;
      edge_demands[id] = IloRange(env, 1, 1);

      IloNumColumn col = obj(K) + edge_demands[id](1);
      dem_violation[id] = IloNumVar(col);
      dem_violation[id].setBounds(0, IloInfinity);
      model.add(dem_violation[id]);
    }
    edge_demands.setNames("demand");
    model.add(edge_demands);

    // bandwidth capacities
    // TODO use get(edge_index, g, edge_descriptor);
    IloNumArray capacities(env, num_edges(substrate));
    for (tie(e, eend) = edges(substrate); e != eend; ++e) {
      const size_t id = substrate[*e].id;
      capacities[id] = substrate[*e].bandwidth;
    }
    edge_capacities = IloRangeArray(env, -IloInfinity, capacities);
    edge_capacities.setNames("capacity");
    model.add(edge_capacities);

    // add capactiy violations to the edges to obtain a valid initial solution
    cap_violation = IloNumVarArray(env, num_edges(substrate));
    for (tie(e, eend) = edges(substrate); e != eend; ++e) {
      const size_t id = substrate[*e].id;
      IloNumColumn col = obj(K) + edge_capacities[id](-1);
      cap_violation[id] = IloNumVar(col);
      cap_violation[id].setBounds(0, IloInfinity);
      model.add(cap_violation[id]);
    }
    
    // no flow through invalid paths
    flow_fx = IloArray<IloRangeArray>(env, num_vertices(virtual1));
    for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
      flow_fx[v] = IloRangeArray(env, num_vertices(substrate), -IloInfinity, 0);
      for (Vertex s = 0; s < num_vertices(substrate); ++s) {
        flow_fx[v][s].setExpr(-x[v][s] * 
            static_cast<float>(in_degree(v, virtual1)));
      }
      model.add(flow_fx[v]);
    }
    
    edges_hosted = IloArray<IloRangeArray>(env, num_edges(virtual1));
    EdgesIterator k, kend;
    for (tie(k, kend) = edges(virtual1); k != kend; ++k) {
      const size_t k_id = virtual1[*k].id;
      assert(k_id < num_edges(virtual1));
      edges_hosted[k_id] = IloRangeArray(env, num_edges(substrate));
      EdgesIterator e, eend;
      for (tie(e, eend) = edges(substrate); e != eend; ++e) {
        const size_t e_id = substrate[*e].id;
        assert(e_id < num_edges(substrate));
        edges_hosted[k_id][e_id] = IloRange(env, -IloInfinity, 1);
      }
      edges_hosted[k_id].setNames("hosted");
      model.add(edges_hosted[k_id]);
    }
}

bool vne::CGModel::isFeasible(const Graph& substrate, const Graph& virtual1) 
      const {
  EdgesIterator e, eend;
  for (tie(e, eend) = edges(substrate); e != eend; ++e) {
    const size_t id = substrate[*e].id;
    if (!DBL_EQL(cplex.getValue(cap_violation[id]), 0.0))
      return false;
  } 
  for (tie(e, eend) = edges(virtual1); e != eend; ++e) {
    const size_t id = virtual1[*e].id;
    if (!DBL_EQL(cplex.getValue(dem_violation[id]), 0.0))
      return false;
  } 
  return true;
}

bool vne::CGModel::checkPaths(const Graph& substrate) const {
  for (const auto& paths_k : paths) {
    for (const auto& pathVar : paths_k) {
      for (const Edge& edge : pathVar.path.path_) {
        if (substrate[edge].id >= num_edges(substrate)) { 
          cout << "wrong id " << substrate[edge].id << endl;
          return false;
        }
      }
    }
  }
  return true;
}

bool vne::CGModel::doesPathExistUsingAuxEdge(const size_t id, const Vertex v,
      const Vertex s) {
  for (const PathVar& pathvar : paths[id]) {
    const Path& path = pathvar.path;
    if ((path.v_ == v && path.vs_ == s 
          && auxEdgeCanBeUsed(path.w_, path.ws_))
        || (path.w_ == v && path.ws_ == s
          && auxEdgeCanBeUsed(path.v_, path.vs_)))
        return true;
  }
  return false;
}

// finds a path from s to w, return false iff no path is found
bool vne::CGModel::getPathThatCoversAuxEdge(const Graph& aux,
    const Graph& substrate, const CGModel& model, Path& path) {
  // the edges of the substrate graph have a weight one, so 
  // the path with the least hops is the minimum cost path
  vector<CplexValue> cost_vec(num_edges(aux), 1);
  // initially all auxiliary edges are removed
  std::fill(cost_vec.begin() + num_edges(substrate), cost_vec.end(),
      numeric_limits<CplexValue>::max());

  // the edge (v,s) is fixed
  cost_vec[aux[
      lookup_edge(num_vertices(substrate)+path.v_, path.vs_, aux).first].id]
      = 0;

  for (Vertex s = 0; s < num_vertices(substrate); ++s) {
    // if the vertex s is fixed on another virtual node, do not add it
    if (model.isNodeFixed(s) && model.node_fixed[s] != path.w_)
      continue;
    Edge f;
    bool edgeExists;
    tie(f, edgeExists) = lookup_edge(num_vertices(substrate)+path.w_, s,
        aux);
    if (edgeExists) {
        cost_vec[aux[f].id] = 0;
    }
  }
  CplexValue cost = vne::min_path_aux_graph(substrate, aux, cost_vec,
      num_vertices(substrate) + path.v_,
      num_vertices(substrate) + path.w_, path.vs_, path.ws_,
      path.vcost_, numeric_limits<CplexValue>::max(), path.path_);

  return path.size() > 0;
}

Mapping vne::CGModel::getCurrentMapping(const Graph& substrate,
    const Graph& virtual1) {
  Mapping m(substrate, virtual1);
  for (Vertex v = 0; v < num_vertices(virtual1); v++) {
    for (Vertex s = 0; s < num_vertices(substrate); s++)
      if (DBL_GE(cplex.getValue(x[v][s]),1.0))
        m.map_vertex(v,s);

    EdgeIterator e, eend;
    for (tie(e, eend) = out_edges(v, virtual1); e != eend; ++e) {
      const Vertex w = target(*e, virtual1);
      if (w > v)
        continue;
      const size_t k = virtual1[*e].id;
      
      for (const PathVar& pv : paths[k]) {
        if (DBL_GE(cplex.getValue(pv.var), 1.0)) {
          m.map_path(virtual1, *e, pv.path.path_);
          break;
        }
      }
    }
  }

  return m;
}

bool vne::CGModel::isNodeFixed(Vertex s) const {
  return node_fixed[s] != numeric_limits<Vertex>::max();
}

void vne::CGModel::fixNode(Vertex v, Vertex s) {
  x[v][s].setLB(1.0);
  node_fixed[s] = v;
}

void vne::CGModel::freeNode(Vertex v, Vertex s) {
  x[v][s].setLB(0.0);
  node_fixed[s] = numeric_limits<Vertex>::max();
}

bool vne::CGModel::auxEdgeCanBeUsed(Vertex v, Vertex s) const {
  return !isNodeFixed(s) || node_fixed[s] == v;
}
