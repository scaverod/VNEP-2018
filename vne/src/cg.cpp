#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <queue>

#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include <ilcplex/ilocplex.h>
#pragma GCC diagnositc pop

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/lookup_edge.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>

#include "../include/cg.h"
//#include "../include/PriorityQueue.h"
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
using vne::PathVar;
using vne::CplexValue;

extern float primal_time, pricing_time, addc_time, aux_time; 

CplexValue cg_iteration(IloEnv& env, Graph& aux,
    Graph& substrate, const Graph& virtual1,
    const Timer<std::chrono::milliseconds>& global_timer,
    vne::CGModel& cg_model, vector<Path>& paths) {
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

// creates num_vertices(virtual1) vertices in an auxiliary graph
// vertex n+v corresponds to the vertex v
void buildImprovedAuxiliaryGraph(const Graph& substrate, const Graph& virtual1,
  DirGraph& aux_out) {
  // copy substrate graph nodes, create extra nodes
  for (Vertex s = 0; s < num_vertices(substrate); s++) {
    add_vertex(aux_out);
    add_vertex(aux_out);
  }
  const size_t S = num_vertices(substrate);
  // copy substrate links
  EdgesIterator e, eend;
  for(tie(e, eend) = edges(substrate); e != eend; e++) {
    add_edge(source(*e, substrate), target(*e, substrate), aux_out);
    add_edge(target(*e, substrate), source(*e, substrate), aux_out);
    add_edge(source(*e, substrate), S + target(*e, substrate), aux_out);
    add_edge(target(*e, substrate), S + source(*e, substrate), aux_out);
  }

  const size_t V = num_vertices(aux_out);
  size_t edges = num_edges(substrate);
  for (Vertex v = 0; v < num_vertices(virtual1); v++) {
    add_vertex(aux_out);
    for (Vertex s = 0; s < S; ++s) {
      if(aux_out[s].cpu >= virtual1[v].cpu) {
        add_edge(V + v, s, aux_out);
        add_edge(S + s, V + v, aux_out);

        EdgeResources eres {numeric_limits<Bandwidth>::max(), edges++};
        //Edge f;
        //bool success;
        //tie(f, success) = add_edge(S+v, s, eres, aux_out);
      }
    }
  }
}

bool vne::cg(Graph& substrate, const Graph& virtual1, const Parameters param,
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

  DirGraph aux;
  buildImprovedAuxiliaryGraph(substrate, virtual1, aux);

  //cg_iteration(env, aux, substrate, virtual1, timer,
  // cg_model, paths);
  env.end();

  // the solution is not necessarily integral
  return false;
}

