// Copyright 2013 Leonardo Moura
#include <iostream>
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include <ilcplex/ilocplex.h>
#pragma GCC diagnositc pop

#include <boost/graph/lookup_edge.hpp>
#include "../include/cplex.h"
ILOSTLBEGIN

using std::cerr;
using std::cout;
using std::endl;

ILOINCUMBENTCALLBACK1(IncumbentCallback, Timer<std::chrono::milliseconds>&, timer) {
  if (timer.isRunning())
    timer.stop();
}

template <typename VMAP, typename EMAP>
void cplexResult2Mapping(Graph& substrate, const Graph& virtual1, 
    IloCplex& cplex, VMAP& M, EMAP& D,
    Mapping& out_mapping) {
  using std::cout;
  using std::endl;
  Instance instance {substrate, virtual1};
  out_mapping = Mapping(instance);
  // out_mapping.cost = cplex.getObjValue();
  cout << "sol = " << out_mapping.cost << endl;
  for (Vertex s = 0; s < num_vertices(substrate); ++s) {
    for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
      if (cplex.getValue(M[s][v]) > 0.0) {
        out_mapping.vertex[v] = s;
        cout << v << " was mapped to " << s << " - > " <<
            cplex.getValue(M[s][v]) << endl;
      }
    }
  }
  int subedges {0};
  for (Vertex v = 1; v < num_vertices(virtual1); ++v) {
    for (Vertex k = 0; k < v; ++k) {
      bool exists;
      Edge virtual_edge;
      tie(virtual_edge, exists) = boost::edge(v, k, virtual1);
      if (!exists)
        continue;
      // map the path from v to k
      Vertex curr = out_mapping.vertex[v],
             dest = out_mapping.vertex[k];
      while (curr != dest) {
        for (Vertex s = 0; s < num_vertices(substrate); ++s) {
          if (cplex.getValue(D[curr][s][v][k])) {
            subedges++;
            out_mapping.cost += virtual1[virtual_edge].bandwidth;
            cout << " add " << v << " - " << k << endl;
             out_mapping.edge_map_ref(v, k)->push_back(
                boost::edge(s,curr,substrate).first);
            curr = s;
            break;
          }
        }
      }
    }
  }
  cout << "@cpedges:" << subedges << endl;
}

template <class T>
bool cplex_generic(Graph& substrate, const Graph& virtual1, const Parameters param,
    Mapping& out_mapping) {
  IloEnv env;
  typedef IloArray<T> cplex_vertex_mapping;
  typedef IloArray<IloArray<IloArray<T>>> cplex_edge_mapping;

#ifndef NDEBUG
  cout << "DEBUG on" << endl;
#endif
  // normalizing was taking a long time
  env.setNormalizer(IloFalse);
  Timer<std::chrono::milliseconds> timer;
  timer.start();
  try {
    // Virtual nodes' cpu demand / Substrate nodes' cpu capacity
    IloNumArray cpuv(env, num_vertices(virtual1));
    IloNumArray cpus(env, num_vertices(substrate));
    // Bandwidth
    IloArray<IloNumArray> bv(env, num_vertices(virtual1));
    IloArray<IloNumArray> bs(env, num_vertices(substrate));
    // M[s][v] = 1 if virtual node v is mapped to substrate node s, 0 otherwise
    cplex_vertex_mapping M(env, num_vertices(substrate));
    // D[s][r][u][v] = 1 if substrate link (s,r) hosts the virtual link (u,v)
    IloArray<IloArray<IloArray<T>>> D(env, num_vertices(substrate));

    // set parameters and initialize vars
    for (Vertex v = 0; v < num_vertices(virtual1); v++) {
      cpuv[v] = virtual1[v].cpu;
      bv[v] = (IloNumArray(env, num_vertices(virtual1)));
      EdgeIterator e, eend;
      for (std::tie(e, eend) = out_edges(v, virtual1); e != eend; e++) {
        const Vertex w = target(*e, virtual1);
        if (v > w) {
          bv[v][w] = virtual1[*e].bandwidth;
        }
      }
    }
    for (Vertex s = 0; s < num_vertices(substrate); s++) {
      M[s] = T(env, num_vertices(virtual1), 0, 1);
      D[s] = IloArray<IloArray<T>>(env, num_vertices(substrate));
      for (Vertex w = 0; w < num_vertices(substrate); w++) {
        D[s][w] = IloArray<T>(env, num_vertices(virtual1));
        for (Vertex k = 0; k < num_vertices(virtual1); k++) {
          D[s][w][k] = T(env, num_vertices(virtual1), 0, 1); 
        }
      }
      cpus[s] = substrate[s].cpu;
      bs[s] = IloNumArray(env, num_vertices(substrate));
      EdgeIterator e, eend;
      for (std::tie(e, eend) = out_edges(s, substrate); e != eend; e++) {
        const Vertex w = target(*e, substrate);
        bs[s][w] = substrate[*e].bandwidth;
      }
    }

    //cout << "obj fun" << endl;
    // set objective function
    IloExpr bandUsed(env);
    for (Vertex v = 1; v < num_vertices(virtual1); v++)
      for (Vertex l = 0; l < num_vertices(substrate); l++)
        for (Vertex m = 0; m < num_vertices(substrate); m++)
          bandUsed += IloScalProd(D[l][m][v], bv[v]);

    IloObjective obj(env);
    obj.setExpr(bandUsed);
    obj.setSense(IloObjective::Minimize); 
    IloModel model(env);
    model.add(obj);
    //cout << "node contraints" << endl;
    // set constraints
    // * node mapping constraints
    // * * each substrate node host at most one virtual node
    for (Vertex s = 0; s < num_vertices(substrate); s++) {
      model.add(IloScalProd(M[s], cpuv) <= cpus[s]);
      model.add(IloSum(M[s]) <= 1);
    }
    // * * each virtual node is mapped to exactly one substrate node
    for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
      IloExpr virtual_nodes(env);
      for (Vertex s = 0; s < num_vertices(substrate); s++)
        virtual_nodes += M[s][v];
      model.add(virtual_nodes == 1);
    }

    //cout << "path constraints" << endl;
    // * path constraints
    // * * each virtual link is mapped
    for (Vertex s = 0; s < num_vertices(substrate); s++) {
      for (Vertex k = 0; k < num_vertices(virtual1); k++) {
        for (Vertex l = 0; l < k; l++) {
          IloExpr mapped_out_minus_in(env);

          for (Vertex j = 0; j < num_vertices(substrate); j++) {
            // Invalid edges are not used
            // model.add(D[s][j][k][l] <= bs[s][j]);
            mapped_out_minus_in += D[s][j][k][l] - D[j][s][k][l];
          }
          model.add(mapped_out_minus_in == M[s][k] - M[s][l]);
        }
      }
    }
    // * bandwidth constraints
    for (Vertex s = 0; s < num_vertices(substrate); s++) {
      for (Vertex j = 0; j < num_vertices(substrate); j++) {
        IloExpr sum(env);
        for (Vertex k = 0; k < num_vertices(virtual1); k++) {
          for (Vertex l = 0; l < k; l++) {
            sum += (D[j][s][k][l] + D[s][j][k][l]) * bv[k][l];
          }
        }
        model.add(sum <= bs[s][j]);
      }
    }

    Timer<std::chrono::milliseconds> int_timer;

    IloCplex cplex(model);
    cplex.use(IncumbentCallback(env, int_timer));
    //cplex.setParam(IloCplex::PreInd, false);
    // cplex.exportModel("test.lp");
    cplex.setParam(IloCplex::IloCplex::Threads, 1);
    cplex.setParam(IloCplex::ClockType, 2);
    cplex.setParam(IloCplex::TiLim, param.timelimit_in_s);
    cplex.solve();
    timer.stop();

    cout << "@time:" << timer.total() / 1000.0 << endl;
    cout << "@optimal:"
        << (cplex.getStatus() == IloAlgorithm::Optimal) << endl;

    bool unknown = cplex.getStatus() == IloAlgorithm::Unknown;
    cout << "@unknown:" << unknown << endl;
    bool infeasible = cplex.getStatus() != IloAlgorithm::Optimal &&
        cplex.getStatus() != IloAlgorithm::Feasible;
    cout << "@infeasible:" << infeasible << endl;
    if (!infeasible)
      cout << "@firstint:" << int_timer.total() / 1000.0 << endl;
    if (infeasible || unknown)
      return false;

    cout << "@cost:"<< cplex.getObjValue() << endl;

    /*
    for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
      for (Vertex s = 0; s < num_vertices(substrate); ++s) {
        cout << "x[" << s << ", " << v << "] = " << cplex.getValue(M[s][v]) << endl;
      }
    }
    EdgesIterator e, eend;
    for (tie(e, eend) = edges(virtual1); e != eend; ++e) {
      Vertex v = source(*e, virtual1), w = target(*e, virtual1);
      EdgesIterator f, fend;
      if (v < w)
        swap(v, w);
      //for (tie(f, fend) = edges(substrate); f != fend; ++f) {
        //Vertex s = source(*f, substrate), t = target(*f, substrate);
        // if (s < t)
        //  swap(s, t);
      for (Vertex s = 0; s < num_vertices(substrate); ++s) {
        for (Vertex t = 0; t < num_vertices(substrate); ++t) {
          cout << v << ", " << w << " -> " << s << "," << t << " = ";
          cout << cplex.getValue(D[s][t][v][w]) << endl;
        }
      }
    }*/

  const bool relaxed_cplex = std::is_same<T, IloIntVarArray>::value;
  if (relaxed_cplex)
    cplexResult2Mapping(substrate, virtual1, cplex, M, D, out_mapping);
  } catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
  } catch (...) {
    cerr << "Unknown exception caught" << endl;
    cout << "@time:" << timer.total() / 1000.0 << endl;
    env.end();
    return false;
  }
  env.end();
  return true;
}

//template <class T>
// nothing done yet
bool cplex_improved(Graph& substrate, const Graph& virtual1, const Parameters param,
    Mapping& out_mapping) {
  IloEnv env;
  typedef IloIntVarArray T;
  typedef IloArray<T> cplex_vertex_mapping;
  typedef IloArray<IloArray<IloArray<T>>> cplex_edge_mapping;

#ifndef NDEBUG
  cout << "DEBUG on" << endl;
#endif
  // normalizing was taking a long time
  env.setNormalizer(IloFalse);
  Timer<std::chrono::milliseconds> timer;
  timer.start();
  try {
    // Virtual nodes' cpu demand / Substrate nodes' cpu capacity
    IloNumArray cpuv(env, num_vertices(virtual1));
    IloNumArray cpus(env, num_vertices(substrate));
    // Bandwidth
    IloArray<IloNumArray> bv(env, num_vertices(virtual1));
    IloArray<IloNumArray> bs(env, num_vertices(substrate));
    // M[s][v] = 1 if virtual node v is mapped to substrate node s, 0 otherwise
    cplex_vertex_mapping M(env, num_vertices(substrate));
    // D[s][r][u][v] = 1 if substrate link (s,r) hosts the virtual link (u,v)
    IloArray<IloArray<IloArray<T>>> D(env, num_vertices(substrate));

    // set parameters and initialize vars
    for (Vertex v = 0; v < num_vertices(virtual1); v++) {
      cpuv[v] = virtual1[v].cpu;
      bv[v] = (IloNumArray(env, num_vertices(virtual1)));
      EdgeIterator e, eend;
      for (std::tie(e, eend) = out_edges(v, virtual1); e != eend; e++) {
        const Vertex w = target(*e, virtual1);
        if (v > w) {
          bv[v][w] = virtual1[*e].bandwidth;
        }
      }
    }
    for (Vertex s = 0; s < num_vertices(substrate); s++) {
      M[s] = T(env, num_vertices(virtual1), 0, 1);
      D[s] = IloArray<IloArray<T>>(env, num_vertices(substrate));
      for (Vertex w = 0; w < num_vertices(substrate); w++) {
        D[s][w] = IloArray<T>(env, num_vertices(virtual1));
        for (Vertex k = 0; k < num_vertices(virtual1); k++) {
          D[s][w][k] = T(env, num_vertices(virtual1), 0, 1); 
        }
      }
      cpus[s] = substrate[s].cpu;
      bs[s] = IloNumArray(env, num_vertices(substrate));
      EdgeIterator e, eend;
      for (std::tie(e, eend) = out_edges(s, substrate); e != eend; e++) {
        const Vertex w = target(*e, substrate);
        bs[s][w] = substrate[*e].bandwidth;
      }
    }

    //cout << "obj fun" << endl;
    // set objective function
    IloExpr bandUsed(env);
    for (Vertex v = 1; v < num_vertices(virtual1); v++)
      for (Vertex l = 0; l < num_vertices(substrate); l++)
        for (Vertex m = 0; m < num_vertices(substrate); m++)
          bandUsed += IloScalProd(D[l][m][v], bv[v]);

    IloObjective obj(env);
    obj.setExpr(bandUsed);
    obj.setSense(IloObjective::Minimize); 
    IloModel model(env);
    model.add(obj);
    //cout << "node contraints" << endl;
    // set constraints
    // * node mapping constraints
    // * * each substrate node host at most one virtual node
    for (Vertex s = 0; s < num_vertices(substrate); s++) {
      model.add(IloScalProd(M[s], cpuv) <= cpus[s]);
      model.add(IloSum(M[s]) <= 1);
    }
    // * * each virtual node is mapped to exactly one substrate node
    for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
      IloExpr virtual_nodes(env);
      for (Vertex s = 0; s < num_vertices(substrate); s++)
        virtual_nodes += M[s][v];
      model.add(virtual_nodes == 1);
    }

    //cout << "path constraints" << endl;
    // * path constraints
    // * * each virtual link is mapped
    for (Vertex s = 0; s < num_vertices(substrate); s++) {
      for (Vertex k = 0; k < num_vertices(virtual1); k++) {
        for (Vertex l = 0; l < k; l++) {
          IloExpr mapped_out_minus_in(env);

          for (Vertex j = 0; j < num_vertices(substrate); j++) {
            // Invalid edges are not used
            // model.add(D[s][j][k][l] <= bs[s][j]);
            mapped_out_minus_in += D[s][j][k][l] - D[j][s][k][l];
          }
          model.add(mapped_out_minus_in == M[s][k] - M[s][l]);
        }
      }
    }
    // * bandwidth constraints
    for (Vertex s = 0; s < num_vertices(substrate); s++) {
      for (Vertex j = 0; j < num_vertices(substrate); j++) {
        IloExpr sum(env);
        for (Vertex k = 0; k < num_vertices(virtual1); k++) {
          for (Vertex l = 0; l < k; l++) {
            sum += (D[j][s][k][l] + D[s][j][k][l]) * bv[k][l];
          }
        }
        model.add(sum <= bs[s][j]);
      }
    }

    Timer<std::chrono::milliseconds> int_timer;

    IloCplex cplex(model);
    cplex.use(IncumbentCallback(env, int_timer));
    //cplex.setParam(IloCplex::PreInd, false);
    // cplex.exportModel("test.lp");
    cplex.setParam(IloCplex::IloCplex::Threads, 1);
    cplex.setParam(IloCplex::ClockType, 2);
    cplex.setParam(IloCplex::TiLim, param.timelimit_in_s);
    cplex.solve();
    timer.stop();

    cout << "@time:" << timer.total() / 1000.0 << endl;
    cout << "@optimal:"
        << (cplex.getStatus() == IloAlgorithm::Optimal) << endl;

    bool unknown = cplex.getStatus() == IloAlgorithm::Unknown;
    cout << "@unknown:" << unknown << endl;
    bool infeasible = cplex.getStatus() != IloAlgorithm::Optimal &&
        cplex.getStatus() != IloAlgorithm::Feasible;
    cout << "@infeasible:" << infeasible << endl;
    if (!infeasible)
      cout << "@firstint:" << int_timer.total() / 1000.0 << endl;
    if (infeasible || unknown)
      return false;

    cout << "@cost:"<< cplex.getObjValue() << endl;

  const bool relaxed_cplex = std::is_same<T, IloIntVarArray>::value;
  if (relaxed_cplex)
    cplexResult2Mapping(substrate, virtual1, cplex, M, D, out_mapping);
  } catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
  } catch (...) {
    cerr << "Unknown exception caught" << endl;
    cout << "@time:" << timer.total() / 1000.0 << endl;
    env.end();
    return false;
  }
  env.end();
  return true;
}

bool vne::cplex(Graph& substrate, const Graph& virtual1, const Parameters param,
    Mapping& out_mapping) {
  const bool result = cplex_generic<IloIntVarArray>(substrate, virtual1, param, out_mapping);
  return result;
}

bool vne::cplexrelax(Graph& substrate, const Graph& virtual1, const Parameters param,
    Mapping& out_mapping) {
  return cplex_generic<IloNumVarArray>(substrate, virtual1, param, out_mapping);
}
