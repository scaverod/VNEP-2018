// Copyright 2013 Leonardo Moura
#include <iostream>
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include <ilcplex/ilocplex.h>
#pragma GCC diagnositc pop

#include <boost/graph/lookup_edge.hpp>
#include "../include/fix.h"
ILOSTLBEGIN

using std::cerr;
using std::cout;
using std::endl;

typedef IloArray<IloIntVarArray> cplex_vertex_mapping;
typedef IloArray<IloArray<IloArray<IloIntVarArray>>> cplex_edge_mapping;


bool vne::fix(Graph& substrate, const Graph& virtual1, const Parameters param,
    Mapping& out_mapping) {
  IloEnv env;
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
    IloArray<IloArray<IloArray<IloIntVarArray>>> D(env, num_vertices(substrate));

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
      M[s] = IloIntVarArray(env, num_vertices(virtual1), 0, 1);
      D[s] = IloArray<IloArray<IloIntVarArray>>(env, num_vertices(substrate));
      for (Vertex w = 0; w < num_vertices(substrate); w++) {
        D[s][w] = IloArray<IloIntVarArray>(env, num_vertices(virtual1));
        for (Vertex k = 0; k < num_vertices(virtual1); k++) {
          D[s][w][k] = IloIntVarArray(env, num_vertices(virtual1), 0, 1); 
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

    // set objective function
    /*IloExpr bandUsed(env);
    for (Vertex v = 1; v < num_vertices(virtual1); v++)
      for (Vertex l = 0; l < num_vertices(substrate); l++)
        for (Vertex m = 0; m < num_vertices(substrate); m++)
          bandUsed += IloScalProd(D[l][m][v], bv[v]);
    */
    IloObjective obj(env);
    // obj.setExpr(bandUsed);
    obj.setSense(IloObjective::Minimize); 
    IloModel model(env);
    model.add(obj);

    // set constraints
    // * node mapping constraints
    // * * each substrate node host at most one virtual node
    IloExprArray substrateNodeUse(env, num_vertices(substrate));
    IloExprArray substrateNodeCpuUse(env, num_vertices(substrate));
    for (Vertex s = 0; s < num_vertices(substrate); s++) {
      substrateNodeUse[s]    = IloSum(M[s]);
      substrateNodeCpuUse[s] = IloScalProd(M[s], cpuv);
    }
    IloRangeArray substrateNodesUsedOnce(env, 0, substrateNodeUse, 1);
    IloRangeArray substrateNodesCpuUse(env, 0, substrateNodeCpuUse, cpus);
    model.add(substrateNodesUsedOnce);
    model.add(substrateNodesCpuUse);

    // * * each virtual node is mapped to exactly one substrate node
    IloExprArray virtualNodeUse(env, num_vertices(virtual1));
    for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
      virtualNodeUse[v] = IloExpr(env);
      for (Vertex s = 0; s < num_vertices(substrate); s++)
        virtualNodeUse[v] += M[s][v];
    }
    IloRangeArray virtualNodesUse(env, 1, virtualNodeUse, 1);
    model.add(virtualNodesUse);

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
    // IloExprArray substrateEdgeBandwidthUse(env,
    for (Vertex s = 1; s < num_vertices(substrate); s++) {
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
    IloCplex cplex(model);
    // cplex.exportModel("test.lp");
    cplex.setParam(IloCplex::ClockType, 2);
    cplex.setParam(IloCplex::TiLim, param.timelimit_in_s);
    if (!cplex.solve()) {
      env.error() << "Failed to optimize LP" << endl;
      throw(-1);
    }
    cout << "Status: " << cplex.getStatus() << endl;
    if (cplex.getStatus() != IloAlgorithm::Optimal) {
      if (cplex.getStatus() != IloAlgorithm::Feasible)
        cout << "CostBest:" << cplex.getObjValue() << endl;
      return false;
    }

    cplexResult2Mapping2(substrate, virtual1, cplex, M, D, out_mapping);
  } catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
  } catch (...) {
    cerr << "Unknown exception caught" << endl;
    env.end();
    return false;
  }
  env.end();
  return true;
}

