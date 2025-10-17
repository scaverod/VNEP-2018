/*
// Copyright 2013 Leonardo Moura
#include <iostream>
#include <memory>
#include <stack>
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include <ilcplex/ilocplex.h>
#pragma GCC diagnositc pop
#include "../include/bb.h"
#include "../include/cplex.h"
ILOSTLBEGIN

using std::cerr;
using std::cout;
using std::endl;

typedef IloIntVarArray cplex_var_type;
typedef IloArray<cplex_var_type> cplex_vertex_mapping;
typedef IloArray<IloArray<IloArray<cplex_var_type>>> cplex_edge_mapping;

ILOSOLVECALLBACK1(UserSolve, cplex_vertex_mapping, M) {
  getEnv().out() << " ---------- Hello World ----------- " << endl;
  getEnv().out() << "teste ------------------- " << endl;
}

ILOBRANCHCALLBACK1(MyBranch, cplex_vertex_mapping, M) {
  if ( getBranchType() != BranchOnVariable ) {
    getEnv().out() << "another type of branch" << endl;
    return;
  }
  getEnv().out() << " ---------- BRANCH ----------- " << endl;
  IloNumArray m_s(getEnv());
  IntegerFeasibilityArray feasibility_s(getEnv());
  try {
    for (IloInt s = 0; s < M.getSize(); ++s) {
      getValues(m_s, M[s]);
      getFeasibilities(feasibility_s, M[s]);
      for (IloInt v = 0; v < M[s].getSize(); ++v) {
        if(feasibility_s[v] == Infeasible) {
          makeBranch(M[s][v] <= 0, getObjValue());
          makeBranch(M[s][v] >= 1, getObjValue());
          return;
        }
      }
    }
  } catch (...) {
    throw;
  }
  m_s.end();
  feasibility_s.end();
}

bool vne::bb(Graph& substrate, const Graph& virtual1, const Parameters param,
    Mapping& out_mapping) {
  IloEnv env;
  try {
    IloModel model(env);
    IloNumArray cpuv(env, num_vertices(virtual1));
    IloNumArray cpus(env, num_vertices(substrate));
    IloArray<IloNumArray> bv(env, num_vertices(virtual1));
    IloArray<IloNumArray> bs(env, num_vertices(substrate));
    cplex_vertex_mapping M(env, num_vertices(substrate));
    cplex_edge_mapping D(env, num_vertices(substrate));

    // set params and vars
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
      M[s] = cplex_var_type(env, num_vertices(virtual1), 0, 1);
      D[s] = IloArray<IloArray<cplex_var_type>>(env, num_vertices(substrate));
      for (Vertex w = 0; w < num_vertices(substrate); w++) {
        D[s][w] = IloArray<cplex_var_type>(env, num_vertices(virtual1));
        for (Vertex k = 0; k < num_vertices(virtual1); k++) {
          D[s][w][k] = cplex_var_type(env, num_vertices(virtual1), 0, 1); 
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
    IloExpr bandUsed(env);
    for (Vertex v = 1; v < num_vertices(virtual1); v++)
      for (Vertex l = 0; l < num_vertices(substrate); l++)
        for (Vertex m = 0; m < num_vertices(substrate); m++) {
          bandUsed += IloScalProd(D[l][m][v], bv[v]);
        }
    IloObjective obj(env);
    obj.setExpr(bandUsed);
    obj.setSense(IloObjective::Minimize); 
    model.add(obj);

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
    // * path constraints
    // * * each virtual link is mapped
    for (Vertex s = 0; s < num_vertices(substrate); s++) {
      for (Vertex k = 0; k < num_vertices(virtual1); k++) {
        for (Vertex l = 0; l < num_vertices(virtual1); l++) {
          IloExpr mapped_out_minus_in(env);

          for (Vertex j = 0; j < num_vertices(substrate); j++) {
            // Invalid edges are not used
            model.add(D[s][j][k][l] <= bs[s][j]);
            mapped_out_minus_in += D[s][j][k][l] - D[j][s][k][l];
          }
          model.add(mapped_out_minus_in == M[s][k] - M[s][l]);
        }
      }
    }
    // * bandiwidth constraints
    for (Vertex s = 1; s < num_vertices(substrate); s++) {
      for (Vertex j = 0; j < s; j++) {
        IloExpr sum(env);
        for (Vertex k = 0; k < num_vertices(virtual1); k++) {
          for (Vertex l = 0; l < k; l++) {
            sum += D[s][j][k][l] * bv[k][l];
          }
        }
        model.add(sum <= bs[s][j]);
      }
    }

    env.out() << "starting........." << endl;
    IloCplex cplex(env);
    cplex.use(UserSolve(env, M));
    cplex.use(MyBranch(env, M));
    cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
    // cplex.setParam(IloCplex::PreInd, 0);
    cplex.extract(model);
    cplex.solve();

    cplexResult2Mapping(substrate, virtual1, cplex, M, D, out_mapping);
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
*/
