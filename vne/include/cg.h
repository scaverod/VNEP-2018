#pragma once
#include <iostream>
#include <stack>

#include <ilcplex/ilocplex.h>
#include "../include/util.h"
#include "../include/vnedefs.h"

//#define PRINT_DUALS
//#define EXPORT_MODEL
//#define PER_NODE_FULL_OUTPUT
//#define CHECK_REP_PATHS
//#define VERBOSE

namespace vne{
typedef IloNumVarArray DecisionVariableArray;
typedef double CplexValue;

bool cg(Graph& substrate, const Graph& virtual1, const Parameters param,
    Mapping& outMapping);
bool cgunsplittable(Graph& substrate, const Graph& virtual1,
    const Parameters param, Mapping& outMapping);

struct Path {
  size_t id_;
  Vertex v_, vs_, w_, ws_;
  Bandwidth vcost_;
  std::vector<Edge> path_;

  size_t size() const; 
  void push_back(Edge e);
  Bandwidth cost(const Graph& g, const std::vector<CplexValue>& cvec) const; 
  bool operator==(const Path&) const;
};

struct PathVar {
    Path path;
    IloNumVar var;
};

struct CGModel {
  const CplexValue M;
  const CplexValue K;
  
  IloObjective obj;
  IloModel model;
  IloCplex cplex;

  // decision variables
  IloArray<DecisionVariableArray> x;
  std::vector<std::vector<PathVar>> paths;

  // constraints
  IloRangeArray vmap;
  IloRangeArray smap;

  IloRangeArray edge_demands;
  IloRangeArray edge_capacities;
  IloArray<IloRangeArray> edges_hosted;
  IloNumVarArray cap_violation;
  IloNumVarArray dem_violation;
  IloArray<IloRangeArray> flow_fx;

  std::vector<Vertex> node_fixed;

  CGModel(IloEnv& env, Bandwidth _M, Bandwidth _K) : M(_M), K(_K), obj(env),
      model(env), cplex(model) {
    model.add(obj);
  }
  /** create variables and add constraints to the model */
  void init(IloEnv& env, const Graph& substrate, const Graph& virtual1);
  bool isFeasible(const Graph& substrate, const Graph& virtual1) const;
  bool checkPaths(const Graph& substrate) const;
  bool doesPathExistUsingAuxEdge(const size_t id, const Vertex v,
      const Vertex s);
  bool getPathThatCoversAuxEdge(const Graph& aux,
      const Graph& substrate, const CGModel& model, Path& path);
  Mapping getCurrentMapping(const Graph& substrate, const Graph& virtual1);
  
  // cuts management
  bool isNodeFixed(Vertex s) const;
  bool auxEdgeCanBeUsed(Vertex v, Vertex s) const;
  void fixNode(Vertex v, Vertex s);
  void freeNode(Vertex v, Vertex s);
};

// debug methods
void print(const Path& path, const Graph& g);
bool loopless_path(const Path& path, const Graph& g);

// return false if some virtual link has no available path
bool get_initial_paths(const Graph&, const Graph&, std::vector<Path>&);
bool get_initial_paths_from_aux(const Graph& substrate, const Graph& virtual1,
    const Graph& aux, std::vector<Path>&);

void buildAuxiliaryGraph(const Graph& substrate, const Graph& virtual1,
    Graph& aux_out);
// restrict mappings based on flow
void buildAuxiliaryGraphRestrict(const Graph& substrate, const Graph& virtual1,
    Graph& aux_out);

bool find_min_path(const Graph& g, const Bandwidth demand, const Vertex s,
    const Vertex t, std::vector<Edge>& out_path);
CplexValue min_path_aux_graph(const Graph& substrate,
    const Graph& aux, const std::vector<CplexValue>& cost,
    const Vertex s, const Vertex t, Vertex& ss, Vertex& ts,
    const Bandwidth demand, const CplexValue maxCost, std::vector<Edge>& out_path);
CplexValue min_path_aux_graph2(const Graph& substrate,
    const Graph& aux, const std::vector<CplexValue>& cost,
    const Vertex s, const Vertex t, Vertex& ss, Vertex& ts,
    const Bandwidth demand, const CplexValue maxCost, std::vector<Edge>& out_path);
// the vector initial_paths is discarded
CplexValue cgunsplittable_iteration(IloEnv& env, Graph& aux, Graph& substrate,
    const Graph& virtual1, const Timer<std::chrono::milliseconds>& timer,
    CGModel& cg_model, std::vector<Path>& paths);
}
