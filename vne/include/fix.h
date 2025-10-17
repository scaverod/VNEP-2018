#pragma once 
#include <iostream>

#include <ilcplex/ilocplex.h>
#include "../include/util.h"
#include "../include/vnedefs.h"

namespace vne{
bool fix(Graph& substrate, const Graph& virtual1, const Parameters param, Mapping& outMapping);

template <typename VMAP, typename EMAP>
void cplexResult2Mapping2(Graph& substrate, const Graph& virtual1, 
    IloCplex& cplex, VMAP& M, EMAP& D,
    Mapping& out_mapping) {
  using std::cout;
  using std::endl;
  Instance instance {substrate, virtual1};
  out_mapping = Mapping(instance);
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
  for (Vertex v = 1; v < num_vertices(virtual1); ++v) {
    for (Vertex k = 0; k < v; ++k) {
      Edge e;
      bool exists;
      tie(e, exists) = boost::edge(v, k, virtual1);
      if(!exists)
        continue;
      // map the path from v to k
      Vertex curr = out_mapping.vertex[v],
             dest = out_mapping.vertex[k];
      while (curr != dest) {
        for (Vertex s = 0; s < num_vertices(substrate); ++s) {
          if (cplex.getValue(D[curr][s][v][k])) {
            cout << " add " << v << " - " << k << endl;
             out_mapping.edge_map_ref(v, k)->push_back(
                boost::edge(curr,s,substrate).first);
             out_mapping.cost += virtual1[e].bandwidth;
            curr = s;
            break;
          }
        }
      }
    }
  }
}
}
