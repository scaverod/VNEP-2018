// Copyright 2013 Leonardo Moura
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include "../include/util.h"
#include "../include/vnedefs.h"

using std::cout;
using std::endl;
using std::stringstream;
using boost::tuples::tie;
Mapping::Mapping(const Mapping& mapping) : vertex(mapping.vertex),
    substrate_mapped(mapping.substrate_mapped), cost(mapping.cost),
    edge(mapping.edge.size()) {
  for (Vertex v = 0; v < mapping.edge.size(); ++v) {
    edge[v] = std::vector<std::vector<Edge>*>(v+1, nullptr);
    for (Vertex w = 0; w < mapping.edge[v].size(); ++w) {
      if (mapping.edge[v][w] != nullptr) {
        edge[v][w] = new std::vector<Edge>(mapping.edge[v][w]->begin(),
            mapping.edge[v][w]->end());
      }
    }
  }
}

void read_file(const char* filename, Graph* g) {
  const int LINE_SIZE = 255;
  std::fstream filestr(filename, std::fstream::in);
  char chrline[LINE_SIZE];
  size_t edges = 0;
  while (filestr.good()) {
    filestr.getline(chrline, LINE_SIZE);
    stringstream sstream(chrline);
    switch (chrline[0]) {
      case 'G':
        break;
      case 'V': {
          VertexResources vres;
          int n;
          char dummy;
          sstream >> dummy >> n >> vres.cpu;
          add_vertex(vres, *g);
          assert(num_vertices(*g) > 0);
                }
        break;
      case 'E': {
        int u, v;
        char dummy;
        Bandwidth bandwidth;
        sstream >> dummy >> u >> v >> bandwidth;
        EdgeResources eres {bandwidth, edges++};
        add_edge(u, v, eres, *g);
        break;
      }
    }
  }
  filestr.close();
}

void print(const Graph& g, const char prefix) {
  VertexIterator v, vend;
  for (tie(v, vend) = vertices(g); v != vend; ++v) {
    cout << prefix << *v << endl;
    EdgeIterator e, eend;
    for (tie(e, eend) = out_edges(*v, g); e != eend; ++e) {
      if (target(*e, g) > *v)
        continue;
      cout << prefix << *v << " -- " << prefix << target(*e, g) << " (" <<
        g[*e].bandwidth << ");" << endl;
    }
  }
}

void print(const Graph& s, const Graph& v) {
  cout << "graph G { " << endl;
  cout << " subgraph substrate { " << endl;
  cout << " style = \"filled\";" << endl;
  print(s, 's');
  cout << "label = \"substrate\";" << endl;
  cout << "} " << endl;
  cout << " subgraph virtual { " << endl;
  cout << " style = \"filled\";" << endl;
  print(v, 'v');
  cout << "label = \"virtual\";" << endl;
  cout << "   } ";
  cout << "}";
}

void print(const Mapping& mapping) {
  for (Vertex v = 0; v != mapping.vertex.size(); ++v) {
    if (mapping.vertex[v] == NONE)
      continue;
    cout << "v" << v << " -- s" << mapping.vertex[v] << " ;" << endl;
    for (Vertex w = 0; w < v; ++w) {
      std::vector<Edge>* const path = mapping.edge_map_ref(v, w);
      if (path == nullptr)
        continue;
      cout << "v " << w << "= ";
      for (const Edge e : *path)
        cout << e << " ";
      cout << endl;
    }
  }
}

int sumCpu(const Graph& g) {
  return std::accumulate(vertices(g).first,
                         vertices(g).second, 0,
                         [&g](int sum, Vertex v) {
                             return sum + g[v].cpu; });
}

Bandwidth sum_bandwidth(const Graph& g) {
  Bandwidth sum {0};
  for (Vertex v = 0; v != num_vertices(g); ++v) {
    EdgeIterator e, eend;
    for (tie(e, eend) = out_edges(v, g); e != eend; ++e) {
      if (source(*e, g) > target(*e, g))
        sum += g[*e].bandwidth;
    }
  }
  return sum;
}

Bandwidth out_bandwidth(const Graph& g, const Vertex v) {
  Bandwidth sum {0};
  EdgeIterator e, eend;
  for (tie(e, eend) = out_edges(v, g); e != eend; ++e)
    sum += g[*e].bandwidth;
  return sum;
}

Bandwidth maxBandwidth(const Graph& g) {
  Bandwidth maxb {0};
  for (Vertex v = 0; v != num_vertices(g); ++v) {
    EdgeIterator e, eend;
    for (tie(e, eend) = out_edges(v, g); e != eend; ++e) {
        maxb = std::max(maxb, g[*e].bandwidth);
    }
  }
  return maxb;
}

bool verify(const Mapping& mapping, Graph substrate,
            const Graph& virtual1) {
  cout << "verifying" << endl;
  std::set<Vertex> substrateVertices;
  Bandwidth cost {0};
  for (Vertex v = 0; v < num_vertices(virtual1); ++v) {
    Vertex s {mapping.vertex[v]};
    substrateVertices.insert(s);
    // verify if substrate vertex has enough resources
    assert(substrate[s].cpu >= virtual1[v].cpu);
    assert(mapping.vertex[v] != NONE);
    // verify if all edges are mapped and all paths can hold the bandwidth
    EdgeIterator e, eend;
    tie(e, eend) = out_edges(v, virtual1);
      for (; e != eend; ++e) {
      Vertex vt = target(*e, virtual1);
      if (v < vt)
        continue;
      const Vertex st = mapping.vertex[vt];
      std::vector<Edge>* const path = mapping.edge_map_ref(v, vt);
      assert(path->size() > 0);
      // verify if the path is correct
      assert(boost::target(path->front(), substrate) == s ||
          boost::target(path->front(), substrate) == st);
      assert(boost::source(path->back(), substrate) == s  ||
          boost::source(path->back(), substrate) == st);
      cout << "check " << v << ", " << vt << endl;
      for (const Edge& ep : *path) {
        cout << source (ep, substrate) << ", " << target(ep, substrate) << " - " << substrate[ep].bandwidth << endl;
        cost += virtual1[*e].bandwidth;
        substrate[ep].bandwidth -= virtual1[*e].bandwidth;
        cout << "b = " << substrate[ep].bandwidth << endl;
        assert(substrate[ep].bandwidth >= 0);
      }
    }
  }
  // verify if each virtual vertex was mapped to a different substrate vertex
  assert(substrateVertices.size() == num_vertices(virtual1));
  assert(cost == mapping.cost);
  return true;
}
