#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include "util.h"

using namespace std;

void read_file(const char* filename, Graph& g) {
    //std::locale l = std::cout.getloc();
    //l.name();
    std::fstream filestr (filename, std::fstream::in);
    char chrline[255];
    while(filestr.good()) {
      filestr.getline(chrline,255);
      stringstream sstream(chrline);
      switch(chrline[0]){
        case 'G':
          break;
        case 'V': {
            VertexResources vres;
            int n;
            char dummy;
            sstream >> dummy >> n >> vres.x >> vres.y >> vres.cpu >> vres.memory 
              >> vres.access;
            add_vertex(vres, g);
            assert(num_vertices(g) > 0);
                  }
          break;
        case 'E': {
          int u, v;
          char dummy;
          EdgeResources eres;
          eres.id = num_edges(g);
          sstream >> dummy >> u >> v >> eres.bandwidth;
          add_edge(u, v, eres, g);
          break;
        }
      }
    }
    filestr.close();
}

void print(Graph& g, char suffix) {
  VertexIterator v, vend;
  for(tie(v,vend) = vertices(g); v != vend; ++v) {
    cout << suffix << *v << " [ pos=\"" << g[*v].x << "," << g[*v].y << "!\"]" << endl; 
    EdgeIterator e, eend;
    for(tie(e, eend) = out_edges(*v, g); e != eend; ++e) {
      if(target(*e, g) > *v)
        continue;
      cout << suffix << *v << " -- " << suffix << target(*e, g) << ";" << endl;
    }
  }
}

void print(Graph& s, Graph& v) {
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

void print(Mapping& mapping) {
  for(size_t v = 0; v != mapping.vMap_.size(); ++v) {
    if(mapping.vMap_[v] == NONE)
      continue;
    cout << "v" << v << " -- s" << mapping.vMap_[v] << " ;" << endl;
   for(auto& w: mapping.eMap_[v]) {
    cout << "v " << w.first << ": ";
    for(Edge e: w.second)
      cout << e << " ";
    cout << endl;
   }
  }
}

int sumAccessCpu(Graph& g) {
 return accumulate(vertices(g).first, 
        vertices(g).second, 0, 
        [&g](int sum, Vertex v) { 
          return sum + (g[v].access ? g[v].cpu : 0); });
}

Bandwidth sumBandwidth(Graph& g) {
  Bandwidth sum {0};
  for(Vertex v = 0; v != num_vertices(g); ++v){
    EdgeIterator e, eend;
    for(tie(e,eend) = out_edges(v, g); e != eend; ++e){
      if(source(*e, g) > target(*e, g))
        sum += g[*e].bandwidth;
    }
  }
  return sum;
}

Bandwidth maxBandwidth(Graph& g) {
  Bandwidth maxb {0};
  for(Vertex v = 0; v != num_vertices(g); ++v){
    EdgeIterator e, eend;
    for(tie(e,eend) = out_edges(v, g); e != eend; ++e){
        maxb = max(maxb, g[*e].bandwidth);
    }
  }
  return maxb;
}

bool verify(Mapping& mapping, Graph& substrate, Graph& virtual1) {
  set<Vertex> substrateVertices;
  Bandwidth cost {0};
  LOG("VERIFY") << num_vertices(virtual1) << endl;
  for(size_t v = 0; v < num_vertices(virtual1); ++v) {
    LOG("VERIFY") << v << endl;
    Vertex s {mapping.vMap_[v]};
    substrateVertices.insert(s);
    // verify if substrate vertex has enough resources
    assert(substrate[s].memory >= 0 &&
      substrate[s].cpu >= 0);
    assert(mapping.vMap_[v] != NONE);
    // verify if all edges are mapped and all paths can hold the bandwidth
    EdgeIterator e, eend;
    tie(e, eend) = out_edges(v, virtual1);
    for(; e != eend; ++e) {
      Vertex vt = target(*e, virtual1);
      if(v > vt)
        continue;
      Vertex st = mapping.vMap_[vt];
      vector<Edge>& path = mapping.edge_map_ref(v, vt);
      assert(path.size() > 0);
      // verify if the path is correct
      assert(boost::target(path.front(),substrate) == s 
          || boost::target(path.front(), substrate) == st);
      assert(boost::source(path.back(),substrate) == s 
          || boost::source(path.back(),substrate) == st);
      for(const Edge& ep: path) {
        assert(substrate[ep].bandwidth >= 0);
        cost += virtual1[*e].bandwidth;
      }
    }
  }
  // verify if each virtual vertex was mapped to a different substrate vertex
  assert(cost == mapping.cost);
  assert(substrateVertices.size() == num_vertices(virtual1));
  return true;
}

double truncateBetween(double x, double xmin, double xmax) {
  return max(min(x, xmax), xmin);
}
