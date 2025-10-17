#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <limits>
#include <boost/graph/adjacency_list.hpp>

using namespace boost;
using namespace std;


struct VertexResources {
  int x, y, cpu, memory;
  bool access;
  bool operator<(VertexResources& res) {
    if(cpu == res.cpu)
      return memory < res.memory;
    return cpu < res.cpu;
  }
  bool hasEnoughResources(VertexResources& vres) {
    return cpu >= vres.cpu && memory >= vres.memory;
  }
  void allocate(VertexResources& vres) {
    cpu -= vres.cpu;
    memory -= vres.memory;
  }
  void deallocate(VertexResources& vres) {
    cpu += vres.cpu;
    memory += vres.memory;
  }
};

typedef int Bandwidth;
struct EdgeResources {
  Bandwidth bandwidth;
  int id;
};

typedef adjacency_list<listS, vecS, undirectedS,
  VertexResources, EdgeResources> Graph;
typedef Graph::vertex_descriptor  Vertex;
typedef Graph::vertex_iterator    VertexIterator;
typedef Graph::out_edge_iterator  EdgeIterator;
typedef Graph::edge_descriptor    Edge;
enum VTypes { NONE = numeric_limits<Vertex>::max()};

struct Mapping {
  vector<Vertex> vMap_; // virtual vertex -> substrate vertex
  vector<map<Vertex,vector<Edge>>> eMap_; // virtual link -> substrate links
  vector<bool> substrateMapped;
  Bandwidth cost;
  
  explicit Mapping(Bandwidth _cost = 0) : cost(_cost) {}
  Mapping(const size_t numVerticesS, const size_t numVerticesV) : 
    vMap_(numVerticesV, NONE), eMap_(numVerticesV), 
    substrateMapped(numVerticesS, 0), cost(0)
  { }
  vector<Edge>& edge_map_ref(Vertex u, Vertex v) {
    if(u > v)
      swap(u,v);
    return eMap_[u][v];
  }
};

struct Parameters {
  double Ro;
  double Phi;
  double Tmin;
  double Tmax;
  int IterMax;
  int NAnts;
  int Hops;
  double D;
  bool onlyTheBestReinforces;
};

void read_file(const char* filename, Graph& g);
void print(Graph& g);
void print(Graph& s, Graph& v);
void print(Mapping& mapping);
int sumAccessCpu(Graph& g);
Bandwidth sumBandwidth(Graph& g);
Bandwidth maxBandwidth(Graph& g);
bool verify(Mapping& mapping, Graph& substrate, Graph& virtual1);
double truncateBetween(double x, double xmin, double xmax);

#define LOG(str) std::cout << "[" << (str) << "]"

#endif

