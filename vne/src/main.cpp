// Copyright 2013 Leonardo Moura
#include <cassert>
#include <cstdio>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <random>
#include <list>
#include <numeric>
#include <string>
#include "../include/util.h"
#include "../include/backtrack.h"
#include "../include/cplex.h"
//#include "../include/bb.h"
#include "../include/cg.h"
#include "../include/fix.h"
//#include "../include/bbcgu.h"
#include "../include/bprice.h"
#include "../include/ac.h"

using std::cout;
using std::endl;
using std::vector;

typedef bool (*VneMethod)(Graph&, const Graph&, const Parameters, Mapping&);

int main(const int argc, char** argv) {
  if (argc != 3) {
    printf("error: correct usage: vne.out method params.ini       \n");
    printf("method -> 0 = backtrack, 1 = cplex, 2 = cplex root    \n");
    printf("       -> 10 = cg, 11 = bp rec, 12 = bp dfs, 13 = bfs \n");
    return 1;
  }
  const int method_code = atoi(argv[1]);
  char*& params_file = argv[2];

  vector<VneMethod> methods(40);
  methods[0] = &vne::backtrack<vne::Backtrack>;
  methods[1] = &vne::cplex;
  methods[2] = &vne::cplexrelax;
  methods[3] = &vne::backtrack<vne::HopBacktrack>;
  methods[4] = &vne::backtrack<vne::GreedyInitialBacktrack>;
  methods[5] = &vne::backtrack<vne::HopBacktrack, vne::Sort<vne::MoreCpu>>;
  methods[6] = &vne::backtrack<vne::HopBacktrack, vne::Sort<vne::GreaterDegree>>;
  //methods[7] = &vne::cg;
  methods[8] = &vne::fix;
  //methods[9] = &vne::bb;
  methods[10] = &vne::cgunsplittable;
  methods[11] = &vne::bprice_rec;
  methods[12] = &vne::bprice<vne::BPStrategy::DFS>;
  methods[13] = &vne::bprice<vne::BPStrategy::BFS>;
  methods[14] = &vne::bprice<vne::BPStrategy::BEST_PROJECTION>;
  methods[15] = &vne::bprice<vne::BPStrategy::DFS, &vne::branch_on_nodes_dichotomic>;
  methods[16] = &vne::bprice<vne::BPStrategy::BFS, &vne::branch_on_nodes_dichotomic>;
  methods[17] = &vne::bprice<vne::BPStrategy::BEST_PROJECTION, &vne::branch_on_nodes_dichotomic>;
  // branch and bound guided by column generation
  //methods[13] = &vne::bbcgu;
  methods[18] = &vne::ac;

  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(params_file, pt);

  srand48(time(NULL));
  srand(time(NULL));

  Parameters param {pt.get<int>("timelimit_in_s"),
    pt.get_optional<float>("ant_colony.ro").get_value_or(0.0),
    pt.get_optional<float>("ant_colony.phi").get_value_or(0),
    pt.get_optional<float>("ant_colony.tmin").get_value_or(0),
    pt.get_optional<float>("ant_colony.tmax").get_value_or(0),
    pt.get_optional<int>("ant_colony.maxiter").get_value_or(0),
    pt.get_optional<int>("ant_colony.ants").get_value_or(0),
    pt.get_optional<int>("ant_colony.hops").get_value_or(0),
    pt.get_optional<int>("ant_colony.d").get_value_or(0),
    pt.get_optional<bool>("ant_colony.onlythebestreinforces").get_value_or(0)};


  Graph substrate, virtual1;
  read_file(pt.get<std::string>("substrate_file").c_str(), &substrate);
  read_file(pt.get<std::string>("virtual_file").c_str(), &virtual1);

  Mapping mapping;
  if ((*methods[method_code])(substrate, virtual1, param, mapping)) {
    print(mapping);
    assert(verify(mapping, substrate, virtual1));
  } 

  return 0;
}
