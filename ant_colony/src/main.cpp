// Copyright 2013 Leonardo Moura
#include <random>
#include <numeric>
#include <iostream>
#include <list>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "util.h"
#include "vmap.h"

using namespace std;
using namespace boost;

double maxBand;

double averageUsage(Graph& g) {
  double sumBand = sumBandwidth(g);
  //double maxBand = maxBandwidth(g);
  return sumBand / (static_cast<double>(num_edges(g)) * maxBand );
}
double averageUsage(Graph& g, double totalBand) {
  return totalBand / (static_cast<double>(num_edges(g)) * maxBandwidth(g) );
}

int main(const int argc, char** argv) {
  if(argc != 2) {
    cout << "wrong number of arguments" << endl;    
    cout << "usage: " << endl;
    cout << "vne.out parameter_file.ini" << endl;
  }
  property_tree::ptree pt;
  property_tree::ini_parser::read_ini(argv[1], pt);
  Parameters p { pt.get<double>("ro"), pt.get<double>("phi"),
    pt.get<double>("tmin"), pt.get<double>("tmax"),pt.get<int>("maxiter"),
    pt.get<int>("ants"), pt.get<int>("hops"), pt.get<double>("d"),
    pt.get<bool>("onlythebestreinforces")};

  auto seed = time(NULL);
  srand48(seed);
  srand(seed);
  if(!pt.get<bool>("simulation.random_seed")) {
    seed = pt.get<int>("simulation.seed");
  }
  std::minstd_rand0 generator(seed);

  string virtualFilePrefix = pt.get<string>("simulation.virtual_file_prefix");
  Graph substrate;
  read_file(pt.get<string>("simulation.substrate_file").c_str(), substrate);

  maxBand = maxBandwidth(substrate);
  
  exponential_distribution<double> 
    lifetimeDist(1 / pt.get<double>("simulation.expectedlifetime"));
  poisson_distribution<int>
    requestDist(pt.get<int>("simulation.requestsperblock"));
  const int N_REQ = pt.get<int>("simulation.nreq");
  const int blocksize = pt.get<int>("simulation.timeunitsperblock");

  struct Request {
    int id;
    int expiration;
    Graph network;
    Mapping mapping;
    bool accepted;
    bool operator<(Request& req) const {
      return expiration < req.expiration;
    }
  };
  //Bandwidth before = sumBandwidth(substrate);  
  int time {0}, totalRequests {0}, rejects {0};
  int totalCpu {0}, ocCpu {0}, sCpu{0}, impossible {0};
  Bandwidth totalCost {0};
  Bandwidth totalBand {0};
  list<Request> activeRequests;
  LOG("MAIN") << "total cpu: " << (sCpu = sumAccessCpu(substrate)) << endl;
  while(true) {
    // free expired requests
    LOG("MAIN") << time << endl;
    while(!activeRequests.empty() &&
        activeRequests.front().expiration <= time) {
      //int ssum = sumAccessCpu(substrate);
      //int sum = sumAccessCpu(activeRequests.front().network);
      int netcpu = sumAccessCpu(activeRequests.front().network);
      totalCpu -= netcpu;
      if(!activeRequests.front().accepted) {
        activeRequests.pop_front();
        continue;
      }
      ocCpu -= netcpu;
      totalBand -= activeRequests.front().mapping.cost;
      VNEAC::free_vmap(substrate, activeRequests.front().network, 
        activeRequests.front().mapping);
      //assert(ssum + sum == sumAccessCpu(substrate));
      std::cout << "[MAIN]" << " free " << activeRequests.front().id << endl;
      activeRequests.pop_front();
    }

    // generate requests for this block and map them
    int requests = requestDist(generator);
    for(int i = 0; i < requests && totalRequests < N_REQ; ++i){
      ++totalRequests;
      Request r;
      r.id = totalRequests;
      r.expiration = time + lifetimeDist(generator);
      LOG("MAIN") << "request " << r.id << " expires at " << r.expiration << endl;

      stringstream filename;
      filename << virtualFilePrefix << totalRequests << ".gtt";
      read_file(filename.str().c_str(), r.network);
      LOG("MAIN") << "map " <<  filename.str() << endl;
      impossible += (ocCpu > sCpu);
      int netcpu = sumAccessCpu(r.network);
      totalCpu += netcpu; 
      if(VNEAC::vmap(substrate, r.network, p, r.mapping)) {
        ocCpu += netcpu;
        totalBand += r.mapping.cost;
        totalCost += r.mapping.cost * (r.expiration - time);
        r.accepted = true;
        activeRequests.push_back(r);
        LOG("MAIN") << "accepted" << endl;
        //assert(old[r.mapping.vMap_[0]].cpu > substrate[r.mapping.vMap_[0]].cpu); 
        //LOG("DIFF") << old[r.mapping.vMap_[0]].cpu << " > " << substrate[r.mapping.vMap_[0]].cpu << endl;
      
      } else {
        r.accepted = false;
        activeRequests.push_back(r);
        //assert(bbefore == sumBandwidth(substrate));
        LOG("MAIN") << "rejected" << endl;
        rejects++;
      }
    }
    activeRequests.sort();
    LOG("GRAPH") << "rejrate:" << static_cast<double>(rejects)/totalRequests 
      << ":" << time << endl;
    LOG("GRAPH") << "cpu:" << totalCpu << ":" << time << endl;
    LOG("GRAPH") << "occcpu:" << ocCpu << ":" << time << endl;
    LOG("GRAPH") << "usage:" << averageUsage(substrate) << ":" << time << endl;
    LOG("GRAPH") << "cost:" << totalCost << ":" << time << endl;
    if(activeRequests.empty() && totalRequests >= N_REQ)
      break;
    time += blocksize;

    //LOG("sum")<<"v: "<< accumulate(activeRequests.begin(), activeRequests.end(), 0,
    //    [](int sum, Request& r) { return sum + sumAccessCpu(r.network); }) << endl;
    //LOG("sum") << sumAccessCpu(substrate) << endl;
  }
  //assert(before == sumBandwidth(substrate));
  LOG("MAIN") << "impossible:" << impossible << endl;
  cout << "[MAIN] " << "rejectrate:" << static_cast<double>(rejects)/N_REQ << endl;
  cout << "[MAIN] " << "finalcost:" << totalCost << endl;
  return 0;
}
