#pragma once
#include <vector>
#include <queue>

#include "../include/util.h"
#include "../include/vnedefs.h"

namespace vne{
bool bprice_rec(Graph& substrate, const Graph& virtual1, const Parameters param, Mapping& outMapping);

enum BPStrategy { DFS, BFS, BEST_PROJECTION };

struct BPNode {
  size_t depth;
  int father;

  size_t  c1;
  size_t  c2;
  CplexValue dual;
  float infeasibilities;

  BPNode(size_t _depth, int _father,
      size_t _c1, size_t _c2, CplexValue _dual, float _infeasibilities);

  virtual void apply(CGModel& model) const = 0;
  virtual void undo(CGModel& model) const = 0;
};


typedef std::priority_queue<size_t, std::vector<size_t>,
    std::function<bool(size_t, size_t)>> QueueType;

struct BPTree {
  std::vector<BPNode*> nodes;
  QueueType Q;

  void push(BPNode* node);
  size_t pop();
};

void branch_on_nodes(const CGModel&, const Vertex sel_v, const Vertex sel_s,
    const size_t current_i, BPTree&);

void branch_on_nodes_dichotomic(const CGModel& model, const Vertex sel_v,
    const Vertex sel_s, const size_t current_i, BPTree&); 

typedef void (*BranchOnNodesFun)(const CGModel&, const Vertex, const Vertex,
    const size_t, BPTree&); 

template <BPStrategy Strategy,
    BranchOnNodesFun branch_on_nodes_fun = &branch_on_nodes>
bool bprice(Graph& substrate, const Graph& virtual1,
    const Parameters param, Mapping& out_mapping);
}

