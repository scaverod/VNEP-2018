#pragma once
#include <iostream>

#include <ilcplex/ilocplex.h>
#include "../include/util.h"
#include "../include/vnedefs.h"

namespace vne{
bool cplex(Graph& substrate, const Graph& virtual1, const Parameters param, Mapping& outMapping);

bool cplex_improved(Graph& substrate, const Graph& virtual1, const Parameters param, Mapping& outMapping);

bool cplexrelax(Graph& substrate, const Graph& virtual1, const Parameters param, Mapping& outMapping);
}
