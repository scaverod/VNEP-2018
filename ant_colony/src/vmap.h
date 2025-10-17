#ifndef VMAP_H_
#define VMAP_H_

#include "util.h"

namespace VNEAC {
  bool vmap(Graph& substrate, Graph& virtual1, const Parameters& param, Mapping& outMapping);
  void free_vmap(Graph& substrate, Graph& virtual1, Mapping& mapping);
}

#endif
