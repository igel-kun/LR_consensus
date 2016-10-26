
#pragma once

#include <unordered_set>

#include "vector2d.h"
//#include "array2d.h"

#include "NodeInfos.h"

using namespace bpp;

//! a 2d vector assigning a pair of PREORDER NUMBERS to a set of STIDs such that the three leaves form a conflict
//typedef symmetric_array2d<unordered_set<unsigned> > MatrixTriplets;
typedef symmetric_vector2d<unordered_set<unsigned> > MatrixTriplets;

void add_triple(MatrixTriplets& trip, const unsigned x, const unsigned y, const unsigned z)
{
  trip[{x, y}].insert(z);
}

bool is_not_triple(const MatrixTriplets& trip, const unsigned x, const unsigned y, const unsigned z)
{
  const auto& conflict_set = trip[{x ,y}];
  return conflict_set.find(z) == conflict_set.end();
}

