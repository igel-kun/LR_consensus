
#pragma once

#include <unordered_set>

#include "vector2d.h"
//#include "array2d.h"

#include "NodeInfos.h"
#include "MyTree.h"

using namespace bpp;

//! a 2d vector assigning a pair of PREORDER NUMBERS to a set of STIDs such that the three leaves form a conflict
//typedef symmetric_array2d<unordered_set<unsigned> > MatrixTriplets;
struct MatrixTriplets: public symmetric_vector2d<unordered_set<unsigned> >
{
  typedef symmetric_vector2d<unordered_set<unsigned> > Parent;
  using Parent::symmetric_vector2d;

  void add_triple(const unsigned x, const unsigned y, const StId z)
  {
    assert((x <= z) && (y <= z));
    Parent::operator[]({x, y}).insert(z);
  }

  bool is_triple(const unsigned x, const unsigned y, const StId z) const
  {
    assert((x <= z) && (y <= z));
    const auto& conflict_set = Parent::operator[]({x ,y});
    return conflict_set.find(z) == conflict_set.end();
  }

  void add_conflicts_at_LCA_of(const StId x, const StId y, const MyTree& T1, const MyTree& T2, const unsigned n)
  {
    const MyNode* const T1_xy_lca = T1.getLCA(T1[x], T1[y]);
    const MyNode* const T2_xy_lca = T2.getLCA(T2[x], T2[y]);
    for(StId z = max(x,y) + 1; z < n; ++z)
      if(T1.displays_triple(T1_xy_lca, z) != T2.displays_triple(T2_xy_lca, z))
        add_triple(x, y, z);
  }

  //! fill the triplet matrix with conflict triples of T1 & T2
  void fill_conflicts(const MyTree& T1, const MyTree& T2)
  {
    const unsigned n = T1.num_leaves();
    assert(n == T2.num_leaves());
    for(StId x = 0; x < n; ++x)
      for(StId y = x + 1; y < n; ++y)
        add_conflicts_at_LCA_of(x, y, T1, T2, n);
  }

  void remove_conflicts_involving(const StId x)
  {
    const unsigned n = size().first;
    // remove (x,y,z) and (y,x,z)
    for(StId y = 0; y < n; ++y) operator[]({x,y}).clear();
    // remove (y,z,x)
    for(StId y = 0; y < n; ++y)
      for(StId z = y + 1; z < n; ++z)
        operator[]({y,z}).erase(x);
  }


  void add_conflicts_involving(const StId x, const MyTree& T1, const MyTree& T2)
  {
    const unsigned n = size().first;
    // add (x,y,z) and (y,x,z)
    for(StId y = 0; y < n; ++y)
      add_conflicts_at_LCA_of(x, y, T1, T2, n);
    // add (y,z,x)
    for(StId y = 0; y < x; ++y)
      for(StId z = 0; z < x; ++z)
        if(T1.displays_triple(y, z, x) != T2.displays_triple(y, z, x))
          add_triple(y, z, x);
  }

};



