
#ifndef _AGREEMENT_KERNEL_HPP
#define _AGREEMENT_KERNEL_HPP

#include <set>

#include "hslinkern/hslinkern.hpp"
#include "MyTree.h" 
#include "MatrixTriplets.h"

// wrapper for the HittingSet kernel
set<StId>* kernelize(const Hypergraph& G, const unsigned k)
{
  // step 1: run kernel
  Graphstats stats(get_stats(G));
  Hypergraph* H(kernelize(G, stats, k));

  if(H){
    // step 2: collect remaining vertices
    set<StId>* result = new set<StId>();
    for(const auto& edge: *H)
      for(const StId& v: edge)
        result->insert(v);
    return result;
  } else return NULL;
}


// the conflict hypergraph is the graph containing all triplets abc of IDs of T1 on which T1 and T2 disagree
/** NOTE: this assumes that all clades have been set up */
void construct_conflict_hypergraph(const MatrixTriplets& conflicts, Hypergraph& G)
{
  const size_t dim = conflicts.size().first;
  for(StId x = 0; x < dim; ++x)
    for(StId y = x + 1; y < dim; ++y)
      for(StId z: conflicts[{x, y}])
        G.emplace_back(vector<StId>{x, y, z});
}

// return the vertices of the HittingSet kernel of the conflict hypergraph of triplets in T1 & T2
set<StId>* disagreement_kernel(const MatrixTriplets& conflicts, const unsigned k)
{
  // step 1: build disagreement hypergraph
  Hypergraph G;
  construct_conflict_hypergraph(conflicts, G);
  //cout << "conflict hypergraph contains "<<G.size()<<" edges" <<endl;

  // step 2: run kernelization
  return kernelize(G, k);
}

#endif
