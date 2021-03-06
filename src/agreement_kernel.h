
#ifndef _AGREEMENT_KERNEL_HPP
#define _AGREEMENT_KERNEL_HPP

#include <set>

#include "hslinkern/hslinkern.hpp"
#include "MyTree.h" 

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
void construct_conflict_hypergraph(const MyTree& T1, const MyTree& T2, Hypergraph& G)
{
  assert(T1.triplets_set_up());
  assert(T2.triplets_set_up());
  const MatrixTriplets& tripletTree1 = T1.get_triplets(); 
  const MatrixTriplets& tripletTree2 = T2.get_triplets();

  const size_t dim = tripletTree1.size().first;
  for (unsigned i = 0; i < dim; i++){
    // translate preorder numbers of T1 into preorder numbers of T2
    const StId i_id = stid(T1.leaf_by_po_num(i));
    const unsigned i_in_T2 = T2[i_id]->getInfos().leaf_po_num();
  	for (unsigned j = i + 1; j < dim; j++){
      const StId j_id = stid(T1.leaf_by_po_num(j));
      const unsigned j_in_T2 = T2[j_id]->getInfos().leaf_po_num();
      
      // ATTENTION: i and j are leaf-postorder numbers whereas z is an stid
      const auto& trip1_ij = tripletTree1[{i,j}];
      const auto& trip2_ij = tripletTree2[{i_in_T2,j_in_T2}];
  		for(const StId z: trip1_ij)
  			if(trip2_ij.find(z) == trip2_ij.end())
  				G.emplace_back(vector<StId>{i_id, j_id, z});
  	}// for all j
  }	// for all i
}

// return the vertices of the HittingSet kernel of the conflict hypergraph of triplets in T1 & T2
set<StId>* disagreement_kernel(const MyTree& T1, const  MyTree& T2, const unsigned k)
{
  // step 1: build disagreement hypergraph
  Hypergraph G;
  construct_conflict_hypergraph(T1, T2, G);
  //cout << "conflict hypergraph contains "<<G.size()<<" edges" <<endl;

  // step 2: run kernelization
  return kernelize(G, k);
}

#endif
