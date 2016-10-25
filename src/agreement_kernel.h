
#ifndef _AGREEMENT_KERNEL_HPP
#define _AGREEMENT_KERNEL_HPP

#include <set>

#include "hslinkern/hslinkern.hpp"
#include "MyTree.h" 

// wrapper for the HittingSet kernel
set<StId> kernelize(const Hypergraph& G, const unsigned k)
{
  // step 1: run kernel
  Graphstats stats(get_stats(G));
  Hypergraph H(kernelize(G, stats, k));

  // step 2: collect remaining vertices
  set<StId> result;
  for(const auto& edge: H)
    for(const StId& v: edge)
      result.insert(v);

  return result;
}


// the conflict hypergraph is the graph containing all triplets abc of IDs of T1 on which T1 and T2 disagree
/** NOTE: this assumes that all clades have been set up */
void construct_conflict_hypergraph(const MyTree& T1, const MyTree& T2, Hypergraph& G)
{
  #warning improve the handeling of already computed triplet sets
  assert(T1.node_infos_set_up());
  assert(T2.node_infos_set_up());
  const MatrixTriplets* const tripletTree1 = T1.get_triplets(); 
  const MatrixTriplets* const tripletTree2 = T2.get_triplets();

  const size_t dim = tripletTree1->size().first;
  for (unsigned i = 0; i < dim; i++){
    // translate preorder numbers of T1 into preorder numbers of T2
    const StId i_id = stid(T1.leaf_by_po_num(i));
    const unsigned i_in_T2 = T2[i_id]->getInfos().leaf_po_num();
  	for (unsigned j = i + 1; j < dim; j++){
      const StId j_id = stid(T1.leaf_by_po_num(j));
      const unsigned j_in_T2 = T2[j_id]->getInfos().leaf_po_num();
  		for (const unsigned z: tripletTree1->at({i, j}))
  			if(is_conflict(*tripletTree2, i_in_T2, j_in_T2, z))
  				G.emplace_back(vector<StId>{i_id, j_id, z}); //stIds
  	}// for all j
  }	// for all i
  delete tripletTree1;
  delete tripletTree2;
}

// return the vertices of the HittingSet kernel of the conflict hypergraph of triplets in T1 & T2
const set<StId> disagreement_kernel(const MyTree& T1, const  MyTree& T2, const unsigned k)
{
  // step 1: build disagreement hypergraph
  Hypergraph G;
  construct_conflict_hypergraph(T1, T2, G);
  //cout << "conflict hypergraph contains "<<G.size()<<" edges" <<endl;

  // step 2: run kernelization
  const set<StId> kernel = kernelize(G, k);

  //cout << "kernel has size "<<kernel.size()<<endl;
  // step 3: return result
  return kernel;
}

#endif
