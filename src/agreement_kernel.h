
#ifndef _AGREEMENT_KERNEL_HPP
#define _AGREEMENT_KERNEL_HPP

#include <set>

#include "hslinkern/hslinkern.hpp"
#include "MyTree.h" 

// wrapper for the HittingSet kernel
set<int> kernelize(const Hypergraph& G, const unsigned k)
{
  // step 1: run kernel
  Graphstats stats(get_stats(G));
  Hypergraph H(kernelize(G, stats, k));

  // step 2: collect remaining vertices
  set<int> result;
  for(const auto& edge: H)
    for(const int& v: edge)
      result.insert(v);

  return result;
}


// the conflict hypergraph is the graph containing all triplets abc of IDs of T1 on which T1 and T2 disagree
void construct_conflict_hypergraph( MyTree& T1,  MyTree& T2, Hypergraph& G)
{
  if(! T1.isPreprocessed()){
  	T1.setClades();
  }	
  #warning improve the handeling of already computed triplet sets
  MatrixTriplets * tripletTree1 = T1.setTriplets(); 
  
  if(! T2.isPreprocessed()){
  	T2.setClades();
  }
  
  MatrixTriplets * tripletTree2 = T2.setTriplets();

  int dim = tripletTree1->getDim();
  for (int i=0;i< dim;i++){
	for (int j=0;j< dim;j++){
		for (int z:tripletTree1->getSet(i,j)){
			if(! tripletTree2->isTriplet(i,j,z))
				G.emplace_back(vector<int>{i,j,z}); //stIds
		}
	}
  }	
  tripletTree1->deleteMatrix();
  tripletTree2->deleteMatrix();
  delete tripletTree1;
  delete tripletTree2;

	
}

// return the vertices of the HittingSet kernel of the conflict hypergraph of triplets in T1 & T2
const set<int> disagreement_kernel( MyTree& T1,  MyTree& T2, const unsigned k)
{
  // step 1: build disagreement hypergraph
  Hypergraph G;
  construct_conflict_hypergraph(T1, T2, G);

  // step 2: run kernelization
  const set<int> kernel = kernelize(G, k);
  
  // step 3: return result

  return kernel;
}

#endif
