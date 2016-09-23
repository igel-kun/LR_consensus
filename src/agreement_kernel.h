
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

typedef pair<const MyNode*, const MyNode*> LeafPair;
typedef vector<LeafPair> LeafCorrespondance;
typedef map<int, LeafPair> IdToLeafCorrespondance;

// return whether the 3 given leaf-correspondances form a conflict triple between T1 and T2
bool is_conflict(const MyTree& T1,
                 const MyTree& T2,
                 const LeafPair& leaf1,
                 const LeafPair& leaf2,
                 const LeafPair& leaf3)
{
#warning TODO
  return true;
}

// compute a correspondance between leaves of T1 and T2
void construct_leaf_correspondance(const MyTree& T1, const MyTree& T2, IdToLeafCorrespondance& result)
{
  // get the leaves of the trees
  const vector<const MyNode*> T1_leaves = T1.getLeaves();
  const vector<const MyNode*> T2_leaves = T2.getLeaves();
  
  // construct a mapping of names to nodes of T2
  map<string, const MyNode*> name_to_node;
  for(const auto& leaf: T2_leaves)
    name_to_node.emplace(leaf->getName(), leaf);

  // construct the resulting leaf correspondance
  for(const auto& leaf: T1_leaves){
    const auto T2_leaf = name_to_node.find(leaf->getName());
    if(T2_leaf != name_to_node.cend())
      result.emplace(leaf->getId(), LeafPair(leaf, T2_leaf->second));
  }
}

// the conflict hypergraph is the graph containing all triplets abc of IDs of T1 on which T1 and T2 disagree
void construct_conflict_hypergraph(const MyTree& T1, const MyTree& T2, const IdToLeafCorrespondance& id_to_leaves, Hypergraph& G)
{
  // consider all triplets to build the hypergraph
  for(auto leaf1 = id_to_leaves.begin(); leaf1 != id_to_leaves.end(); ++leaf1)
    for(auto leaf2 = next(leaf1); leaf2 != id_to_leaves.end(); ++leaf2)
      for(auto leaf3 = next(leaf2); leaf3 != id_to_leaves.end(); ++leaf3)
        if(is_conflict(T1, T2, leaf1->second, leaf2->second, leaf3->second))
          G.emplace_back(vector<int>{leaf1->first, leaf2->first, leaf3->first});
}

// return the vertices of the HittingSet kernel of the conflict hypergraph of triplets in T1 & T2
vector<LeafPair> disagreement_kernel(const MyTree& T1, const MyTree& T2, const unsigned k)
{
  // step 1: build a leaf-correspondance from T1 to T2
  IdToLeafCorrespondance corr;
  construct_leaf_correspondance(T1, T2, corr);
  // step 1: build disagreement hypergraph
  Hypergraph G;
  construct_conflict_hypergraph(T1, T2, corr, G);

  // step 2: run kernelization
  const set<int> kernel = kernelize(G, k);

  // step 3: construct result
  vector<LeafPair> result;
  result.reserve(kernel.size());
  for(const int i: kernel)
    result.emplace_back(corr.at(i));
  return result;
}

#endif
