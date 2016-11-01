

#pragma once

#include <list>
#include "MyTree.h"
#include "mastDP.h"
#include "MatrixTriplets.h"
#include "agreement_kernel.h"

class CandidateTree;



// a regraft operation is a pair of vertices: the sibling of the origin and the sibling of the target
// Thus, to apply the operation, regraft the sibling of the first on top of the second
// To undo the operation, regraft the sibling of the second on top of the first
typedef StId RegraftOp;

class RegraftCandidates: public set<StId>
{
public:
  const StId sibling; //< the sibling of the leaf which is being regrafted

  RegraftCandidates(const StId _sibling):
    set<StId>::set(),
    sibling(_sibling)
  {}

};


//! a relation between the candidate tree and another tree
class TreeRelation
{
public:
  const CandidateTree& candidate;
  const MyTree& tree;
protected:
  MatrixTriplets conflicts;
  TracableDPTable MAST_table;
public:


  TreeRelation(const CandidateTree& _candidate, const MyTree& _tree);

  unsigned get_mast(list<StId>* const leaf_list = NULL) const 
  {
    if(leaf_list) MAST_table.trace(*leaf_list);
    return MAST_table.get_root_score();
  }

  unsigned get_mast_dist(list<StId>* const leaf_list = NULL) const;

  const set<StId>* disagreement_kernel(const unsigned k) const
  {
    Hypergraph G;
    construct_conflict_hypergraph(conflicts, G);
    return kernelize(G, k);
  }

  //! react to a regraft in the candidate tree (the regraft has been done before)
  void react_to_regraft(const MyNode* leaf, const MyNode* const old_sibling)
  {
    // Step 1: update conflicts
    conflicts.remove_conflicts_involving(stid(leaf));
    conflicts.add_conflicts_involving(stid(leaf), (const MyTree&)candidate, tree);
    // Step 2: update MAST table
    MAST_table.react_to_regraft(leaf, old_sibling);
  }
};

//! a candidate tree is a special kind of tree that carries tables of conflicting triples and MAST DP tables
class CandidateTree: public MyTree
{
private:
  using MyTree::MyTree;
  using MyTree::StId_to_node;

  //! keep track of the leaves we already regrafted; those are considered "fix" to avoid re-regrafting one
  vector<bool> fixed_stids; //TODO: use boost::dynamic_bitset for this to save space?

  list<TreeRelation> relations;

public:

  //! setup traversal numbers, subtree sizes, and StIds
  template<const bool change_StIds = false>
  MyTree* setup_node_infos()
  {
    MyTree::setup_node_infos<change_StIds>();
    if(change_StIds) fixed_stids.resize(StId_to_node.size(), false);
    return this;
  }

  //! creates a relation for the given tree, after synching StIds and setting up LCA queries
  void add_relation(MyTree& tree)
  {
    tree.setup_node_infos<true>();
    tree.sync_leaf_stids(*this);
    tree.lca_preprocess();
    relations.emplace_back(*this, tree);
  }

  //! get the tree relations
  const list<TreeRelation>& get_relations() const
  {
    return relations;
  }

  //! regraft a leaf above a node and update all relations
  void regraft_leaf_above(MyNode* const x, MyNode* const v)
  {
    const MyNode* const old_sibling = get_sibling(x);
    MyTree::regraft_leaf_above(x, v);
    // update clades, nodes_po, and leaves_po
    setup_node_infos();
    lca_preprocess();
    for(TreeRelation& rel: relations) rel.react_to_regraft(x, old_sibling);
  }

  void regraft_leaf_above(const StId x, const StId y){
    regraft_leaf_above(StId_to_node[x], StId_to_node[y]);
  }

};



TreeRelation::TreeRelation(const CandidateTree& _candidate, const MyTree& _tree):
  candidate(_candidate),
  tree(_tree),
  conflicts(_candidate.num_leaves()),
  MAST_table(_candidate, _tree)
{
  MAST_table.fill();
  conflicts.fill_conflicts(candidate, tree);
}

unsigned TreeRelation::get_mast_dist(list<StId>* const leaf_list) const
{
  return candidate.num_leaves() - get_mast(leaf_list);
}


