// Created by: Celine Scornavacca

#ifndef MASTRL_H_
#define MASTRL_H_

#include <cassert>
#include "agreement_kernel.h"
#include "MyTree.h"
#include "mast_SW93.h"
#include "op_factory.h"


//! compute the MAST-RL distance of the given trees
/** any two leaves with the same name should have the same StId
 * max_dist is the maximum allowed distance between the solution and any of the given trees
 * max_moves is the maximum number of leaf-prune-and-regraft steps that can be applied to candidate to find the solution
 **/
bool mastRL(CandidateTree& candidate,
               const unsigned max_dist,
               const unsigned max_moves)
{
  // step 1: get a tree Ti at maximal distance to candidate
  const TreeRelation* max_dist_rel = NULL;
  unsigned max_dist_to_candidate = 0;
  DEBUG3(cout << "Checking candidate"<<endl; candidate.pretty_print(); cout << "with "<<candidate.num_nodes()<<" nodes against "<<candidate.get_relations().size()<<" others"<<endl);

  for(const TreeRelation& rel: candidate.get_relations()){
    const unsigned dist = rel.get_mast_dist();
    if(dist > max_dist_to_candidate){
      max_dist_rel = &rel;
      max_dist_to_candidate = dist;

      // if any tree is too far from the candidate, return failure
      if(dist > max_dist + max_moves) {
        DEBUG3(cout << "distance "<<dist<<" too far, quitting"<<endl);
        return false;
      }
    }
  }
  DEBUG3(cout << "distant tree has distance "<<max_dist_to_candidate<<"; we've got "<<max_dist<<" max dist and "<<max_moves<<" moves left"<<endl);
  // if all trees are close to T1, then return a copy of the candidate as solution
  if(max_dist_to_candidate <= max_dist) return true;
  // otherwise, go on with a tree that is furthest from the candidate
  const MyTree& distant_tree = max_dist_rel->tree;

  // step 2: construct disagreement kernel between candidate and Ti
  const set<StId>* const kernel = max_dist_rel->disagreement_kernel(max_dist + max_moves);

  if(kernel){
    // step 3: enumerate all first (LPR-) steps on the way from candidate to a solution and recurse for each of them
    DEBUG3(cout << "---- level: "<<max_moves<<" ---> kernel: " << *kernel <<endl);
    for(const StId& x_id: *kernel){
      RegraftCandidates regrafts(get_regraft_ops(candidate, distant_tree, x_id, max_dist, max_moves));
      MyNode* const x = candidate[x_id];
      MyNode* const x_sibling = candidate.get_sibling(x);
      DEBUG3(cout << "=== LEVEL: "<<max_moves<<" ==== kernel leaf: "<< x_id<<" ("<<x->getName()<<") ==== "<<regrafts.size()<<" candidates ==="<<endl);
      cout << "got sibling " << regrafts.sibling<<endl;
      for(const RegraftOp& v: regrafts){
#warning TODO: cache-check here if we have already seen this candidate with at least this distance
        // do the regraft
        candidate.regraft_leaf_above(x, candidate[v]);
        if(mastRL(candidate, max_dist, max_moves - 1)) return true;
        // undo the previous regraft
        candidate.regraft_leaf_above(x, x_sibling);
      }// for each candidate produced by regrafting x
      DEBUG3(cout << "done processing candidates"<<endl);
    }// for each node x in the kernel
  }// if the kernelization didn't conclude that it's a no-instance
  return false;
}


//! without a d, find a consensus tree that minimizes d
void mastRL(CandidateTree& candidate)
{
  // the smallest possible d is half the MAST distance between the candidate and anyone else
  unsigned d = 1;
  for(const TreeRelation& rel: candidate.get_relations()) d = max(d, (rel.get_mast_dist()+1) / 2 );
  --d;
  do{
    ++d;
    DEBUG1(cout << "computing consensus for "<<d<<endl);
  } while( !mastRL(candidate, d, d));
}



#endif /*MASTRL_H_*/
