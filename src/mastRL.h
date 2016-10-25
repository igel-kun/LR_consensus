// Created by: Celine Scornavacca

#ifndef MASTRL_H_
#define MASTRL_H_

#include <cassert>
#include "agreement_kernel.h"
#include "MyTree.h"
#include "mast_SW93.h"
#include "candidate_factory.h"


//! compute the MAST-RL distance of the given trees
/** any two leaves with the same name should have the same StId
 * max_dist is the maximum allowed distance between the solution and any of the given trees
 * max_moves_in_T0 is the maximum number of leaf-prune-and-regraft steps that can be applied to T0 to find the solution
 **/
MyTree* mastRL(MyTree& T0,
               const vector<MyTree*>& trees,
               const unsigned max_dist,
               const unsigned max_moves_in_T0)
{
  // step 1: get a tree Ti at maximal distance to T0
  unsigned max_dist_index = 0;
  unsigned max_dist_to_T0 = 0;
  T0.setup_node_infos(false);
  const unsigned T0_leaves = T0.num_leaves();

  //cout << "Checking candidate"<<endl;
  //T0.pretty_print();
  //cout << "with "<<T0.num_nodes()<<" nodes"<<endl;

  for(unsigned i = 0; i < trees.size(); ++i){
    assert(trees[i]->node_infos_set_up());
    //cout << "mast against"<<endl;
    //trees[i]->pretty_print();
    const unsigned dist = T0_leaves - mast(T0, *trees[i]);
    if(dist > max_dist_to_T0){
      max_dist_index = i;
      max_dist_to_T0 = dist;

      // if any tree is too far from T1, return failure
      if(max_dist_to_T0 > max_dist + max_moves_in_T0) return NULL;
    }
  }
  //cout << "distant tree has index "<<max_dist_index<<" and distance "<<max_dist_to_T0<<"; we've got "<<max_dist<<" max dist and "<<max_moves_in_T0<<" moves left"<<endl;
  // if all trees are close to T1, then return a copy of T1 as solution
  if(max_dist_to_T0 <= max_dist) return new MyTree(T0);
  // otherwise, go on with a tree that is furthest from T1
  MyTree& distant_tree = *trees[max_dist_index];

  //cout << "computing disagreement kernel between "<<endl;
  //T0.pretty_print();
  //cout << "and (index "<<max_dist_index<<", distance "<<max_dist_to_T0<<")"<<endl;
  //distant_tree.pretty_print();
  // step 2: construct disagreement kernel between T0 and Ti
  const set<StId> kernel = disagreement_kernel(T0, distant_tree, max_dist);

  // step 3: enumerate all first (LPR-) steps on the way from T0 to a solution and recurse for each of them
  //cout << "---- level: "<<max_moves_in_T0<<" ---> kernel: " << kernel <<endl;
  for(const StId& x: kernel){
    const CandidateFactory fac(T0, distant_tree, x, max_dist, max_moves_in_T0);
    //cout << "===== LEVEL: "<<max_moves_in_T0<<" ==== kernel leaf: "<< x<<" ====== "<<fac.size()<<" candidates ====="<<endl;
    for(MyTree& candidate: fac){
      //cout << "candidate: "<<endl;
      //candidate.pretty_print();
#warning TODO: cache-check here if we have already seen this candidate with at least this distance
      MyTree* solution = mastRL(candidate, trees, max_dist, max_moves_in_T0 - 1);
      if(solution != NULL) {
        //cout << "found a solution!"<<endl;
        return solution;
      } else;// cout << "no solution"<<endl;
    }
  }

  return NULL;
}


// why does gcc not let me set the "local variable" max_dist as default for max_moves? This declaration below can be easily automated...
MyTree* mastRL(MyTree& T0,
               const vector<MyTree*>& trees,
               const unsigned max_dist)
{
  return mastRL(T0, trees, max_dist, max_dist);
}

#endif /*MASTRL_H_*/
