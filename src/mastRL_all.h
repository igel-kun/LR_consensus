// Created by: Celine Scornavacca

#ifndef MASTRL_H_
#define MASTRL_H_

#include <climits>
#include <cassert>
#include "agreement_kernel.h"
#include "MyTree.h"
#include "mast_SW93.h"
#include "candidate_factory.h"


typedef pair<unsigned, unsigned> MAST_score; //< the mast score is the max and, afterwards, the sum of mast distances
//! a solution keeps track of its score, a representative tree and how many equivalent ones we've seen
typedef pair<MAST_score, pair<MyTree*, unsigned> > Solution;

#define NO_SCORE MAST_score({UINT_MAX, UINT_MAX})
const Solution NO_SOLUTION = {NO_SCORE, {NULL, 0}};

MAST_score operator+(const MAST_score& sc1, const MAST_score& sc2)
{
  return make_pair(sc1.first + sc2.first, sc1.second + sc2.second);
}

void register_new_solution(const MyTree& tree, const MAST_score& score, Solution& s)
{
  //cout << "got offered a solution with score "<<score<<endl;
  if(score == s.first)
    s.second.second++;// if they have equal scores, increase the count of solutions
  else if(score < s.first){
    s.first = score;
    if(s.first != NO_SCORE) delete s.second.first;
    s.second.first = new MyTree(tree);
    s.second.second = 1;
  }
}

//! compute the MAST-RL distance of the given trees
/** any two leaves with the same name should have the same StId
 * max_dist is the maximum allowed distance between the solution and any of the given trees
 * moves_left is the maximum number of leaf-prune-and-regraft steps that can be applied to candidate to find the solution
 **/
void mastRL(const MyTree& candidate,
               const vector<MyTree*>& trees,
               Solution& best_solution,
               const unsigned past_moves,
               set<StId>& forbidden)
{
  // step 1: get a tree Ti at maximal distance to candidate
  unsigned max_dist_index = 0;
  MAST_score current_score(0 ,0);
  const MAST_score& best_score = best_solution.first;
  if(past_moves > best_score.first) return; // if we already took too long to get here, return failure
  const unsigned moves_left = best_score.first - past_moves;
  //cout << "Checking candidate"<<endl;
  //candidate.pretty_print();
  //cout << "with "<<candidate.num_nodes()<<" nodes"<<endl;

  for(unsigned i = 0; i < trees.size(); ++i){
    assert(trees[i]->node_infos_set_up());
    //cout << "mast against tree "<<i<<":"<<endl;
    //trees[i]->pretty_print();
    const unsigned dist = candidate.num_leaves() - mast(candidate, *trees[i]);
    current_score.second += dist;
    //cout << "got distance "<<dist<<" vs. "<<max_dist + moves_left<<endl;
    if(dist > current_score.first){
      max_dist_index = i;
      current_score.first = dist;

      // if this candidate cannot beat our current best, return failure
      if(current_score > best_score + MAST_score(moves_left, 0)) return;
    }
  }
  //! don't forget to take the distance to the original candidate into account
  current_score.first = max(current_score.first, past_moves);
  current_score.second += past_moves;

  //cout << "distant tree has index "<<max_dist_index<<" and distance "<<max_dist_to_candidate<<"; we've got "<<max_dist<<" max dist and "<<moves_left<<" moves left"<<endl;
  // if all trees are close to T1, then add a copy of the candidate as solution
  //NOTE: as best_score is a reference, this automatically updates the best_score!
  register_new_solution(candidate, current_score, best_solution);

  //  go on with a tree that is furthest from the candidate
  const MyTree& distant_tree = *trees[max_dist_index];

  //cout << "computing disagreement kernel between "<<endl;
  //candidate.pretty_print();
  //cout << "and (index "<<max_dist_index<<", distance "<<max_dist_to_candidate<<")"<<endl;
  //distant_tree.pretty_print();
  // step 2: construct disagreement kernel between candidate and Ti
#warning TODO: use the forbidden leaves to restrict the kernel before computing it
  set<StId>* const kernel = disagreement_kernel(candidate, distant_tree, 2 * best_score.first - past_moves);

  if(kernel){
    // erase forbidden StIds
    for(StId id: forbidden) kernel->erase(id);
    // step 3: enumerate all first (LPR-) steps on the way from candidate to a solution and recurse for each of them
    //cout << "---- level: "<<moves_left<<" ---> kernel: " << *kernel <<endl;
    for(const StId& x: *kernel){
      forbidden.insert(x);
      CandidateFactory fac(candidate, distant_tree, x, best_score.first, best_score.first - past_moves);
      //cout << "==== kernel leaf: "<< x<<" ("<<candidate[x]->getName()<<") ====== "<<fac.size()<<" candidates ====="<<endl;
      for(const MyTree& new_candidate: fac){
#warning TODO: cache-check here if we have already seen this candidate with at least this distance
        mastRL(new_candidate, trees, best_solution, past_moves + 1, forbidden);
      }// for each candidate produced by regrafting x
      //cout << "done processing candidates"<<endl;
      //NOTE: don't re-allow x here, as this may lead to exploring both (regraft x)+(regraft y) and (regraft y)+(regraft x)
    }// for each node x in the kernel
    // when backtracking, re-allow the leaves that were forbidden in this sub-search-tree
    for(const StId id: *kernel) forbidden.erase(id);
  }// if the kernelization didn't conclude that it's a no-instance
}


// why does gcc not let me set the "local variable" max_dist as default for max_moves? This declaration below can be easily automated...
Solution mastRL(const MyTree& candidate, const vector<MyTree*>& trees)
{
  Solution result = NO_SOLUTION;
  set<StId> forbidden_nodes;
  mastRL(candidate, trees, result, 0u, forbidden_nodes);
  return result;
}


#endif /*MASTRL_H_*/
