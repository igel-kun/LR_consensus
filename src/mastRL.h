// Created by: Celine Scornavacca

#ifndef MASTRL_H_
#define MASTRL_H_

#include <cassert>
#include <set>
#include "agreement_kernel.h"
#include "MyTree.h"
#include "MultiGraph.h"
#include "mast.h"

unsigned mast(const MyTree& T0, const MyTree& T1)
{
}
//! get the leaves to remove from T0 to break all conflicts between T0 and T1
vector<StId> mast_solution(const MyTree& T0, const MyTree& T1)
{
}

//! find all highest descendants of x such that   max_size >= |clade(x)| >= min_size 
void large_clade_nodes_below(const MyTree& T, const MyNode* x, const unsigned min_size, const unsigned max_size, set<StId>& result)
{
  const unsigned x_clade_size = x->getInfos().getClade().size();
  if(x_clade_size >= min_size){
    if(x_clade_size <= max_size)
      result.insert(x->getInfos().getStId());
    else
      for(const auto& child: T.get_children(*x))
        large_clade_nodes_below(T, &child, min_size, max_size, result);
  }
}

//! compute information that limits where x can be regrafted onto Ti - Xm (which contains x)
//NOTE: this assumes that clades have been set up!
pair<StId, set<StId> > location_restriction(const MyTree& Ti_minus_Xm,
                                            const MyNode* x,
                                            const unsigned max_dist_and_moves)
{
  set<StId> Z;
  // step 1: walk upwards from x until we have a big enough clade
  const MyNode* z = x;
  unsigned z_clade_size;
  // sorry for assignment-operator comparison causing dense code...
  while((z_clade_size = z->getInfos().getClade().size()) < max_dist_and_moves + 1) z = z->getFather();
  // step 2: walk upwards from z until we have a big enough clade
  const MyNode* y = z->getFather();
  while(y->getInfos().getClade().size() < max_dist_and_moves + z_clade_size) y = y->getFather();

  // get the child of z that is not ancestor of x (whose clade does not contain x's clade)
  const MyNode* z2 = z->getSon(0);
  if(z2->getInfos().getClade().contains(x->getInfos().getClade())) z2 = z->getSon(1);
  const Clade z2_clade = z2->getInfos().getClade();
  // get the vertices z' below z2 such that clade(z') is large and clade(z2)\clade(z') is large
  large_clade_nodes_below(Ti_minus_Xm, z2, max_dist_and_moves, max_dist_and_moves + z2_clade.size(), Z);

  // get the child of y that is not ancestor of x (whose clade does not contain x's clade)
  const MyNode* y2 = y->getSon(0);
  if(y2->getInfos().getClade().contains(x->getInfos().getClade())) y2 = y->getSon(1);
  const Clade y2_clade = y2->getInfos().getClade();
  // get the vertices y' below y2 such that clade(y') is large and clade(y2)\clade(y') is large
  large_clade_nodes_below(Ti_minus_Xm, y2, max_dist_and_moves, max_dist_and_moves + y2_clade.size(), Z);

  return make_pair(y->getInfos().getStId(), Z);
}

// we represent Edges as the StId of the lower vertex
typedef StId TreeEdge;

// get all edges below y but avoid diving into subtrees rooted at Z
//NOTE: this will also add the edge above y at the very end of the vector!
void get_all_below_y_avoiding_Z(const MyTree& T, const StId y, const set<StId>& Z, vector<TreeEdge>& result)
{
  if(Z.find(y) == Z.end()){
    // step 1: put all edges under y
    for(const auto& child: T.get_children(y))
      get_all_below_y_avoiding_Z(T, child.getInfos().getStId(), Z, result);
  }
  // step 2: put the edge above y
  result.push_back(y);
}

// collect all StIds in T going upwards from x, until y is reached or we reach a vertex that is already in the set, do not add y
//NOTE: this assumes that y lies above x in T
void collect_upwards(const MyTree& T, StId x, const StId y, set<StId>& result)
{
  while(x != y){
    // bail if we've seen x before
    if(!result.insert(x).second) return;
    const MyNode* x_node = T.getNodeWithStId(x);
    assert(x_node->hasFather());
    x = x_node->getFather()->getInfos().getStId();
  }
}

// translate a set of regraft candidates in Tm := T0 - Xm to a set of regraft candidates in T0
void translate_regraft_candidates(const vector<TreeEdge>& Tm_candidates,
                                  set<TreeEdge>& T0_candidates,
                                  const MyTree& Tm,
                                  const MyTree& T0,
                                  const vector<StId>& Xm)
{
  T0_candidates.clear();
  for(const TreeEdge& uv: Tm_candidates){
    const StId v_id = uv; // remember that uv is represented by v's StId
    const StId u_id = Tm.getNodeWithStId(v_id)->getFather()->getInfos().getStId();
    // get u & v in T0
    const MyNode* u0 = T0.getNodeWithStId(u_id);
    // get the edges between u & v in T0
    collect_upwards(T0, v_id, u_id, T0_candidates);
    // get the edges between u & certain leaves l in Xm
    for(const StId& l: Xm){
      // check whether l has an ancestor v_prime between v0 and u0
      const MyNode* v_prime = T0.getLCA(v_id, l);
      // if so, add everything between l & u0 (that has not been added before)
      if(T0.getLCA(v_prime, u0) == u0)
        collect_upwards(T0, l, u_id, T0_candidates);
    }
  }
}

// forward declaration of the factory
class CandidateFactory;

//! an iterator produced by a CandidateFactory to allow iterating over candidates produced by the factory
class CandidateIterator: public set<TreeEdge>::iterator
{
  typedef typename set<TreeEdge>::iterator ParentClass;

  const CandidateFactory& factory;
  MyTree* tree;

public:
  CandidateIterator(const CandidateFactory& _factory):
    ParentClass(),
    factory(_factory),
    tree(NULL)
  {}
  CandidateIterator(const CandidateFactory& _factory, const typename set<TreeEdge>::iterator _it):
    ParentClass(_it),
    factory(_factory),
    tree(NULL)
  {}
  CandidateIterator(const CandidateIterator& _it):
    ParentClass(_it),
    factory(_it.factory),
    tree(_it.tree)
  {}


  // on dereference, ask the factory to create the tree corresponding to the candidate StId
  MyTree* operator->();
  MyTree& operator*()
  {
    return *(operator->());
  }
  //! increment operator
  CandidateIterator& operator++() 
  {
    delete tree; tree = NULL;
    ParentClass::operator++();
    return *this;
  }
  //! post-increment
  CandidateIterator operator++(int i) 
  {
    CandidateIterator result(*this);
    ++(*this);
    return result;
  }
  //! addition
  CandidateIterator& operator=(const CandidateIterator& it) 
  {
    tree = it.tree;
    ParentClass::operator=(it);
    return *this;
  }

};

//! a factory producing candidate trees
/** given a tree T0, a tree Ti, and a leaf x contained in both, try multiple regrafts of x in T0
 * the choice of regrafts is limited by the structure of Ti
 **/
class CandidateFactory
{
  const MyTree* T0;
  const MyTree& Ti;
  const StId x;
  const unsigned max_dist;
  const unsigned max_moves;
  set<TreeEdge> regraft_candidates_T0;

  void init()
  {
    MyTree* Ti_minus_x = Ti - x;
    // here, Xm is actually Xm' in the paper, but we'll add x later to make it Xm
    vector<StId> Xm = mast_solution(*T0, *Ti_minus_x);
    delete Ti_minus_x;

    // compute Tm by removing Xm from T0 (= _T0 - x)
    MyTree* T0_minus_Xm = new MyTree(*T0);
    for(const StId id: Xm) T0_minus_Xm -= id;

    // compute Tm' by removing Xm from Ti
    MyTree* Ti_minus_Xm = new MyTree(Ti);
    for(const StId id: Xm) Ti_minus_Xm -= id;

    //NOTE: from here on, Ti_minus_Xm - x and T0_minus_Xm are isomorph
    // hence, we can copy the StIds from T0-Xm to Ti_minus_Xm using MyTree::sync_stids
    // (this will give an arbitrary StId to x and its parent, but we don't care about that)
    const MyNode* x_in_Ti_minus_Xm = Ti_minus_Xm->getNodeWithStId(x);
    Ti_minus_Xm->sync_stids_from(*T0_minus_Xm);

    // use Ti - Xm to limit the set of edges to which we can regraft x
    //NOTE: due to the previous sync, these will be represented by StIds that are IDENTICAL to the StIds in T0 - Xm
    Ti_minus_Xm->setClades();
    pair<StId, set<StId> > restriction = location_restriction(*Ti_minus_Xm, x_in_Ti_minus_Xm, max_dist + max_moves);

    vector<TreeEdge> regraft_candidates_Tm;
    get_all_below_y_avoiding_Z(*T0_minus_Xm, restriction.first, restriction.second, regraft_candidates_Tm);
    regraft_candidates_Tm.pop_back(); // remember to remove the last element as this is y

    // translate this set of edges of T0 - Xm to a set of edges of T0
    Xm.push_back(x);
    translate_regraft_candidates(regraft_candidates_Tm, regraft_candidates_T0, *T0_minus_Xm, *T0, Xm);

    delete T0_minus_Xm;
    delete Ti_minus_Xm;
  }

public:
  // NOTE: x should be the StId of the leaf x in both _T0 and _Ti
  CandidateFactory(const MyTree& _T0, const MyTree& _Ti, const StId _x, const unsigned _max_distance, const unsigned _max_moves_in_T0):
    T0(_T0 - x), Ti(_Ti), x(_x), max_dist(_max_distance), max_moves(_max_moves_in_T0)
  {
    init();
  }

  ~CandidateFactory()
  {
    delete T0;
  }

  CandidateIterator begin() const
  {
    return CandidateIterator(*this, regraft_candidates_T0.begin());
  }

  CandidateIterator end() const
  {
    return CandidateIterator(*this, regraft_candidates_T0.end());
  }

  MyTree* create_candidate_tree(const TreeEdge& uv) const { 
  }

};

MyTree* CandidateIterator::operator->(){ 
  if(tree == NULL) tree = factory.create_candidate_tree(ParentClass::operator*());
  return tree;
}


CandidateFactory get_candidate_trees(const MyTree& T0, const MyTree& Ti, const StId x)
{
}


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
  for(unsigned i = 0; i < trees.size(); ++i){
    const unsigned dist = mast(T0, *trees[i]);
    if(dist > max_dist_to_T0){
      max_dist_index = i;
      max_dist_to_T0 = dist;
    }
  }
  // if all trees are close to T1, then return a copy of T1 as solution
  if(max_dist_to_T0 <= max_dist) return new MyTree(T0);
  // if any tree is too far from T1, return failure
  if(max_dist_to_T0 > max_dist + max_moves_in_T0) return NULL;
  // otherwise, go on with a tree that is furthest from T1
  MyTree& distant_tree = *trees[max_dist_index];

  // step 2: construct disagreement kernel between T0 and Ti
  const set<StId> kernel = disagreement_kernel(T0, distant_tree, max_dist);

  // step 3: enumerate all first (LPR-) steps on the way from T0 to a solution and recurse for each of them
  for(const StId& x: kernel)
    for(MyTree& candidate: get_candidate_trees(T0, distant_tree, x)){
#warning TODO: cache-check here if we've already seen this candidate with at least this distance
      MyTree* solution = mastRL(candidate, trees, max_dist, max_moves_in_T0 - 1);
      if(solution != NULL) return solution;
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
