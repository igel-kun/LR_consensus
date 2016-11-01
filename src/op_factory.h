

#pragma once

#include <cassert>
#include <set>

#include "utils.h"
#include "mast_SW93.h"
#include "MyTree.h"
#include "candidate_tree.h"

// we represent Edges as the StId of the lower vertex
typedef StId TreeEdge;









//! find all highest descendants of x such that   max_size >= |clade(x)| >= min_size 
void large_clade_nodes_below(const MyTree& T, const MyNode* x, const unsigned min_size, const unsigned max_size, set<StId>& result)
{
  const unsigned x_clade_size = clade(x).size();
  if(x_clade_size >= min_size){
    if(x_clade_size <= max_size)
      result.insert(stid(x));
    else
      for(const auto& child: get_children(*x))
        large_clade_nodes_below(T, &child, min_size, max_size, result);
  }
}

//! find all highest descendants of x such that   max_size >= |clade(x)| >= min_size 
const MyNode* large_clade_node_above(const MyNode* x, const unsigned min_size)
{
  while((clade(x).size() < min_size) && (x->hasFather())) x = x->getFather();
  return x;
}


//! compute information that limits where x can be regrafted onto T0 - Xm using information of Ti - Xm (which contains x)
//NOTE: this assumes that clades have been set up!
pair<const MyNode*, set<StId> > location_restriction(const MyTree& Ti_minus_Xm,
                                            const MyNode* const x,
                                            const unsigned max_dist_and_moves)
{
  set<StId> Z;
  // step 1: walk upwards from x until we have a big enough clade
  const MyNode* const z = large_clade_node_above(x, max_dist_and_moves + 1);
  const unsigned z_clade_size = clade(z).size();
  const MyNode* const y = z->hasFather() ? large_clade_node_above(z->getFather(), max_dist_and_moves + z_clade_size) : NULL;

  // get the child of z that is not ancestor of x (whose clade does not contain x's clade)
  const MyNode* z2 = z->getSon(0);
  if(clade(z2).contains(clade(x))) z2 = z->getSon(1);
  const Clade z2_clade = clade(z2);

  // get the vertices z' below z2 such that clade(z') is large and clade(z2)\clade(z') is large
  large_clade_nodes_below(Ti_minus_Xm, z2, max_dist_and_moves, max_dist_and_moves + z2_clade.size(), Z);

  // get the child of y that is not ancestor of x (whose clade does not contain x's clade)
  if(y){
    const MyNode* y2 = y->getSon(0);
    if(clade(y2).contains(clade(x))) y2 = y->getSon(1);
    const Clade y2_clade = clade(y2);
    // get the vertices y' below y2 such that clade(y') is large and clade(y2)\clade(y') is large
    large_clade_nodes_below(Ti_minus_Xm, y2, max_dist_and_moves, max_dist_and_moves + y2_clade.size(), Z);
  }

  return make_pair(y, Z);
}


// get all edges below y but avoid diving into subtrees rooted at Z
//NOTE: this will also add the edge above y at the very end of the vector!
//NOTE: we use dont_add to avoid adding the topmost vertex... there might be a more elegant way to do that tho
template<bool dont_add = true>
void get_all_below_y_avoiding_Z(const MyTree& T, const MyNode* const y, const set<StId>& Z, vector<TreeEdge>& result)
{
  const StId y_id = stid(y);
  if(Z.find(y_id) == Z.end()){
    // step 1: put all edges under y
    for(const auto& child: get_children(*y))
      get_all_below_y_avoiding_Z<false>(T, &child, Z, result);
  }
  // step 2: put the edge above y
  if(!dont_add) result.push_back(y_id);
}


// collect all StIds in T going upwards from x, until y is reached or we reach a vertex that is already in the set, do not add y
//NOTE: this assumes that y lies above x in T
void collect_upwards(const MyTree& T, StId x, const StId y, RegraftCandidates& result)
{
  const MyNode* x_node = T[x];
  while(x != y){
    // bail if we've seen x before
    if(!result.insert(x).second) return;
    assert(x_node->hasFather());
    x_node = x_node->getFather();
    x = stid(x_node);
  }
}


//! translate a set of regraft candidates in Tm := T0 - Xm to a set of regraft candidates in T0
void translate_regraft_candidates(const vector<TreeEdge>& Tm_candidates,
                                  RegraftCandidates& T0_candidates,
                                  const MyTree& Tm,
                                  const MyTree& T0,
                                  const list<StId>& Xm)
{
  T0_candidates.clear();
  DEBUG3(cout << "translating "<<Tm_candidates<<" with Xm = "<<Xm<<endl);
  for(const TreeEdge& uv: Tm_candidates){
    const StId v_id = uv; // remember that uv is represented by v's StId
    //NOTE: v_id could be the root of Tm, in which case there is no upper bound
    StId u_id;
    if(v_id == stid(Tm.getRootNode())){
      u_id = stid(T0.getRootNode());
      T0_candidates.insert(u_id);
    } else {
      u_id = stid(Tm[v_id]->getFather());
    }
    // get the edges between u & v in T0
    collect_upwards(T0, v_id, u_id, T0_candidates);
    // get the edges between u & certain leaves l in Xm
    const MyNode* const u0 = T0[u_id];
    for(const StId& l: Xm){
      // check whether l has an ancestor v_prime STRICTLY between v0 and u0 (TODO: check with Mark if it is indeed strict)
      const MyNode* const v_prime = T0.getLCA(T0[v_id], T0[l]);
      // if so, add everything between l & u0 (that has not been added before)
      if((v_prime != u0) && (T0.getLCA(v_prime, u0) == u0))
        collect_upwards(T0, l, u_id, T0_candidates);
    }// for all leaf IDs in Xm
  }
  DEBUG3(cout << "done translating, candidates are: "<<T0_candidates<<endl);
}




//! a factory producing regraft operations
/** given a tree T0, a tree Ti, and a leaf x contained in both, try multiple regrafts of x in T0
 * the choice of regrafts is limited by the structure of Ti
 **/


// NOTE: x should be the StId of the leaf x in both candidate and _Ti
RegraftCandidates get_regraft_ops(const CandidateTree& candidate,
                                  const MyTree& Ti,
                                  const StId x,
                                  const unsigned max_dist,
                                  const unsigned max_moves)
{
  MyTree* T0(candidate - x);
  const MyNode* const x_node = candidate[x];
  const string x_name(x_node->getName());
  RegraftCandidates result(stid(get_sibling(x_node)));
    
  T0->setup_node_infos();
  T0->lca_preprocess();


  list<StId> mast_list;
  MyTree* Ti_minus_x = Ti - x;
  Ti_minus_x->setup_node_infos();

  cout << "masting"<<endl; T0->pretty_print(); cout << "and"<<endl; Ti_minus_x->pretty_print();
  mast(*T0, *Ti_minus_x, &mast_list);
  //NOTE: mast()s response uis the leaves of a mast in post-order of T0
  delete Ti_minus_x;
  DEBUG3(cout << "2-mast: "<<mast_list<<" in "<<endl;T0->pretty_print();)

  // here, Xm is actually Xm' in the paper, but we'll add x later to make it Xm
  list<StId> Xm;
  DEBUG5(cout << "checking out "<<T0->num_leaves()<<" leaves"<<endl;)
  // to compute Xm, invert mast_list, making use of the fact that mast_list is in postorder with respect to T0
  auto mast_list_iter = mast_list.begin();
  for(unsigned i = 0; i < T0->num_leaves(); ++i){
    const MyNode* const u = T0->leaf_by_po_num(i);
    if(u){
      const StId id = stid(u);
      if((mast_list_iter == mast_list.end()) || (id != *mast_list_iter)) Xm.push_back(id); else ++mast_list_iter;
    }
  }
  
  DEBUG5(cout << "getting Tm=T0-Xm"<<endl;)
  // compute Tm = T0 - Xm (where T0 = candidate - x)
  const MyTree* const T0_minus_Xm = T0->induced_subtree(mast_list);
  
  DEBUG5(cout << "getting Tm'=Ti-Xm"<<endl;)
  // compute Tm_prime = Ti - Xm
  //NOTE: since T0 and Ti agree on Xm and Xm is in post-order with respect to T0, it is also in post-order with respect to Ti
  MyTree* Ti_minus_Xm = new MyTree(Ti);
  for(const StId id: Xm) *Ti_minus_Xm -= id;

  //NOTE: from here on, Ti_minus_Xm - x and T0_minus_Xm are isomorph
  // hence, we can copy the StIds from T0-Xm to Ti_minus_Xm using MyTree::sync_stids
  // (this will give an arbitrary StId to x and its parent, but we don't care about that)
  const MyNode* const x_in_Ti_minus_Xm = (*Ti_minus_Xm)[x];
  
  Ti_minus_Xm->sync_stids_from(*T0_minus_Xm);
  Ti_minus_Xm->setup_node_infos();

  // use Ti - Xm to limit the set of edges to which we can regraft x
  //NOTE: due to the previous sync, these will be represented by StIds that are IDENTICAL to the StIds in T0 - Xm
  pair<const MyNode*, const set<StId> > restriction = location_restriction(*Ti_minus_Xm, x_in_Ti_minus_Xm, max_dist + max_moves);
  DEBUG3(cout << "the restriction is "<<(string)(restriction.first ? to_string(stid(restriction.first)) : "NULL")<<", "<<restriction.second<<" in "<<endl;Ti_minus_Xm->pretty_print();)

  DEBUG5(cout << "now for"<<endl;T0_minus_Xm->pretty_print();)
  vector<StId> regraft_candidates_Tm;
  if(restriction.first){
    get_all_below_y_avoiding_Z(*T0_minus_Xm, (*T0_minus_Xm)[stid(restriction.first)], restriction.second, regraft_candidates_Tm);
  } else {
    const MyNode* const y = T0_minus_Xm->getRootNode();
    get_all_below_y_avoiding_Z<false>(*T0_minus_Xm, y, restriction.second, regraft_candidates_Tm);
  }
  DEBUG3(cout << "got "<<regraft_candidates_Tm.size()<<" regraft candidates: "<<regraft_candidates_Tm<<endl;)
  DEBUG5(cout << "in"<<endl;T0_minus_Xm->pretty_print(); cout << "translating to"<<endl; T0->pretty_print();)

  // translate this set of edges of T0 - Xm to a set of edges of T0
  translate_regraft_candidates(regraft_candidates_Tm, result, *T0_minus_Xm, *T0, Xm);

  delete T0_minus_Xm;
  delete Ti_minus_Xm;
  delete T0;
  cout << "returning with sibling stid "<<result.sibling<<endl;
  return result;
}




