// Created by: Celine Scornavacca

#ifndef MASTRL_H_
#define MASTRL_H_

#include <cassert>
#include <set>
#include "agreement_kernel.h"
#include "MyTree.h"
#include "mast_SW93.h"


// we represent Edges as the StId of the lower vertex
typedef StId TreeEdge;

// forward declaration of the factory
class CandidateFactory;


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
void collect_upwards(const MyTree& T, StId x, const StId y, set<StId>& result)
{
  while(x != y){
    // bail if we've seen x before
    if(!result.insert(x).second) return;
    const MyNode* x_node = T[x];
    assert(x_node->hasFather());
    x = stid(x_node->getFather());
  }
}


//! translate a set of regraft candidates in Tm := T0 - Xm to a set of regraft candidates in T0
void translate_regraft_candidates(const vector<TreeEdge>& Tm_candidates,
                                  set<TreeEdge>& T0_candidates,
                                  const MyTree& Tm,
                                  const MyTree& T0,
                                  const list<StId>& Xm)
{
  T0_candidates.clear();
  //cout << "translating "<<Tm_candidates<<" with Xm = "<<Xm<<endl;
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
  //cout << "done translating, candidates are: "<<T0_candidates<<endl;
}



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
  //! assignment
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
  MyTree* T0;
  vector<MyNode*> nodes_T0;

  const MyTree* const Ti;
  const StId x;
  const string x_name;
  const unsigned max_dist;
  const unsigned max_moves;
  set<TreeEdge> regraft_candidates_T0;

  void init()
  {    
    list<StId> mast_list;
    MyTree* Ti_minus_x = *Ti - x;
    vector<MyNode*> nodes_Ti_minus_x;
    Ti_minus_x->setup_node_infos(false, &nodes_Ti_minus_x);
    mast(*T0, nodes_T0, *Ti_minus_x, nodes_Ti_minus_x, &mast_list);
    //NOTE: mast()s response is the leaves of a mast in post-order of T0
    delete Ti_minus_x;
    //cout << "2-mast: "<<mast_list<<" in "<<endl;
    //T0->pretty_print(cout, true);

    // here, Xm is actually Xm' in the paper, but we'll add x later to make it Xm
    list<StId> Xm;
    // to compute Xm, invert mast_list, making use of the fact that mast_list is in postorder with respect to T0
    
    auto mast_list_iter = mast_list.begin();
    for(const auto& u: T0->postorder_traversal()) {
      if(u.isLeaf()) {
        if((mast_list_iter == mast_list.end()) || (stid(&u) != *mast_list_iter)){
          //cout << "adding stid "<<stid(&u)<<" to Xm"<<endl;
          Xm.push_back(stid(&u));
        } else ++mast_list_iter;
      }
    }
    
    //cout << "getting Tm=T0-Xm"<<endl;
    // compute Tm = T0 - Xm (where T0 = _T0 - x)
    const MyTree* const T0_minus_Xm = T0->induced_subtree(mast_list);
    
    //cout << "getting Tm'=Ti-Xm"<<endl;
    // compute Tm_prime = Ti - Xm
    //NOTE: since T0 and Ti agree on Xm and Xm is in post-order with respect to T0, it is also in post-order with respect to Ti
    MyTree* Ti_minus_Xm = new MyTree(*Ti);
    for(const StId id: Xm) *Ti_minus_Xm -= id;

    //NOTE: from here on, Ti_minus_Xm - x and T0_minus_Xm are isomorph
    // hence, we can copy the StIds from T0-Xm to Ti_minus_Xm using MyTree::sync_stids
    // (this will give an arbitrary StId to x and its parent, but we don't care about that)
    const MyNode* const x_in_Ti_minus_Xm = (*Ti_minus_Xm)[x];
    
    Ti_minus_Xm->sync_stids_from(*T0_minus_Xm);
    Ti_minus_Xm->setup_node_infos(false);

    // use Ti - Xm to limit the set of edges to which we can regraft x
    //NOTE: due to the previous sync, these will be represented by StIds that are IDENTICAL to the StIds in T0 - Xm
    pair<const MyNode*, const set<StId> > restriction = location_restriction(*Ti_minus_Xm, x_in_Ti_minus_Xm, max_dist + max_moves);
    //cout << "the restriction is "<<(string)(restriction.first ? to_string(stid(restriction.first)) : "NULL")<<", "<<restriction.second<<" in "<<endl;
    //Ti_minus_Xm->pretty_print();

    //TODO: ask Mark what to do when y is above the root of T0-Xm, for now, we just set it to the root

    //cout << "now for"<<endl;
    //T0_minus_Xm->pretty_print();
    vector<TreeEdge> regraft_candidates_Tm;
    if(restriction.first){
      get_all_below_y_avoiding_Z(*T0_minus_Xm, (*T0_minus_Xm)[stid(restriction.first)], restriction.second, regraft_candidates_Tm);
    } else {
      const MyNode* const y = T0_minus_Xm->getRootNode();
      get_all_below_y_avoiding_Z<false>(*T0_minus_Xm, y, restriction.second, regraft_candidates_Tm);
    }
    //cout << "got "<<regraft_candidates_Tm.size()<<" regraft candidates: "<<regraft_candidates_Tm<<endl;
    //cout << "in"<<endl;
    //T0_minus_Xm->pretty_print();
    //cout << "translating to"<<endl;
    //T0->pretty_print();
    //cout << "("<<T0->num_stids()<<" stids)"<<endl;

    // translate this set of edges of T0 - Xm to a set of edges of T0
    translate_regraft_candidates(regraft_candidates_Tm, regraft_candidates_T0, *T0_minus_Xm, *T0, Xm);

    delete T0_minus_Xm;
    delete Ti_minus_Xm;
    //cout << "done building factory"<<endl;
  }

public:
  // NOTE: x should be the StId of the leaf x in both _T0 and _Ti
  CandidateFactory(const MyTree& _T0, MyTree& _Ti, const StId _x, const unsigned _max_distance, const unsigned _max_moves_in_T0):
    T0(_T0 - _x), Ti(new MyTree(_Ti)), x(_x), x_name(_T0[_x]->getName()), max_dist(_max_distance), max_moves(_max_moves_in_T0)
  {
    //cout << "building a candidate factory with x = "<<x<<endl;
    T0->sync_leaf_stids(*Ti);
    //cout << "trees"<<endl;
    //T0->pretty_print(cout, true);
    //cout << " ("<<T0->num_stids()<<" stids) & "<<endl;
    //_Ti.pretty_print();
    //cout << "("<<_Ti.num_stids()<<" stids)"<<endl;
    
    T0->setup_node_infos(false, &nodes_T0);
    T0->lca_preprocess();

//    for(unsigned i = 0; i < T0->num_leaves()-1; ++i) cout << "LCA("<<stid(T0->leaf_by_po_num(i))<<"["<<T0->leaf_by_po_num(i)->getId()<<"], "<<stid(T0->leaf_by_po_num(i+1))<<"["<<T0->leaf_by_po_num(i+1)->getId()<<"]) = "<<stid(T0->getLCA(T0->leaf_by_po_num(i), T0->leaf_by_po_num(i+1)))<<"["<<T0->getLCA(T0->leaf_by_po_num(i), T0->leaf_by_po_num(i+1))->getId()<<"]"<<endl;
    init();
    // cout << "produced candidates: "<<regraft_candidates_T0<<endl;
  }

  ~CandidateFactory()
  {
    delete T0;
    delete Ti;
  }

  CandidateIterator begin() const
  {
    CandidateIterator i(*this, regraft_candidates_T0.begin());
    return i;
  }

  CandidateIterator end() const
  {
    CandidateIterator i(*this, regraft_candidates_T0.end());
    return i;
  }

  // create a new candidate tree from T0 by grafting x on the uv
  MyTree* create_candidate_tree(const TreeEdge& uv) const 
  {
#warning TODO: can we skip this tree-copy and return a modified version of T0 instead?
    //cout << "creating candidate tree for TreeEdge above "<<uv<<" in "<<endl;
    //T0->pretty_print(cout, true);

    MyTree* result = new MyTree(*T0);
    const StId v_id = uv;
    MyNode* const new_x = result->graft_leaf_above(v_id);
    new_x->setName(x_name);
    result->assign_StId(new_x, x);
    result->consolidate_StIds();
    //cout << "new tree:"<<endl;
    //result->pretty_print(cout, true);

    return result;
  }

  size_t size() const
  {
    return regraft_candidates_T0.size();
  }
};


MyTree* CandidateIterator::operator->(){ 
  if(tree == NULL) tree = factory.create_candidate_tree(ParentClass::operator*());
  return tree;
}


//! compute the MAST-RL distance of the given trees
/** any two leaves with the same name should have the same StId
 * max_dist is the maximum allowed distance between the solution and any of the given trees
 * max_moves_in_T0 is the maximum number of leaf-prune-and-regraft steps that can be applied to T0 to find the solution
 **/
MyTree* mastRL(MyTree& T0,
               const vector<MyTree*>& trees,
               const vector<vector<MyNode*> >& nodes_PO,
               const unsigned max_dist,
               const unsigned max_moves_in_T0)
{
  // step 1: get a tree Ti at maximal distance to T0
  unsigned max_dist_index = 0;
  unsigned max_dist_to_T0 = 0;
  vector<MyNode*> nodes_T0;
  T0.setup_node_infos(false, &nodes_T0);
  const unsigned T0_leaves = T0.num_leaves();

  //cout << "Checking candidate"<<endl;
  //T0.pretty_print();
  //cout << "with "<<nodes_T0.size()<<" nodes"<<endl;

  for(unsigned i = 0; i < trees.size(); ++i){
    //cout << "mast against"<<endl;
    //trees[i]->pretty_print();
    const unsigned dist = T0_leaves - mast(T0, nodes_T0, *trees[i], nodes_PO[i]);
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
      MyTree* solution = mastRL(candidate, trees, nodes_PO, max_dist, max_moves_in_T0 - 1);
      if(solution != NULL) {
        //cout << "found a solution!"<<endl;
        return solution;
      } else;// cout << "no solution"<<endl;
    }
  }

  return NULL;
}

MyTree* mastRL(MyTree& T0,
               const vector<MyTree*>& trees,
               const unsigned max_dist,
               const unsigned max_moves_in_T0)
{
  vector<vector<MyNode*> > nodes_PO;
  for(const auto& t: trees){
    nodes_PO.push_back(vector<MyNode*>());
    vector<MyNode*>& current = nodes_PO.back();
    for(auto& u: t->postorder_traversal()) current.push_back(&u);
  }
  return mastRL(T0, trees, nodes_PO, max_dist, max_moves_in_T0);
}


// why does gcc not let me set the "local variable" max_dist as default for max_moves? This declaration below can be easily automated...
MyTree* mastRL(MyTree& T0,
               const vector<MyTree*>& trees,
               const unsigned max_dist)
{
  return mastRL(T0, trees, max_dist, max_dist);
}

#endif /*MASTRL_H_*/
