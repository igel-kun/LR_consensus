

#pragma once

#include <set>

#include "utils.h"
#include "vector2d.h"

#include "MyTree.h"

using namespace bpp;

//! in order to reconstruct the resulting MAST, we need to be able to traceback on the DP table
/** each trace item with index [p,q] is a pair (x,y) of integers such that,
 * if x == p, then we've got the maximum from [p,y]
 * if y == q, then we've got the maximum from [x,q]
 * otherwise, we've got the maximum from matching the first child of p to x and the second child to y
 * NOTE: this assumes that the order of children does not change between creating the trace and reading it.
 * */
struct TraceItem: public pair<ssize_t, ssize_t>
{
  using pair<ssize_t, ssize_t>::pair;
  // we can deduce wether this item corresponds to a matching item or not
  bool from_matching() const
  {
    return (first != -1) && (second != -1);
  }
  // we can deduce that the item comes from matching two leaves together if the trace points to itself
  bool from_leaf() const
  {
    return (first == -1) && (second == -1);
  }
};

typedef pair<unsigned, TraceItem> DPEntry;

template<typename Entry>
class _DPTable: public vector2d<Entry>
{
  using vector2d<Entry>::vector2d;

  virtual unsigned score(const Entry& e) const = 0;
  virtual Entry construct_entry(const unsigned score, const TraceItem& trace) const = 0;

  unsigned score(const StId i, const StId j) const
  {
    return score(operator[]({i,j})); 
  }

  void update_entry(Entry& e, const StId u_id, const StId v_id, const TraceItem& trace)
  {
    unsigned existing_score = score(u_id, v_id);
    if(existing_score > score(e))
      e = construct_entry(existing_score, trace);
  }

  unsigned leaf_entry(const MyNode& u, const MyNode& v)
  {
    // check whether q is an ancestor of the twin of p in the T2
    const unsigned _score = (is_ancestor(v, u) || is_ancestor(u,v));
    const StId u_id = u.getInfos().stId;
    const StId v_id = v.getInfos().stId;
    DEBUG5(cout << "["<<u_id<<", "<<v_id<<"] = "<< _score <<" (leaf)"<<endl);
    operator[]({u_id, v_id}) = construct_entry(_score, {-1, -1});
    return _score;
  }

  void matching_binary(const MyNode& u, const MyNode& v, Entry& target)
  {
    // there are only 2 maximal matchings to compare: left-left + right-right   and   left-right + right-left
    const unsigned u_child_ids[2] = {stid(u.getSon(0)), stid(u.getSon(1))};
    const unsigned v_child_ids[2] = {stid(v.getSon(0)), stid(v.getSon(1))};
    const unsigned left_left = score(u_child_ids[0], v_child_ids[0]) +
                               score(u_child_ids[1], v_child_ids[1]);
    const unsigned left_right = score(u_child_ids[0], v_child_ids[1]) + 
                                score(u_child_ids[1], v_child_ids[0]);
    if(left_left >= left_right)
      target = construct_entry(left_left, {v_child_ids[0], v_child_ids[1]});
    else
      target = construct_entry(left_right, {v_child_ids[1], v_child_ids[0]});
  }

  unsigned set_entry(const MyNode& u, const MyNode& v)
  {
    const StId u_id = stid(&u);
    const StId v_id = stid(&v);
    if(!u.isLeaf()){
      if(!v.isLeaf()){
        Entry& max_mast = operator[]({u_id, v_id});
        matching_binary(u, v, max_mast); // TODO: enable non-binary
        DEBUG5(cout << "\t" << max_mast << " from matching the sons "<<stid(u.getSon(0))<<" & "<<stid(u.getSon(1)));
        for(const MyNode& u_child: get_children(u))
          update_entry(max_mast, stid(&u_child), v_id, {stid(&u_child), -1});
        for(const MyNode& v_child: get_children(v)) 
          update_entry(max_mast, u_id, stid(&v_child), {-1, stid(&v_child)});

        DEBUG5(cout << endl << "["<< u_id <<", "<< v_id <<"] = "<<max_mast<<endl);
        return score(max_mast);
      } else return leaf_entry(u, *T1[v_id]); 
    } else 
      return leaf_entry(*T2[u_id], v);
  }

public:
  const MyTree& T1;
  const MyTree& T2;

  using vector2d<Entry>::operator[];

  _DPTable(const MyTree& _T1, const MyTree& _T2):
    vector2d<Entry>(_T1.max_StId() + 1, _T2.max_StId() + 1),
    T1(_T1), T2(_T2)
  {
    assert(T1.node_infos_set_up()); // check that the postorder numbers of T1 have been set up
    assert(T2.node_infos_set_up());
  }

  void fill()
  {
    DEBUG5(std::cout << "constructing DP table..."<<std::endl);
    for(PostOrderConstIterator u_iter(T1.getRootNode()); u_iter.is_valid(); ++u_iter)
      for(PostOrderConstIterator v_iter(T2.getRootNode()); v_iter.is_valid(); ++v_iter)
        set_entry(*u_iter, *v_iter);
    DEBUG5(std::cout << "done constructing DP table."<<std::endl);
  }

  const Entry& get_root_entry() const
  {
    const StId T1_root_id = stid(T1.getRootNode());
    const StId T2_root_id = stid(T2.getRootNode());
    DEBUG5(cout << "mast table root entry ["<<T1_root_id<<", "<<T2_root_id<<"] = " << score(T1_root_id, T2_root_id)<< endl);
    return operator[]({T1_root_id, T2_root_id});
  }
  const unsigned get_root_score() const
  {
    return score(get_root_entry());
  }

protected:
  //! recompute an entry and return whether it changed
  bool update_entry(const MyNode* const x, const MyNode* const y)
  {
    const unsigned old_score = score(operator[]({stid(x), stid(y)}));
    const unsigned new_score = set_entry(*x, *y);
    cout << "updated "<<old_score<<" to "<<new_score<<endl;
    return (old_score != new_score);
  }
  //! upadte the table entries between x in T1 and upwards from from_T2 in T2
  void update_upwards2(const MyNode* const x, const MyNode* from_T2)
  {
    while(from_T2->hasFather()) {
      if(!update_entry(x, from_T2)) return;
      from_T2 = from_T2->getFather();
    }
    set_entry(*x, *from_T2);
  }
  //! update the table entries upwards in T1 from from_T1 to (including) to_T1, and in T2 from from_T2 to the root
  void update_upwards(const MyNode* from_T1, const MyNode* const to_T1, const MyNode* from_T2)
  {
    if(!from_T1) return;
    while(from_T1->hasFather() && (from_T1 != to_T1)){
      cout << "updating upwards from ("<<stid(from_T1)<<" ,"<<stid(from_T2)<<")"<<endl;
      update_upwards2(from_T1, from_T2);
      from_T1 = from_T1->getFather();
    }
    update_upwards2(from_T1, from_T2);
  }

  //! update the table for a node x of T1 up from the leaves of T2
  void update_from_leaves(const MyNode* const x)
  {
    DepthCompare<true> cmp(T2);
    std::set<StId, DepthCompare<true> > todo_T2(cmp);
    for(StId i = 0; i < T2.num_leaves(); ++i) todo_T2.insert(i);
    while(!todo_T2.empty()){
      const StId y_id = *(todo_T2.begin());
      const MyNode* const y = T2[y_id];
      const MyNode* const parent = (y->hasFather() ? y->getFather() : NULL);
      todo_T2.erase(todo_T2.begin());
      if(update_entry(x, y))
        if(parent) todo_T2.insert(stid(parent));
    }
  }

public:
  //! update the DP table AFTER a regraft of leaf, provided the StIds didn't change!
  //NOTE: the nodes should be in T1
  void react_to_regraft(const MyNode* const leaf, const MyNode* const old_sibling)
  {
    assert(leaf->hasFather());
    const MyNode* const old_grand_parent = (old_sibling->hasFather() ? old_sibling->getFather() : NULL);
    const MyNode* const parent = leaf->getFather();
    const MyNode* const new_grand_parent = (parent->hasFather() ? parent->getFather() : NULL);
    const MyNode* const lca = T1.getLCA(parent, old_sibling);
    const MyNode* const leaf_in_T2 = T2[stid(leaf)];
    cout << "MAST with "<<endl;
    T2.pretty_print();
    cout << "updating DP table after regraft of "<<leaf->getName()<<" with parent id "<<stid(parent)<<" (lca with "<<stid(old_sibling)<<" is "<<stid(lca)<<")"<<endl;
    // what changes with the regraft of x?
    //    1. [parent,v] has to be updated for many a leaf v (and above)
    //    2. x gets removed from all pairings [u,v] where u is strictly below the lca and v is above x
    //    3. x gets added to all pairings [u,v] where u is strictly above parent and v is above x
    // Step 1: update [parent, v] for leaves and above in T2
    update_from_leaves(parent);
    // Step 1: update all pairings [u,v] where u is strictly above parent and v is above x
    if(lca == old_sibling){
      // leaf got regrafted below old_sibling
      update_upwards(new_grand_parent, NULL, leaf_in_T2);
    } else if(lca == parent) {
      // leaf got regrafted above old_sibling
      update_upwards(old_grand_parent, NULL, leaf_in_T2);
    } else {
      // leaf got regrafted into an unrelated subtree
      update_upwards(old_grand_parent, lca, leaf_in_T2);
      update_upwards(new_grand_parent, lca, leaf_in_T2);
      update_upwards(lca, NULL, leaf_in_T2);
    }
  }

};

typedef pair<unsigned, TraceItem> TracableEntry;

class DPTable: public _DPTable<unsigned>{
  using _DPTable<unsigned>::_DPTable;
  virtual unsigned score(const unsigned& e) const { return e; }
  virtual unsigned construct_entry(const unsigned score, const TraceItem& trace) const { return score; }
};

class TracableDPTable: public _DPTable<TracableEntry>{
  using _DPTable<TracableEntry>::_DPTable;
  using _DPTable<TracableEntry>::T1;
  using _DPTable<TracableEntry>::T2;
  virtual unsigned score(const TracableEntry& e) const { return e.first; }
  virtual TracableEntry construct_entry(const unsigned score, const TraceItem& trace) const { return {score, trace}; }

public:
  using _DPTable<TracableEntry>::operator[];

  template<typename Container = list<StId> >
  void trace(Container& result) const
  {
    return trace(result, stid(T1.getRootNode()), stid(T2.getRootNode()));
  }

  template<typename Container = list<StId> >
  void trace(Container& result, const StId u_id, const StId v_id) const
  {
    const DPEntry& entry = operator[]({u_id, v_id});
    const TraceItem& ti = entry.second;
    DEBUG5(cout << "tracing ["<<u_id<<", "<<v_id<<"]"<<endl);
    if(ti.from_leaf()){
      // if one of the vertices is a leaf and the other has the leaf in its clade, add the leaf to the result
      if(entry.first){
        const MyNode* const u = T1[u_id];
        result.push_back( u->isLeaf() ? u_id : v_id);
        DEBUG5(cout << "leaf-stid: "<<result.back()<<endl);
      }
    } else if(ti.from_matching()){
      // if the best MAST for p and q comes from matching their children, recurse to each matched pair
      const MyNode* const u = T1[u_id];
      DEBUG5(unsigned size_before = result.size());
      trace(result, stid(u->getSon(0)), ti.first);
      DEBUG5(unsigned size_middle = result.size());
      DEBUG5(cout << "level ["<<u_id<<", "<<v_id<<"]: (1) gained "<<size_middle-size_before<<" leaves from the entry ["<<stid(u->getSon(0))<<", "<<ti.first<<"]"<<endl);
      trace(result, stid(u->getSon(1)), ti.second);
      DEBUG5(unsigned size_end = result.size());
      DEBUG5(cout << "level ["<<u_id<<", "<<v_id<<"]: (2) gained "<<size_end-size_middle<<" leaves from the entry ["<<stid(u->getSon(1))<<", "<<ti.second<<"]"<<endl);
    } else{
      // in all other cases, the TraceItem tells us where to go next
      DEBUG5(const unsigned size_before = result.size());
      const StId new_u = (ti.first == -1) ? u_id : ti.first;
      const StId new_v = (ti.second == -1) ? v_id : ti.second;
      trace(result, new_u, new_v);
      DEBUG5(const unsigned size_end = result.size());
      DEBUG5(cout << "level ["<<u_id<<", "<<v_id<<"]: gained "<<size_end-size_before<<" leaves from the entry ["<<ti.first<<", "<<ti.second<<"]"<<endl);
    }
  }

};






