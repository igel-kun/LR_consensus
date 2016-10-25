

#pragma once

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

  unsigned score(const unsigned i, const unsigned j) const {
    return score(operator[]({i,j})); 
  }

  void update_entry(Entry& e, const unsigned u_po, const unsigned v_po, const TraceItem& trace)
  {
    unsigned existing_score = score(u_po, v_po);
    if(existing_score > score(e))
      e = construct_entry(existing_score, trace);
  }

  void leaf_entry(const MyNode* const u, const MyNode* const v, const unsigned u_po, const unsigned v_po)
  {
    // check whether q is an ancestor of the twin of p in the T2
    //cout << "["<<u_po<<", "<<v_po<<"] = ("<< is_ancestor(&v_ref,&u_ref)<<", {-1,-1})"<<endl;
    operator[]({u_po, v_po}) = construct_entry(is_ancestor(v, u), {-1, -1});
  }

  void matching_binary(const MyNode& u, const MyNode& v, Entry& target)
  {
    // there are only 2 maximal matchings to compare: left-left + right-right   and   left-right + right-left
    const unsigned u_child_pos[2] = {u.getSon(0)->getInfos().po_num, u.getSon(1)->getInfos().po_num};
    const unsigned v_child_pos[2] = {v.getSon(0)->getInfos().po_num, v.getSon(1)->getInfos().po_num};
    const unsigned left_left = score(u_child_pos[0], v_child_pos[0]) +
                               score(u_child_pos[1], v_child_pos[1]);
    const unsigned left_right = score(u_child_pos[0], v_child_pos[1]) + 
                                score(u_child_pos[1], v_child_pos[0]);
    if(left_left >= left_right)
      target = construct_entry(left_left, {v_child_pos[0], v_child_pos[1]});
    else
      target = construct_entry(left_right, {v_child_pos[1], v_child_pos[0]});
  }

  void set_entry(const MyNode* const u, const unsigned u_po, const MyNode* const v, const unsigned v_po)
  {
    //cout << "filling ["<<u_po<<", "<<v_po<<"]"<<endl;
    if(!u->isLeaf()){
      if(!v->isLeaf()){
        Entry& max_mast = operator[]({u_po, v_po});
        matching_binary(*u, *v, max_mast); // TODO: enable non-binary
        //cout << max_mast << " from matching, ";
        for(const MyNode& u_child: get_children(*u))
          update_entry(max_mast, u_child.getInfos().po_num, v_po, {u_child.getInfos().po_num, -1});
        for(const MyNode& v_child: get_children(*v)) 
          update_entry(max_mast, u_po, v_child.getInfos().po_num, {-1, v_child.getInfos().po_num});

        //cout << endl << "["<< u_po <<", "<< v_po <<"] = "<<max_mast<<endl;
      } else {
        //cout << "getting T1["<<stid(v)<<"]"<<endl;
        leaf_entry(u, T1[stid(v)], u_po, v_po);
      }
    } else {
      //cout << "arrived at leaf "<<u->getName()<<" (po-id "<<u_po<<")"<<endl;
      leaf_entry(T2[stid(u)], v, u_po, v_po);
    }
  }

public:
  const MyTree& T1;
  const MyTree& T2;

  using vector2d<Entry>::operator[];

  _DPTable(const unsigned cols, const unsigned rows, const MyTree& _T1, const MyTree& _T2):
    vector2d<Entry>(cols, rows),
    T1(_T1), T2(_T2)
  {}



  void fill(const vector<MyNode*>& nodes_T1,
            const vector<MyNode*>& nodes_T2)
  {
    //std::cout << "constructing DP table..."<<std::endl;
    for(unsigned u_po = 0; u_po < nodes_T1.size(); ++u_po){
      const MyNode* const u = nodes_T1[u_po];
      for(unsigned v_po = 0; v_po < nodes_T2.size(); ++v_po)
        set_entry(u, u_po, nodes_T2[v_po], v_po);
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
  virtual unsigned score(const TracableEntry& e) const { return e.first; }
  virtual TracableEntry construct_entry(const unsigned score, const TraceItem& trace) const { return {score, trace}; }
};


void trace_table(const vector<MyNode*>& nodes_T1,
                 const vector<MyNode*>& nodes_T2,
                 const TracableDPTable& mast_table,
                 list<StId>& result,
                 const unsigned p,
                 const unsigned q)
{
  const DPEntry& entry = mast_table[{p, q}];
  const TraceItem& trace = entry.second;
  //cout << "considering ["<<p<<", "<<q<<"]"<<endl;
  if(trace.from_leaf()){
    // if one of the vertices is a leaf and the other has the leaf in its clade, add the leaf to the result
    if(entry.first){
      const MyNode* const u = nodes_T1[p];
      result.push_back( u->isLeaf() ? stid(u) : stid(nodes_T2[q]));
      //cout << "leaf-stid: "<<result.back()<<endl;
    }
  } else if(trace.from_matching()){
    // if the best MAST for p and q comes from matching their children, recurse to each matched pair
    const MyNode* const u = nodes_T1[p];
    //unsigned size_before = result.size();
    trace_table(nodes_T1, nodes_T2, mast_table, result, u->getSon(0)->getInfos().po_num, trace.first);
    //unsigned size_middle = result.size();
    //cout << "level ["<<p<<", "<<q<<"]: (1) gained "<<size_middle-size_before<<" leaves from the entry ["<<u->getSon(0)->getInfos().po_num<<", "<<trace.first<<"]"<<endl;
    trace_table(nodes_T1, nodes_T2, mast_table, result, u->getSon(1)->getInfos().po_num, trace.second);
    //unsigned size_end = result.size();
    //cout << "level ["<<p<<", "<<q<<"]: (2) gained "<<size_end-size_middle<<" leaves from the entry ["<<u->getSon(1)->getInfos().po_num<<", "<<trace.second<<"]"<<endl;
  } else{
    // in all other cases, the TraceItem tells us where to go next
    //const unsigned size_before = result.size();
    const unsigned new_u = (trace.first == -1) ? p : trace.first;
    const unsigned new_v = (trace.second == -1) ? q : trace.second;
    trace_table(nodes_T1, nodes_T2, mast_table, result, new_u, new_v);
    //const unsigned size_end = result.size();
    //cout << "level ["<<p<<", "<<q<<"]: gained "<<size_end-size_before<<" leaves from the entry ["<<trace.first<<", "<<trace.second<<"]"<<endl;
  }
}

void trace_table(const vector<MyNode*>& nodes_T1,
                 const vector<MyNode*>& nodes_T2,
                 const TracableDPTable& mast_table,
                 list<StId>& result)
{
  return trace_table(nodes_T1, nodes_T2, mast_table, result, nodes_T1.size() - 1, nodes_T2.size() - 1);
}



