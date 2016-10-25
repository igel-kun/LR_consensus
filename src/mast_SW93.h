
#ifndef MAST_H_
#define MAST_H_

#include <cassert>
#include "mastDP.h"
#include "MyTree.h"

using namespace bpp;



// TODO: preprocess the trees to find same agreeing clades & work on them independently
// TODO: only compare node x of T1 to nodes of T2 in the subtree of T2 induced by the clade of x
//! compute the size of a maximum agreement subtree, optinally returning the leaves in such an agreement forest
/** NOTE: this assumes that all trees' postorder numbers, clades, and StIds have been set up
 * and each two leaves have the same name iff they have the same StId
 * NOTE: we don't need the LCA's for the 2-mast so far
 */
unsigned mast(const MyTree& T1,
              const vector<MyNode*>& nodes_T1,
              const MyTree& T2,
              const vector<MyNode*>& nodes_T2,
              list<StId>* mast_leaves = NULL)
{
  assert(T1.node_infos_set_up()); // check that the postorder numbers of T1 have been set up
  assert(T2.node_infos_set_up());
  // T2->pretty_print();

  //cout << "T1: "; for(const auto& u: nodes_T1) cout << stid(u)<<"("<<u->getInfos().po_num<<") "; cout << endl;
  //cout << "T2: "; for(const auto& u: nodes_T2) cout << stid(u)<<"("<<u->getInfos().po_num<<") "; cout << endl;
  //for(const auto& u: nodes_T1) if(u->isLeaf()) assert(u->getName() == T2[stid(u)]->getName());

  if(mast_leaves){
    // if we want to output the remaining leaves, we need to keep a trace
    TracableDPTable mast_table(nodes_T1.size(), nodes_T2.size(), T1, T2);
    mast_table.fill(nodes_T1, nodes_T2);
    //cout << "done computing table, root entry = "<<mast_table.back().first<<", now tracing..."<<endl;
    trace_table(nodes_T1, nodes_T2, mast_table, *mast_leaves);
    //cout << "traced "<<mast_leaves->size()<<" leaves."<<endl;
    return mast_table.back().first;
  } else {
    DPTable mast_table(nodes_T1.size(), nodes_T2.size(), T1, T2);
    //cout << "building 2-mast table of size "<<mast_table.size()<<" for"<<endl;
    //T1->pretty_print(cout, true);
    //cout << " and "<<endl;
    //T2->pretty_print(cout, true);

    mast_table.fill(nodes_T1, nodes_T2);
    return mast_table.back();
  }
}


#endif /*MAST_H_*/
