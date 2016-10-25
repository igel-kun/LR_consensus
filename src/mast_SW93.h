
#ifndef MAST_H_
#define MAST_H_

#include <cassert>
#include "mastDP.h"
#include "MyTree.h"

using namespace bpp;



// TODO: preprocess the trees to find same agreeing clades & work on them independently
// TODO: only compare node x of T1 to nodes of T2 in the subtree of T2 induced by the clade of x
//! compute the size of a maximum agreement subtree, optinally returning the leaves in such an agreement forest
/** NOTE: this assumes that T1s nodes' postorder numbers have been set up!
 */
unsigned mast(const MyTree* const T1, const vector<MyNode*>& nodes_T1, MyTree* const T2, list<StId>* mast_leaves = NULL)
{
  assert(T1->node_infos_set_up()); // check that the postorder numbers of T1 have been set up
  assert(T2->node_infos_set_up());

  // do T2 preprocessing and get vertices post order
  vector<MyNode*> nodes_T2;
  T2->setup_node_infos(true, &nodes_T2);
  T2->sync_leaf_stids(*T1);
  // T2->pretty_print();

  if(mast_leaves){
    // if we want to output the remaining leaves, we need to keep a trace
    TracableDPTable mast_table(nodes_T1.size(), nodes_T2.size(), *T1, *T2);
    mast_table.fill(nodes_T1, nodes_T2);
    //cout << "done computing table, root entry = "<<mast_table.back().first<<", now tracing..."<<endl;
    trace_table(nodes_T1, nodes_T2, mast_table, *mast_leaves);
    //cout << "traced "<<mast_leaves->size()<<" leaves."<<endl;
    return mast_table.back().first;
  } else {
    DPTable mast_table(nodes_T1.size(), nodes_T2.size(), *T1, *T2);
    //cout << "building 2-mast table of size "<<mast_table.size()<<" for"<<endl;
    //T1->pretty_print(cout, true);
    //cout << " and "<<endl;
    //T2->pretty_print(cout, true);

    mast_table.fill(nodes_T1, nodes_T2);
    return mast_table.back();
  }
}


#endif /*MAST_H_*/
