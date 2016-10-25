



#include <string>
#include <vector>
#include <iostream>
#include "utils.h"
#include "MyTree.h"
#include "mastRL.h"

using namespace std;
using namespace bpp;


void test_nearest_indices()
{
  vector<int> test_vec({9,6,0,5,7,9,1,4});
  nearest_lower<vector<int> > nearest_left_index_above(test_vec);
  nearest_lower<vector<int>, DIR_RIGHT> nearest_right_index_above(test_vec);

  for(unsigned i = 0; i < test_vec.size(); ++i) cout << test_vec[i] << " "; cout << endl;
  for(unsigned i = 0; i < test_vec.size(); ++i) {
    ssize_t left_index = nearest_left_index_above.query(i);
    if(left_index == -1) cout << "- "; else cout << test_vec[left_index] << " "; 
  }
  cout << endl;
  for(unsigned i = 0; i < test_vec.size(); ++i) {
    ssize_t right_index = nearest_right_index_above.query(i);
    if(right_index == -1) cout << "- "; else cout << test_vec[right_index] << " "; 
  }
  cout << endl;
}

void test_LCA(MyTree& t)
{
  //  std::cout << "testing LCAs..." << std::endl;
	vector<MyNode*> leaves = t.getLeaves();
  for(unsigned i = 0; i < leaves.size(); ++i)
    for(unsigned j = i; j < leaves.size(); ++j){
      cout << "LCA("<<leaves[i]->getName()<<", "<<leaves[j]->getName()<<") = "<<std::flush;
      const auto n =t.getLCA(leaves[i], leaves[j]);
      cout <<n->getId()<<" (stid "<<stid(n)<<")"<<endl;
    }

}

void test_induced_subtree(const MyTree& t)
{
  std::cout << "computing induced subtree..." << std::endl;
  MyTree* restricted = t.induced_subtree({3, 4, 9, 12});
  std::cout << TreeTemplateTools::treeToParenthesis(*restricted) << std::endl;
}

int main(int argc, char** argv){
  vector<MyTree*> trees;
	trees = readTrees(argv[1]);
	unsigned d = atoi(argv[2]);
  MyTree* t = trees.back();
  trees.pop_back();


  // get nodes of t in post order; this is important for mast
  cout << "computing PO traversal..."<< endl;
  vector<MyNode*> nodes;
  t->setup_node_infos(true, &nodes);
  //t->pretty_print();
  
  // sync stids of all leaves with that of t
  for(MyTree* T2: trees) {
    T2->setup_node_infos();
    T2->sync_leaf_stids(*t);
  }

  std::cout << "preprocessing LCAs..." << std::endl;
  t->lca_preprocess();


  cout << "nodes: ";
  for(const auto& n: nodes)
    cout << n->getId() <<": "<< (n->hasName() ? n->getName() : (string)" " ) << " ("<<stid(n)<<", "<<n->getInfos().po_num<<")  ";
  cout << endl;
  
//  test_LCA(t);
//  test_nearest_indices();
//  test_induced_subtree(t); 

  MyTree* mRL = mastRL(*t, trees, d);

  cout << "consensus: "<<endl;
  if(mRL) mRL->pretty_print(); else cout << " none"<<endl;

  
  /*
  for(MyTree* const T2: trees){
    list<StId> mast_leaves;
    std::cout << "computing 2-tree mast..." << std::endl;
    unsigned result = mast(&t, nodes, T2, &mast_leaves);
    std::cout << result << " mast leaves between "<< std::endl;
    t.pretty_print();
    std::cout << std::endl << " and " <<std::endl;
    T2->pretty_print();
    std::cout << std::endl;
    for(const auto& id: mast_leaves) std::cout << t[id]->getName() << ", ";
    cout << endl;
  }
*/

  return 0;
}
