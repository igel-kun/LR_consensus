



#include <string>
#include <vector>
#include <iostream>
#include "utils.h"
#include "MyTree.h"
#include "mastRL.h"
// From PhylLib:
#include <Bpp/Phyl/Io/IOTree.h>
#include <Bpp/Phyl/Io/Newick.h>
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

/*
void test_2mast(const MyTree& t, const vector<MyNode*>& nodes, const vector<MyTree*>& trees)
{
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
}
*/

int main(int argc, char** argv){
  vector<MyTree*> trees;
	trees = readTrees(argv[1]);
	unsigned d = atoi(argv[2]);
	string outputPath = "out.txt";
	if(argc>3)
	outputPath= argv[3];


  // NOTE: to call mastRL, we need the following preprocessing steps:
  // Step 1: seperate a candidate t from the trees
  MyTree* t = trees.back();
  trees.pop_back();
  // Step 2: setup t's node infos (depths, clades, stids, ...) & lcas
  t->setup_node_infos(true);
  t->lca_preprocess();
  //t->pretty_print();
  
  // Step 3: setup node infos and sync stids of all leaves of the trees with that of t
  for(MyTree* T2: trees) {
    T2->setup_node_infos();
    T2->sync_leaf_stids(*t);
  }
 
//  test_LCA(t);
//  test_nearest_indices();
//  test_induced_subtree(t); 
//  test_2mast(t, trees);

  MyTree* mRL = mastRL(*t, trees, d);

  cout << "consensus: "<<endl;
  //if(mRL) mRL->pretty_print(); else cout << " none"<<endl;
  
  if(mRL) {
	Newick newick;
  	newick.write(* mRL, outputPath, true);
}
  

  return 0;
}
