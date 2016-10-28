

#include "profiling.hpp"

#include <string>
#include <vector>
#include <iostream>
#include "utils.h"
#include "MyTree.h"
#include "mastRL.h"
// From PhylLib:
#include <Bpp/Phyl/Io/IoTree.h>
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
	string outputPath = "out.txt";
  bool check_mode = false;
	if(argc>2) {
    if((string)(argv[2]) == "check")
      check_mode = true;
    else
	    outputPath= argv[2];
  }


  // NOTE: to call mastRL, we need the following preprocessing steps:
  // Step 1: seperate a candidate t from the trees
  MyTree* t = trees.back();
  trees.pop_back();
  // Step 2: setup t's node infos (depths, clades, stids, ...) & lcas
  t->setup_node_infos(true);
  t->setup_triplets();
  t->lca_preprocess();
  //t->pretty_print();
  
  // Step 3: setup node infos and sync stids of all leaves of the trees with that of t
  for(MyTree* T2: trees) {
    T2->setup_node_infos();
    T2->sync_leaf_stids(*t);
    T2->setup_triplets();
    T2->lca_preprocess();
  }
 
//  test_LCA(t);
//  test_nearest_indices();
//  test_induced_subtree(t); 
//  test_2mast(t, trees);

  timer tm;
  tm.start();

  if(!check_mode){
    MyTree* mRL = mastRL(*t, trees);

    cout << "consensus: "<<endl;
    if(mRL) {
      mRL->pretty_print();
/*
      mRL->setup_node_infos(false);
      trees.push_back(t);
      for(unsigned i = 0; i < trees.size(); ++i){
        //cout << "agreement with:"<<endl;
        //trees[i]->pretty_print();
        list<StId> mast_list;
        mast(*mRL, *trees[i], &mast_list);
        list<string> non_agreement;
        auto mast_list_iter = mast_list.begin();
        for(unsigned j = 0; j < mRL->num_leaves(); ++j){
          const StId id = stid(mRL->leaf_by_po_num(j));
          if((mast_list_iter == mast_list.end()) || (id != *mast_list_iter)) non_agreement.push_back((*mRL)[id]->getName()); else ++mast_list_iter;
        }
        cout << "mast with trees["<<i<<"]: "<<mast_list<<" (disagreeing on "<<non_agreement<<")"<<endl;
      }
      */
    }
    
    if(mRL) {
      cout << " OK"<<endl;
      Newick newick;
      newick.write(* mRL, outputPath, true);
    } else cout << " none"<<endl;
  } else {
    cout << "2-MASTs of "<<endl;
    t->pretty_print();
    size_t max_dist = 0;
    for(MyTree* other: trees){
      cout << "vs"<<endl;
      other->pretty_print();

      list<StId> mast_list;
      unsigned mst = mast(*t, *other, &mast_list);
      list<string> non_agreement;
      auto mast_list_iter = mast_list.begin();
      for(unsigned j = 0; j < t->num_leaves(); ++j){
        const StId id = stid(t->leaf_by_po_num(j));
        if((mast_list_iter == mast_list.end()) || (id != *mast_list_iter))
          non_agreement.push_back((*t)[id]->getName()); else ++mast_list_iter;
      }
      assert(mst + non_agreement.size() == other->num_leaves());
      cout << mst << " mast ("<<mast_count<<"): "<<mast_list<<" (disagreeing on "<<non_agreement<<")"<<endl;
      max_dist = std::max(max_dist, non_agreement.size());
    }
    cout << "max-distance: "<<max_dist<<endl;
  }

  tm.stop();
  if(tm.seconds_passed())
    cout << "computed "<<mast_count<<" 2-masts in "<<tm.seconds_passed()<<"s = "<<((float)mast_count)/((float)tm.seconds_passed())<<" mast/s"<<endl;

  return 0;
}





/*
(((N,C),L),((A,(((B,I),D),J)),(((H,M),E),((O,P),(Q,((R,S),T))))));
*/
