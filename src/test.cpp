

#include "profiling.hpp"

#include <string>
#include <vector>
#include <iostream>
#include "utils.h"
#include "io.h"
#include "MyTree.h"
#include "candidate_tree.h"
#include "mastRL.h"
// From PhylLib:
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
using namespace std;
using namespace bpp;



int main(int argc, char** argv){
  CandidateTree* t = NULL;
  vector<MyTree*> trees;
	readTrees(t, trees, argv[1]);
  bool check_mode = false;
	string outputPath = "out.txt";
	if(argc>2) {
    if((string)(argv[2]) == "check")
      check_mode = true;
    else
	    outputPath= argv[2];
  }


  // NOTE: to call mastRL, we need the following preprocessing steps:
  // Step 1: setup t's node infos (depths, clades, stids, ...) & lcas
  t->setup_node_infos<true>();
  t->lca_preprocess();
  DEBUG3(t->pretty_print());
  
  // Step 2: setup node infos and sync stids of all leaves of the trees with that of t
  for(MyTree* T2: trees) t->add_relation(*T2);
 
  timer tm;
  tm.start();

  if(!check_mode){
    mastRL(*t);

    cout << "consensus: "<<endl;
    t->pretty_print();
    Newick newick;
    newick.write(*t, outputPath, true);
  } else {
    cout << "2-MASTs of "<<endl;
    t->pretty_print();
    size_t max_dist = 0;
    for(const TreeRelation& rel: t->get_relations()){
      cout << "vs"<<endl;
      rel.tree.pretty_print();

      list<StId> mast_list;
      unsigned mst = rel.get_mast(&mast_list);
      list<string> non_agreement;
      auto mast_list_iter = mast_list.begin();
      for(unsigned id = 0; id < t->num_leaves(); ++id){
        if((mast_list_iter == mast_list.end()) || (id != *mast_list_iter))
          non_agreement.push_back((*t)[id]->getName());
        else ++mast_list_iter;
      }
      assert(mst + non_agreement.size() == rel.tree.num_leaves());
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

