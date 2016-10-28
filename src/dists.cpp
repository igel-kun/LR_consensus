

#include "profiling.hpp"

#include <string>
#include <vector>
#include <iostream>
#include "utils.h"
#include "MyTree.h"
//#include "mastRL.h"
#include "mastRL_all.h"
// From PhylLib:
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
using namespace std;
using namespace bpp;

int main(int argc, char** argv){
  vector<MyTree*> trees;
	trees = readTrees(argv[1]);

  // NOTE: to call mastRL, we need the following preprocessing steps:
  // Step 1: seperate a candidate t from the trees
  MyTree* t = trees.back();
  trees.pop_back();
  // Step 2: setup t's node infos (depths, clades, stids, ...) & lcas
  t->setup_node_infos(true);
  t->setup_triplets();
  t->lca_preprocess();
  //t->pretty_print();
  return 0;
}


