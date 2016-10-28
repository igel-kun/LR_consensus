

#include "profiling.hpp"

#include <string>
#include <vector>
#include <iostream>
#include "utils.h"
#include "MyTree.h"
#include "mast_SW93.h"
#include "agreement_kernel.h"
// From PhylLib:
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
using namespace std;
using namespace bpp;


unsigned RF_distance(const MyTree& T1, const MyTree& T2)
{
  return TreeTools::robinsonFouldsDistance(T1, T2, false);
}

unsigned num_bipartitions(const MyTree& T1, const MyTree& T2)
{
  BipartitionList BL1(T1);
  BipartitionList BL2(T2);
  return BL1.getNumberOfBipartitions() + BL2.getNumberOfBipartitions();
}

unsigned MAST_distance(const MyTree& T1, const MyTree& T2)
{
  return T1.num_leaves() - mast(T1, T2);
}

unsigned Trip_distance(const MyTree& T1, const MyTree& T2)
{
  Hypergraph conflict_HG;
  construct_conflict_hypergraph(T1, T2, conflict_HG);
  return conflict_HG.size();
}

unsigned num_triples(const MyTree& T1)
{
  const unsigned n = T1.num_leaves();
  unsigned result = (n * (n-1)) / 2;
  return (result * (n-2)) / 3;
}

int main(int argc, char** argv){
  vector<MyTree*> trees, tmp;
  if(argc < 2){
    cout << "usage: "<< argv[0] << " <infile> [outfile|check]" << endl;
    return 1;
  }
	trees = readTrees(argv[1]);
  if(argc > 1) tmp = readTrees(argv[2]);
  trees.insert(trees.end(), tmp.begin(), tmp.end());
  // NOTE: to call mastRL, we need the following preprocessing steps:
  // Step 1: seperate a candidate t from the trees
  MyTree* t = trees.back();
  trees.pop_back();
  // Step 2: setup t's node infos (depths, clades, stids, ...) & lcas
  t->setup_node_infos(true);
  t->setup_triplets();
  t->lca_preprocess();
  //t->pretty_print();
  cout << "RF\tTrip\tRL\tnrmRF\tnrmTp\tnrmRL" << endl;
  for(MyTree* other: trees){
    other->setup_node_infos();
    other->setup_triplets();

    const unsigned d_rf = RF_distance(*t, *other);
    const unsigned d_trip = Trip_distance(*t, *other);
    const unsigned d_mast = MAST_distance(*t, *other);
    cout << d_rf << "\t"
         << d_trip << "\t"
         << d_mast << "\t"
         << fixed << setprecision(3) << (float)(d_rf)/num_bipartitions(*t, *other) << "\t"
         << fixed << setprecision(3) << (float)(d_trip)/num_triples(*t) << "\t"
         << fixed << setprecision(3) << (float)(d_mast)/t->num_leaves()
         << endl;
  }
  return 0;
}


