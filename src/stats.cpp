

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




struct Statistics 
{
  const vector<double>& data;
  
  Statistics(const vector<double>& _data): data(_data)
  {}
  
  double getMean() const
  {
    double sum = 0.0;
    for(double a : data)
      sum += a;
    return sum/data.size();
  }
  
  double getVariance() const
  {
    double mean = getMean();
    double temp = 0;
    for(double a :data)
      temp += (a-mean)*(a-mean);
    return temp/data.size();
  }

  double getStdDev() const { return std::sqrt(getVariance()); }

};




int main(int argc, char** argv)
{
  vector<MyTree*> trees, tmp;
  if(argc < 2){
    cout << "usage: "<< argv[0] << " <infile> [outfile|check]" << endl;
    return 1;
  }
  cout << "reading trees from "<<argv[1]<<endl;
	trees = readTrees(argv[1]);
  if(argc > 2) {
    cout << "reading trees from "<<argv[2]<<endl;
    tmp = readTrees(argv[2]);
    trees.insert(trees.end(), tmp.begin(), tmp.end());
  }
  // NOTE: to call mastRL, we need the following preprocessing steps:
  // Step 1: seperate a candidate t from the trees
  cout << "mean\tstddiv"<<endl;
  for(MyTree* t: trees){
    t->setup_node_infos();
    t->setup_triplets();

    vector<double> supports;
    
    vector<MyNode*> nodes = t->getNodes();
    for(const MyNode* n: nodes){
      if(n->hasBranchProperty(TreeTools::BOOTSTRAP)) supports.push_back(dynamic_cast<const Number<double>*>(n->getBranchProperty(TreeTools::BOOTSTRAP))->getValue());
    }

/*      if(!n->isLeaf()) 
        if(n->hasName()) {
          cout << "using "<<n->getName()<<endl;
          supports.push_back(std::atof(n->getName().c_str()));
        }
*/
/*    for(MyNode& c: get_children(*t->getRootNode())){
      vector<double> new_supports = TreeTemplateTools::getBranchLengths(c);
      std::copy(new_supports.begin(), new_supports.end(), std::back_inserter(supports));
    }
*/
    Statistics stat(supports);

    cout << fixed << setprecision(2) << stat.getMean() << "\t" << stat.getStdDev()<<endl;
  }
  return 0;
}


