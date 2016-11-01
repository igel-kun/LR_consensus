

#pragma once

#include <fstream>
#include <vector>
#include <string>
#include "MyTree.h"
#include "candidate_tree.h"

using namespace std;

/* this function reads a list of trees written in a file (path) in a newich format, separed by semicolons and returns a list of MyTree*/
void readTrees(CandidateTree*& candidate, vector<MyTree*>& other_trees, const string & path) throw (Exception) {
  assert(candidate == NULL); // the candidate should always be NULL before
  // Checking the existence of specified file
  ifstream file(path.c_str(), ios::in);
  if (! file) { throw IOException ("\nError reading file.\n"); }
    
  string description;
  while(file.good()) {
    const string temp = FileTools::getNextLine(file);
    if(temp.size() != 0){
      const string::size_type index = temp.find_first_of(';');
      if(index == string::npos) throw Exception("readTrees(). Bad format: no semi-colon found.");
      description += temp.substr(0, index + 1);
      // put the first tree as CandidateTree
      if(candidate)
        other_trees.push_back(new MyTree(description));
      else
        candidate = new CandidateTree(description);
      description = temp.substr(index);
    }
  }
  file.close();
};

