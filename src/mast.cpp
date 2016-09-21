


//
// File: mast.cpp
// Created by: Celine Scornavacca, Mathias Weller
// Created on: Sep Wed 21 11:27 2016
//

// From the STL:
#include <iostream>

using namespace std;

// From PhylLib:
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeExceptions.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

// From NumCalc:
#include <Bpp/Numeric/VectorTools.h>

// From Utils:
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Utils/AttributesTools.h>

using namespace bpp;

#include "NodeInfos.h"
#include "MyTree.h"
#include "agreement_kernel.h"

//typedef NodeTemplate<NodeInfos> MyNode;
//typedef TreeTemplate<MyNode> MyTree;

#include "utils.h"


const string help_text = (string)"" +
                         "____________________________________________________________________________\n" +
                         "This program takes as input a path to a list of phylogenetic trees in newick\n" +
                         "format and an integer d and outputs a MAST-RL for the input, if any exists. \n" +
                         "input.list.file            | [path] toward multi-trees file (newick)        \n" +
                         "max.distance               | [integer] toward multi-trees file (newick)     \n" +
                         "output.tree.file           | file where to write the MAST-RL tree           \n" +
                         "___________________________|________________________________________________\n";

void quit(const string& message, const int return_code)
{
  if(return_code != 0)
    cerr << message << endl;
  else
    cout << message << endl;
  exit(return_code);
}

void parse_params(const int args, char** const argv, map<string, string>& params)
{
  // Get the parameters from command line:
  params = AttributesTools::getAttributesMap(AttributesTools::getVector(args, argv), "=");
  
  // Look for a specified file with parameters:
  if(params.find("param") != params.end()) {
    const string filename = params.at("param");
    if(FileTools::fileExists(filename)) {
      const map<string, string> file_params = AttributesTools::getAttributesMapFromFile(filename, "=");
      // Update attributes with ones passed to command line:
      AttributesTools::actualizeAttributesMap(params, file_params);
    } else quit("Parameter file not found.", -1);
  }
}

int main(int args, char ** argv){
	  
  cout << "******************************************************************" << endl;
  cout << "*               MAST-RL, version 0.1.0 *" << endl;
  cout << "* Authors: C. Scornavacca, Mathias Weller   Created     21/09/16 *" << endl;
  cout << "*                                           Last Modif. 21/09/16 *" << endl;
  cout << "******************************************************************" << endl << endl;

	if(args == 1) quit(help_text, 0);
  
		ApplicationTools::startTimer();
		
		cout << "Parsing options:" << endl;
    map<string, string> params;
    parse_params(args, argv, params);
		  
		const string listPath = ApplicationTools::getAFilePath("input.list.file", params);
		ApplicationTools::displayResult("Input list file", listPath);
		if(listPath == "none") throw Exception("You must provide an input tree list file.");
	    int d =  ApplicationTools::getDoubleParameter("max.distance", params, 1);

		       
		string outputPath = ApplicationTools::getAFilePath("output.tree.file", params, false, false, "",true,"MAST-RL.txt");
		ApplicationTools::displayResult("Output file", outputPath);

	  	vector<MyTree *> trees;
	      
	  	//Reading trees to root
		trees = readTrees(listPath);
		
	  	vector < string> leaves[trees.size() +1];
	  		  	
	  	for(unsigned int y=0;y< trees.size();y++){
			leaves[y]= (* trees[y]).getLeavesNames();  
			sort(leaves[y].begin(), leaves[y].end());
		}
		
	  	vector <string> AllLeaves = allLeaves(leaves, trees.size());
	  	
	  	map<string,int> association;
	  	map<string,int>::iterator iter;
	  	
	  	for (unsigned int i =0;i<  AllLeaves.size(); i++){
	  		association.insert( make_pair(  AllLeaves[i], i));
		}  	
	  	vector < int> idLeaves[trees.size() +1];
	  			
	  	for(unsigned int y=0;y< trees.size();y++){
	  		vector <MyNode * >  leavesForId= (* trees[y]).getLeaves();  
			for(unsigned int l=0;l< leavesForId.size();l++){		  
				iter= association.find((* leavesForId[l]).getName());	
				((* leavesForId[l]).getInfos()).setStid(iter->second);
				(idLeaves[y]).push_back(iter->second);
			}	
			sort(idLeaves[y].begin(),idLeaves[y].end());
	  	}
	  	
		
		if(trees.size()>1){
		    //newick.write(* trees[t], outputPath, false);	
			for(unsigned int i = 0; i < trees.size(); i++) delete trees[i];
					
		}
		else{
			cout << "In file " << listPath << " there is only a tree...\n";
		}	
	
	return 0;
};	
	
		
		


