


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

//typedef NodeTemplate<NodeInfos> MyNode;
//typedef TreeTemplate<MyNode> MyTree;

#include "utils.h"



void help()
{
  ApplicationTools::displayMessage("__________________________________________________________________________" );
  ApplicationTools::displayMessage("This program takes as input a path to a list of phylogenetic trees in ");
  ApplicationTools::displayMessage("newick format and an integer d and outputs a MAST-RL for the input, if any exists.");
  ApplicationTools::displayMessage("input.list.file            | [path] toward multi-trees file (newick)      " );
  ApplicationTools::displayMessage("max.distance               | [integer] toward multi-trees file (newick)      " );
  ApplicationTools::displayMessage("output.tree.file           | file where to write the MAST-RL tree" );
  ApplicationTools::displayMessage("___________________________|___________________________________________" );
}


int main(int args, char ** argv){
	  
	cout << "******************************************************************" << endl;
  	cout << "*               MAST-RL, version 0.1.0 *" << endl;
  	cout << "* Authors: C. Scornavacca, Mathias Weller   Created     21/09/16 *" << endl;
  	cout << "*                                           Last Modif. 21/09/16 *" << endl;
  	cout << "******************************************************************" << endl;
  	cout << endl;

	if(args == 1)
  	{
    	help();
    	exit(0);
  	}
  
	try {
	  
		ApplicationTools::startTimer();
		
		cout << "Parsing options:" << endl;
		  
		// Get the parameters from command line:
		map<string, string> cmdParams = AttributesTools::getAttributesMap(AttributesTools::getVector(args, argv), "=");
		
		// Look for a specified file with parameters:
		map<string, string> params;
		if(cmdParams.find("param") != cmdParams.end())
		{
			string file = cmdParams["param"];
			if(!FileTools::fileExists(file))
			{
				cerr << "Parameter file not found." << endl;
			    exit(-1);
			}
			else
			{
				params = AttributesTools::getAttributesMapFromFile(file, "=");
		  	    // Actualize attributes with ones passed to command line:
		  	  AttributesTools::actualizeAttributesMap(params, cmdParams);
			}
		}
		else
		{
			params = cmdParams;
		}
		
		Newick newick;
		string listPath = ApplicationTools::getAFilePath("input.list.file", params);
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
	
	
	}
catch (std::exception&  e){
    	cout << e.what() << endl;
    	exit(-1);
	}
	return (0);
};	
	
		
		


