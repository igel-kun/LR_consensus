

using namespace std;


//
// File: mast.cpp
// Created by: Celine Scornavacca, Mathias Weller
// Created on: Sep Wed 21 11:27 2016
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

// From the STL:
#include <iostream>

using namespace std;

// From PhylLib:
#include <Bpp/Phyl/Io/IOTree.h>
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
	    int d =  ApplicationTools::getDoubleParameter("bootstrap.treshold", params, 0);

		       
		string outputPath = ApplicationTools::getAFilePath("output.tree.file", params, true, false);
		ApplicationTools::displayResult("Output file", outputPath);
		if(outputPath == "none") throw Exception("You must provide an output file.");
		      

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
	
		
		


