// Created by: Celine Scornavacca

#ifndef MAST_H_
#define MAST_H_

#include <cassert>
#include "MyTree.h"

using namespace bpp;



void mast(vector<MyTree *> trees, map<string,unsigned> association)
{
	map<string,unsigned>::iterator iter;
	unsigned maxId = trees[0]->getLeaves().size();

	for(unsigned y=0;y< trees.size();y++){
		vector <MyNode * >  nodes= trees[y]->getNodes();
		  //create a vector of pounsigneders from stIds to nodes
		trees[y]->setCorrespondanceLenghtId(nodes.size());
		for(unsigned l=0;l< nodes.size();l++){
		// twin leaves will have the same stId in the two trees, associated to the name via the mapping association
		// the unsignedernal nodes will have stIds starting from maxId up
			if(nodes[l]->isLeaf()){
				iter= association.find( nodes[l]->getName());
				nodes[l]->getInfos().setStId(iter->second);
				//set the pounsigneder from stId to node
				trees[y]->setNodeWithStId(nodes[l], iter->second);
			}
			else{
				nodes[l]->getInfos().setStId(maxId);
				//set the pounsigneder from stId to node
				trees[y]->setNodeWithStId(nodes[l],maxId);
				maxId++;
			}
		}
		trees[y]->setDepthAndNumberOfDescendents();

		//STEP 1

		trees[y]->getCentroidDecompostion();
		trees[y]->lca_preprocess();
	}

	//STEP 2

	trees[1]->setPreOrder();
	vector < MyNode * > leavesSecondTree = trees[1]->getLeaves();
	unsigned corrMi[leavesSecondTree.size()];

	unsigned orderedIds[leavesSecondTree.size()] ; //leaf ids ordered in preorder in T2
	for(unsigned y=0;y< leavesSecondTree.size();y++){
		orderedIds[leavesSecondTree[y]->getInfos().getPreOrder()] = leavesSecondTree[y]->getInfos().getStId();
	}
	vector < MyNode * > rootsOfMiSubtrees; // the roots of the subtrees that are hanging from the centroid path of the root
	MyNode * currentNode = trees[0]->getRootNode();
	unsigned i=0;
	while(! currentNode->isLeaf()){
	     if(currentNode->getSon(0)->getInfos().getCentroidPathNumber()==0){
	     	currentNode = currentNode->getSon(0);
	     	rootsOfMiSubtrees.push_back( currentNode->getSon(1));
	     }
	     else{
	     	currentNode = currentNode->getSon(1);
	     	rootsOfMiSubtrees.push_back( currentNode->getSon(0));
	     }
	     vector < MyNode * > leavesNode = TreeTemplateTools::getLeaves(* rootsOfMiSubtrees[rootsOfMiSubtrees.size()]);
		 for(unsigned y=0;y< leavesNode.size();y++){
			corrMi[leavesNode[y]->getInfos().getStId()]	=i ; //this leaf is in subtree Mi
		 }
	     i++;

	}


	vector<unsigned> * corrMiT2 = new vector<unsigned>[i-1];

	for(unsigned y=0;y< leavesSecondTree.size();y++){
		corrMiT2[corrMi[orderedIds[y]]].push_back(orderedIds[y]); // to comment
	}

	for(unsigned y=0;y< rootsOfMiSubtrees.size();y++){

		vector<MyTree *> recursiveCallTrees;
		recursiveCallTrees.push_back(new MyTree(* rootsOfMiSubtrees[y])); // subtrees Mi in T1
		recursiveCallTrees.push_back(trees[1]->induced_subtree(corrMiT2[y]));
		mast(recursiveCallTrees, association);
	}
	delete [] corrMiT2;

 }

#endif /*MAST_H_*/
