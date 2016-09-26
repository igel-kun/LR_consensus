// Created by: Celine Scornavacca

#ifndef MAST_H_
#define MAST_H_


#include "MyTree.h"

using namespace bpp;

MyTree * restrictTree(MyTree * tree, vector<int> ids){
	MyTree * restrictedTree ;
	restrictedTree->setCorrispondanceLenghtId(ids.size());
	vector<int> idsForEdges; //used later to find the edge set
	//defining the node set
	for(unsigned int y=0;y< ids.size();y++){ //leaves
		MyNode * newNode =new MyNode();
		newNode->getInfos().setStId(ids[y]);
		newNode->setName(tree->getNodeWithStId(idNode).getName());
		setNodeWithStId(newNode,idNode);	
	}
	for(unsigned int y=0;y< ids.size()-1;y++){//internal nodes  
		idsForEdges.push_back(ids[y]);
		int idNode = tree->getLCA(ids[y],ids[y+1]);
		idsForEdges.push_back(idNode);
		MyNode * newNode =new MyNode();
		newNode->getInfos().setStId(idNode);
		setNodeWithStId(newNode,idNode);	
	}
	idsForEdges.push_back(ids.size());
	
	//defining the edge set

	for(unsigned int y=0;y< idsForEdges.size()-1;y++){
		int vright = -1;
		int vleft = -1;

	}


}

void mast(MyTree * tree1, MyTree * tree2, bool firstCall)
{
	tree1->setDepthAndNumberOfDescendents();	
	tree2->setDepthAndNumberOfDescendents();
	if(firstCall){	
		tree1->getCentroidDecompostion();
		tree2->getCentroidDecompostion();
	}
 //preprocessings, LCA prep here;
}

#endif /*MAST_H_*/
