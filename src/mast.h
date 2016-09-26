// Created by: Celine Scornavacca

#ifndef MAST_H_
#define MAST_H_


#include "MyTree.h"

using namespace bpp;

MyTree * restrictTree(MyTree * tree, vector<int> idsLeaves){
	vector<int> ids; //used later to find the edge set
	
	//defining the node set
	for(unsigned int y=0;y< idsLeaves.size()-1;y++){
		ids.push_back(idsLeaves[y]);
		int idNode = 0;//tree->getLCA(idsLeaves[y],idsLeaves[y+1]); //internal nodes  
		ids.push_back(idNode);
	}
	ids.push_back(ids.size());
	
	MyTree * restrictedTree ;
	restrictedTree->setCorrispondanceLenghtId(ids.size());

	//defining the edge set

	for(unsigned int y=0;y< ids.size();y++){

		MyNode * newNode =new MyNode();
		newNode->getInfos().setStId(ids[y]);
		restrictedTree->setNodeWithStId(newNode,ids[y]);	
		if(tree->getNodeWithStId(ids[y])->isLeaf())
			newNode->setName(tree->getNodeWithStId(ids[y])->getName());
	
	
		int vright = -1;
		int vleft = -1;
		int depthNode = tree->getNodeWithStId(ids[y])->getInfos().getDepth();
		int depthNodeR;
		int depthNodeL;
		
		for(unsigned int l=y-1;l>=0 ;l--){
			depthNodeL = tree->getNodeWithStId(ids[l])->getInfos().getDepth();
			if(depthNodeL < depthNode){
				vleft=ids[l]; //first node on the left of idsForEdges[y] with lower depth
				break;
			}

		}
		
		
		for(unsigned int r=y+1;r< ids.size();r++){
			depthNodeR = tree->getNodeWithStId(ids[r])->getInfos().getDepth();
			if(depthNodeR < depthNode){
				vright=ids[r]; //first node on the right of idsForEdges[y] with lower depth
				break;
			}
		}
		if(vleft==-1 && vright!=-1){
			restrictedTree->getNodeWithStId(vright)->addSon(newNode);
		}
		else if(vleft!=-1 && vright==-1){
			restrictedTree->getNodeWithStId(vleft)->addSon(newNode);

		}
		else if(vleft==-1 && vright==-1){
			restrictedTree->setRootNode(newNode);
			restrictedTree->resetNodesId();
		}
		
		else{
			if(depthNodeR< depthNodeL)
				restrictedTree->getNodeWithStId(vleft)->addSon(newNode);
			else
				restrictedTree->getNodeWithStId(vright)->addSon(newNode);
		}
	}
	return restrictedTree;

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
