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
		
		//the next two loops are not optimal
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

void mast(vector<MyTree *> trees, map<string,int> association)
{
	map<string,int>::iterator iter;
	int maxId = trees[0]->getLeaves().size(); 
	
	for(unsigned int y=0;y< trees.size();y++){
		vector <MyNode * >  nodes= trees[y]->getNodes(); 
		  //create a vector of pointers from stIds to nodes
		trees[y]->setCorrispondanceLenghtId(nodes.size());
		for(unsigned int l=0;l< nodes.size();l++){	
		// twin leaves will have the same stId in the two trees, associated to the name via the mapping association 
		// the internal nodes will have stIds starting from maxId up
			if(nodes[l]->isLeaf()){	  
				iter= association.find( nodes[l]->getName());	
				nodes[l]->getInfos().setStId(iter->second);
				//set the pointer from stId to node
				trees[y]->setNodeWithStId(nodes[l], iter->second);
			}
			else{
				nodes[l]->getInfos().setStId(maxId);
				//set the pointer from stId to node
				trees[y]->setNodeWithStId(nodes[l],maxId);
				maxId++;
			}
		}
		trees[y]->setDepthAndNumberOfDescendents();	

		//STEP 1
	
		trees[y]->getCentroidDecompostion();
	}
	
	//STEP 2
	
	trees[1]->setPreOrder();
	vector < MyNode * > tempLeaves = trees[1]->getLeaves();
	int orderedIds[tempLeaves.size()] ; //leaf ids ordered in preorder in T2
	for(unsigned int y=0;y< tempLeaves.size();y++){
		orderedIds[tempLeaves[y]->getInfos().getPreorder()] = tempLeaves[y]->getInfos().getStId();
	}
	vector < MyNode * > rootsOfMiSubtrees; // the roots of the subtrees that are hanging from the centroid path of the root
	MyNode * currentNode = trees[0]->getRootNode();
	while(! currentNode->isLeaf()){
	     if(currentNode->getSon(0)->getInfos().getCentroidPathNumber()==0){
	     	currentNode = currentNode->getSon(0);
	     	rootsOfMiSubtrees.push_back( currentNode->getSon(1));
	     }	
	     else{
	     	currentNode = currentNode->getSon(1);
	     	rootsOfMiSubtrees.push_back( currentNode->getSon(0));
	     }
	}  
	  
	for(unsigned int y=0;y< rootsOfMiSubtrees.size();y++){
		vector<MyTree *> recoursiveCallTrees;
		recoursiveCallTrees.push_back(new MyTree(* rootsOfMiSubtrees[y])); // subtrees Mi in T1
	}
 
	     
	
	

	
 //preprocessings, LCA prep here;
}

#endif /*MAST_H_*/
