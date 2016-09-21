// Created by: Celine Scornavacca

#ifndef MYTREE_H_
#define MYTREE_H_


#include "NodeInfos.h"

using namespace bpp;

#include <Bpp/Phyl/TreeTemplateTools.h>
typedef NodeTemplate<NodeInfos> MyNode;

class MyTree: public TreeTemplate<MyNode> {
public:
	MyNode** correspondance;
	MyNode** correspondanceId;

	int dimCorr;
	int dimCorrId;
	int dimTS;
 	MyTree(int dim): TreeTemplate<MyNode>(){}

 	MyTree(MyNode& root): TreeTemplate<MyNode>(& root), correspondance(NULL), correspondanceId(NULL){}
 	MyTree(MyNode& root, int dim): TreeTemplate<MyNode>(& root), correspondance(new MyNode*[dim]), correspondanceId(new MyNode*[dim]){} 
 	
 	virtual ~MyTree()
  {
 		if(correspondance)
 			delete [] (correspondance);
 		if(correspondanceId)
 			delete [] (correspondanceId);
 	} 
 	
  void setCorrispondanceNull()
  {
    for(int i = 0; i < dimCorr; ++i)
 		  correspondance[i]= NULL; 	
 	}
 	
 	void setCorrispondanceLenght(int dim){
 		correspondance= new MyNode*[dim];
 		for(int i = 0; i < dim; i++)
 			correspondance[i]= NULL;
 		dimCorr=dim;
 	}
 	
 	 void setCorrispondanceLenghtId(int dim){
 		correspondanceId= new MyNode*[dim];
 		for(int i = 0; i < dim; i++)
 			correspondanceId[i]= NULL;
 		dimCorr=dim;
 	}
 	
 	
 	int getCorrispondanceLenght(){
		return dimCorr;
 	}
 	
	void deleteCorrispondanceLenght(){
 		delete [] (correspondance);
 	}
 	
 	void resetCorrispondanceLenght(int dim){
 		for(int i = 0; i < dim; i++)
 			correspondance[i]= NULL;
 	}
 	
 	 void resetCorrispondanceLenghtId(int dim){
 		for(int i = 0; i < dim; i++)
 			correspondanceId[i]= NULL;
 	}

	 MyNode * getNodeWithPostOrder(int postOrder){
	 	return correspondance[postOrder];
	 }	
	 
	 void setNodeWithPostOrder(MyNode * node, int postOrder){
	 	correspondance[postOrder] = node;
	 }
	 

	 void setNodeWithStId(MyNode * node, int stId){
	 	correspondanceId[stId] = node;
	 }
	 
	 	 
	 MyNode * getNodeWithStId(int stId){
	 	return correspondanceId[stId];
	 }	
	 
	 
	 vector<string> getNames(MyNode & node){
	   	vector<string> stId;
	  	if(node.isLeaf()){
	     	stId.push_back(node.getName());
	  	 }
	  	 for(unsigned int i = 0; i < node.getNumberOfSons(); i++){
	     	vector<string> subStId = getNames(* node.getSon(i));
	    	for(unsigned int j = 0; j < subStId.size(); j++) stId.push_back(subStId[j]);
	   	}
	   	return stId;   
	 }
	
	vector<int> getLeavesStId(MyNode & node){
	   	vector<int> stId;
	  	if(node.isLeaf()){
	     	stId.push_back(node.getInfos().getStid());
	  	 }
	  	 for(unsigned int i = 0; i < node.getNumberOfSons(); i++){
	     	vector<int> subStId = getLeavesStId(* node.getSon(i));
	    	for(unsigned int j = 0; j < subStId.size(); j++) stId.push_back(subStId[j]);
	   	}
	   	return stId;   
	 }
	 
	 
};	


#endif /*MYTREE_H_*/

