// Created by: Celine Scornavacca

#ifndef MYTREE_H_
#define MYTREE_H_


#include "NodeInfos.h"

using namespace bpp;

#include <Bpp/Phyl/TreeTemplateTools.h>
typedef NodeTemplate<NodeInfos> MyNode;

class MyTree: public TreeTemplate<MyNode> {
public:

	MyNode ** correspondanceId ;
	int dimCorrId;

 	MyTree(int dim): TreeTemplate<MyNode>(){}

 	MyTree(MyNode & root): TreeTemplate<MyNode>(& root), correspondanceId(NULL){}
 	MyTree(MyNode & root,int dim): TreeTemplate<MyNode>(& root), correspondanceId(new MyNode*[dim]){} 
 	
 	virtual ~MyTree(){
 		if(correspondanceId)
 			delete [] (correspondanceId);
 	} 
 	
 	 void setCorrispondanceLenghtId(int dim){
 		correspondanceId= new MyNode*[dim];
 		for(int i = 0; i < dim; i++)
 			correspondanceId[i]= NULL;
 		dimCorrId=dim;
 	}
 		
 	 void resetCorrispondanceLenghtId(int dim){
 		for(int i = 0; i < dim; i++)
 			correspondanceId[i]= NULL;
 	}

	 void setNodeWithStId(MyNode * node, int stId){
	 	correspondanceId[stId] = node;
	 }
	 
	 	 
	 MyNode * getNodeWithStId(int stId){
	 	return correspondanceId[stId];
	 }	

	 void setDepth(MyNode * node){
	    setDepth();
	 }	 	 
	 
	 void setNumberOfDescendents(MyNode * node){
	    setNumberOfDescendents();
	 }
	 
	 void setDepth(){
	    setDepth(this->getRootNode());
	 }	
	  	 
	 void setNumberOfDescendents(){
	    setNumberOfDescendents(this->getRootNode());
	 }
	 	 
	 void setDepthAndNumberOfDescendents(){
	    setDepth();
	    setNumberOfDescendents();
	 }
	 

	 
	//not used for now
 
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
