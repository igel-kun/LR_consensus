// Created by: Celine Scornavacca

#ifndef MYTREE_H_
#define MYTREE_H_


#include "NodeInfos.h"
#include "rmq/lca.hpp"

using namespace bpp;

#include <Bpp/Phyl/TreeTemplateTools.h>
typedef NodeTemplate<NodeInfos> MyNode;

class MyTree: public TreeTemplate<MyNode> {
public:

	MyNode ** correspondanceId ;
	int dimCorrId;
	//vector< int > centroidPaths ;

 	MyTree(): TreeTemplate<MyNode>(){}

 	MyTree(MyNode & root): TreeTemplate<MyNode>(& root), correspondanceId(NULL){}
 	MyTree(MyNode & root,int dim): TreeTemplate<MyNode>(& root), correspondanceId(new MyNode*[dim]){} 
 	
 	virtual ~MyTree(){
    if(lca_oracle) delete lca_oracle;
 		if(correspondanceId) delete [] (correspondanceId);
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
	 	if(TreeTemplateTools::isRoot(* node))
	 		node->getInfos().setDepth(0); //root has depth 0
	 	else
	 		node->getInfos().setDepth(node->getFather()->getInfos().getDepth() +1); //otherwise depth of father +1
		
    for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
			setDepth(node->getSon(i)); // recursive calls for sons
	 }	 	 
	 
	 void setNumberOfDescendents(MyNode * node){
     for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
       setNumberOfDescendents(node->getSon(i)); // recursive calls for sons

     if(!node->isLeaf()){
       unsigned tempNumberOfDescendent=0;
       for(unsigned i = 0; i < node->getNumberOfSons(); i++)
         tempNumberOfDescendent += (node->getSon(i)->getInfos().getNumberOfDescendents() +1); //adding up all the descendants of sons + sons 
       node->getInfos().setNumberOfDescendents(tempNumberOfDescendent);	
     } else node->getInfos().setNumberOfDescendents(0); //leaves have no descendants 
	 }
	 
	 	 
	 void setDepthAndNumberOfDescendents(){
	    setDepth(getRootNode());
	    setNumberOfDescendents(getRootNode());
	 }
	 

	 //construct the centroid decomposition of a tree
	 
	 void getCentroidDecompostion(MyNode * node, int maxPartition){
     node->getInfos().setCentroidPathNumber(maxPartition);
     int maxDesc=0;
     int chosenChild=0;
     for(unsigned int i = 0; i < node->getNumberOfSons(); i++){
       int numbDesc = node->getSon(i)->getInfos().getNumberOfDescendents() ;
       if(numbDesc>maxDesc)
         chosenChild=i;	     	
     }
     for(unsigned int i = 0; i < node->getNumberOfSons(); i++){
       if(i== chosenChild){
         //the child with the higher number of sons stay in the same centroid path than the father (ties broken arbitrarly)
         getCentroidDecompostion(node->getSon(i),maxPartition);
       } else 
         getCentroidDecompostion(node->getSon(i),++maxPartition); // the other start new centroid paths
     }
   }
	 
	 	 
	 void getCentroidDecompostion(){
	 	getCentroidDecompostion(getRootNode(),0);  
	 }
	 
	 
	 void setPreOrder(MyNode * node, int & maxPreorder){	    
	 	for(unsigned int i = 0; i < node->getNumberOfSons(); i++){
			setPreOrder(node->getSon(i),maxPreorder); 
		}	
		if(node->isLeaf()){
	 		node->getInfos().setPreOrder(++maxPreorder); 
	 	}
	 }
	 
	 void setPreOrder(){
	 	int zero =0;
	 	setPreOrder(this->getRootNode(),zero);  
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
	     	stId.push_back(node.getInfos().getStId());
	  	 }
	  	 for(unsigned int i = 0; i < node.getNumberOfSons(); i++){
	     	vector<int> subStId = getLeavesStId(* node.getSon(i));
	    	for(unsigned int j = 0; j < subStId.size(); j++) stId.push_back(subStId[j]);
	   	}
	   	return stId;   
	 }
	 


  //======================== LCA code =========================
protected:
  //! access node stId's
  struct access_stId {
    int operator[](const MyNode& u) const { return u.getInfos().getStId(); }
  };
  // access node children
  struct access_children {
    // damn you Bio++ for not giving me access to the ::sons_ vector!!!
    // for now, I'll make my own copy of the sons vector...
    const vector<const MyNode*> operator[](const MyNode& u) const {
      vector<const MyNode*> result;
      const unsigned num_children = u.getNumberOfSons();
      result.reserve(num_children);
      for(unsigned i = 0; i < num_children; ++i)
        result.push_back(u[i]);
      return result;
    }
  };

  typedef lca<MyNode, int, const access_stId, const access_children> LCA_Oracle;

  LCA_Oracle* lca_oracle = NULL;

public:
  void lca_preprocess()
  {
    lca_oracle = new LCA_Oracle(*getRootNode(), dimCorrId, access_stId(), access_children());
  }

  const MyNode* getLCA(const MyNode& u, const MyNode& v) const
  {
    return lca_oracle->query(u, v);
  }
  const MyNode* getLCA(const MyNode * u, const MyNode * v) const
  {
    return lca_oracle->query(*u, *v);
  }
   
};	


#endif /*MYTREE_H_*/
