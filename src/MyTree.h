// Created by: Celine Scornavacca

#ifndef MYTREE_H_
#define MYTREE_H_


#include "rmq/lca.hpp"
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeExceptions.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

#include "NodeInfos.h"

using namespace bpp;

typedef NodeTemplate<NodeInfos> MyNode;

class MyTree: public TreeTemplate<MyNode> {
public:

	MyNode** correspondanceId;
	unsigned dimCorrId;
	//vector< unsigned > centroidPaths ;

 	MyTree(): TreeTemplate<MyNode>(){}
 	MyTree(MyNode& root): TreeTemplate<MyNode>(& root), correspondanceId(NULL){}
 	MyTree(MyNode& root, unsigned dim): TreeTemplate<MyNode>(& root), correspondanceId(new MyNode*[dim]){}

 	virtual ~MyTree(){
    if(lca_oracle) delete lca_oracle;
 		if(correspondanceId) delete [] (correspondanceId);
 	}

  void setCorrispondanceLenghtId(unsigned dim){
 		correspondanceId= new MyNode*[dim];
 		for(unsigned i = 0; i < dim; i++)
 			correspondanceId[i]= NULL;
 		dimCorrId=dim;
 	}

  void resetCorrispondanceLenghtId(unsigned dim){
 		for(unsigned i = 0; i < dim; i++)
 			correspondanceId[i]= NULL;
 	}

	void setNodeWithStId(MyNode * node, unsigned stId){
	 	correspondanceId[stId] = node;
	}

	MyNode * getNodeWithStId(unsigned stId){
	 	return correspondanceId[stId];
	}

	void setDepth(MyNode * node){
	 	if(TreeTemplateTools::isRoot(* node))
	 		node->getInfos().setDepth(0); //root has depth 0
	 	else
	 		node->getInfos().setDepth(node->getFather()->getInfos().getDepth() +1); //otherwise depth of father +1

    for(unsigned i = 0; i < node->getNumberOfSons(); i++)
			setDepth(node->getSon(i)); // recursive calls for sons
	}

	void setNumberOfDescendents(MyNode * node){
    for(unsigned i = 0; i < node->getNumberOfSons(); i++)
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
	void getCentroidDecompostion(MyNode * node, unsigned maxPartition){
    node->getInfos().setCentroidPathNumber(maxPartition);
    unsigned maxDesc=0;
    unsigned chosenChild=0;
    for(unsigned i = 0; i < node->getNumberOfSons(); i++){
      unsigned numbDesc = node->getSon(i)->getInfos().getNumberOfDescendents();
      if(numbDesc>maxDesc)
        chosenChild=i;
    }
    for(unsigned i = 0; i < node->getNumberOfSons(); i++){
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

	void setPreOrder(MyNode * node, unsigned & maxPreorder){
    for(unsigned i = 0; i < node->getNumberOfSons(); i++){
      setPreOrder(node->getSon(i),maxPreorder);
		}
		if(node->isLeaf()){
	 		node->getInfos().setPreOrder(++maxPreorder);
	 	}
  }

  void setPreOrder(){
    unsigned zero =0;
	 	setPreOrder(this->getRootNode(),zero);
  }

	//not used for now
	 vector<string> getNames(MyNode & node){
	   	vector<string> stId;
	  	if(node.isLeaf()){
	     	stId.push_back(node.getName());
	  	 }
	  	 for(unsigned i = 0; i < node.getNumberOfSons(); i++){
	     	vector<string> subStId = getNames(* node.getSon(i));
	    	for(unsigned j = 0; j < subStId.size(); j++) stId.push_back(subStId[j]);
	   	}
	   	return stId;
	 }

	vector<unsigned> getLeavesStId(MyNode & node){
	   	vector<unsigned> stId;
	  	if(node.isLeaf()){
	     	stId.push_back(node.getInfos().getStId());
	  	 }
	  	 for(unsigned i = 0; i < node.getNumberOfSons(); i++){
	     	vector<unsigned> subStId = getLeavesStId(* node.getSon(i));
	    	for(unsigned j = 0; j < subStId.size(); j++) stId.push_back(subStId[j]);
	   	}
	   	return stId;
	 }



  //======================== LCA code =========================
protected:
  //! access node stId's
  struct access_stId {
    unsigned operator[](const MyNode& u) const { return u.getInfos().getStId(); }
  };
  struct access_Id {
    unsigned operator[](const MyNode& u) const { return u.getId(); }
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

  typedef lca<MyNode, unsigned, const access_Id, const access_children> LCA_Oracle;

  LCA_Oracle* lca_oracle = NULL;

public:
  void lca_preprocess()
  {
    lca_oracle = new LCA_Oracle(*getRootNode(), dimCorrId, access_Id(), access_children());
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
