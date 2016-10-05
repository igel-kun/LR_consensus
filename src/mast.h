// Created by: Celine Scornavacca

#ifndef MAST_H_
#define MAST_H_

#include <cassert>
#include "MyTree.h"
#include "MultiGraph.h"

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

	#warning another preoder traversal, eliminate one!
	
	unsigned orderedIds[leavesSecondTree.size()] ; //leaf ids ordered in preorder in T2
	for(unsigned y=0;y< leavesSecondTree.size();y++){
		orderedIds[leavesSecondTree[y]->getInfos().getPreOrder()] = leavesSecondTree[y]->getInfos().getStId();
	}
	vector < MyNode * > rootsOfMiSubtrees; // the roots of the subtrees that are hanging from the centroid path of the root
	MyNode * currentNode = trees[0]->getRootNode();
	unsigned i=0;
	vector < MyNode * > ui;
	
	while(! currentNode->isLeaf()){
	    ui.push_back(currentNode); //the ui are set 
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
	     i++; //number of Mi

	}
	//int up_id =currentNode->getInfos().getStId() ;  // the stid of up
    ui.push_back(currentNode);
    


	vector<unsigned> * corrMiT2 = new vector<unsigned>[i-1];

	for(unsigned y=0;y< leavesSecondTree.size();y++){
		corrMiT2[corrMi[orderedIds[y]]].push_back(orderedIds[y]); // to comment
	}

    vector<MyTree *> Sis;
    
	for(unsigned y=0;y< rootsOfMiSubtrees.size();y++){

		vector<MyTree *> recursiveCallTrees;
		recursiveCallTrees.push_back(new MyTree(* rootsOfMiSubtrees[y])); // subtrees Mi in T1
		MyTree * Si = trees[1]->induced_subtree(corrMiT2[y]);
		recursiveCallTrees.push_back(Si);
		Sis.push_back(Si);
		mast(recursiveCallTrees, association);
	}
	
	//STEP 3
	
	//preprocess trees[1] to store, for each node y of trees[1], the root of the centroid paths to which y belongs
	
	deque < MyNode * > nodesToConsider;
	unsigned corrNodesRootCentroidPaths[trees[1]->getCorrespondanceLenghtId()];
	nodesToConsider.push_back(trees[1]->getRootNode());
	unsigned rootId = trees[1]->getRootNode()->getInfos().getStId();
	corrNodesRootCentroidPaths[rootId]=rootId;
	unsigned numberCentroidPathsT2=1;
	while(! nodesToConsider.empty()){
		MyNode * currentNode =nodesToConsider.front();
    	nodesToConsider.pop_front();
		if(currentNode->hasFather()){
			MyNode * father = currentNode->getFather();
			int idFather= father->getInfos().getStId();
			int idNode= currentNode->getInfos().getStId();
		 	if(currentNode->getInfos().getCentroidPathNumber()==father->getInfos().getCentroidPathNumber())
		 		corrNodesRootCentroidPaths[idNode]=corrNodesRootCentroidPaths[idFather]; //same root of centroid path as its father
		 	else{
		 		corrNodesRootCentroidPaths[idNode]=idNode; //the node is a root of a centroid path in T2
				numberCentroidPathsT2++;
			}
			
		}
		unsigned int numSons= currentNode->getNumberOfSons();
		for(unsigned int j = 0;j < numSons; j++)
			nodesToConsider.push_back(currentNode->getSon(j)); 
	}
	
	//creating all the muligraphs we will need 
	MultiGraph ** Gxs = new MultiGraph*[numberCentroidPathsT2];
	for(unsigned int j = 0;j < numberCentroidPathsT2; j++)
		Gxs[j]= new MultiGraph();
		
	
	//create the Rxs set, with the nodes ordered in preorder 
	nodesToConsider.push_back(trees[1]->getRootNode());
	while(! nodesToConsider.empty()){
		MyNode * currentNode =nodesToConsider.front();
    	nodesToConsider.pop_front();
    	pair< int, bool > v(currentNode->getInfos().getStId(), true);
    	auto vertex = boost::add_vertex(v, * Gxs[currentNode->getInfos().getCentroidPathNumber()]->getGraph());
    	Gxs[currentNode->getInfos().getCentroidPathNumber()]->addRxVertex(currentNode->getInfos().getStId());
		unsigned int numSons= currentNode->getNumberOfSons();
		for(unsigned int j = 0;j < numSons; j++)
			nodesToConsider.push_back(currentNode->getSon(j)); 
	}
	
	// adding vertices and edges involving up
	
	if(! trees[1]->isPreprocessed()){
  		trees[1]->setClades();
  }	
  
	MyNode * twin_up_in_T2 = trees[1]->getNodeWithStId(ui[ui.size()-1]->getInfos().getStId());
	for(unsigned int j = 0;j < numberCentroidPathsT2; j++){
		MyNode * rootOfCP = trees[1]->getNodeWithStId(Gxs[j]->getRxVertices()[0]); //this is true because of the way RxVertices are added
		int index_twin_up_in_T2_leafPreOrder = twin_up_in_T2->getInfos().getClade().first;
		pair <int, int > clade_rootOfCP = rootOfCP->getInfos().getClade();
		if (index_twin_up_in_T2_leafPreOrder >= clade_rootOfCP.first && index_twin_up_in_T2_leafPreOrder <=clade_rootOfCP.second){ // the twin of up in T2 is in L(T2x)
			pair< int, bool > v(ui[ui.size()-1]->getInfos().getStId(), false);
			auto vertex = boost::add_vertex(v, * Gxs[j]->getGraph());
    		Gxs[j]->addRxVertex(vertex);
    		MyNode * vq = trees[1]->getNodeWithStId(Gxs[j]->getRxVertices().back()); //this is true because of the way RxVertices are added
			const MyNode *lca = trees[1]->getLCA(vq,twin_up_in_T2);
			if(lca->getInfos().getCentroidPathNumber()==j){//the LCA is in Rx
				pair< int, bool > vlca(lca->getInfos().getStId(), true);
				boost::add_edge_by_label(v,vlca, EdgeMultiGraph{ 1, "white" }, * Gxs[j]->getGraph() );
				boost::add_edge_by_label(v,vlca, EdgeMultiGraph{ 1, "red" }, * Gxs[j]->getGraph() );
				boost::add_edge_by_label(v,vlca, EdgeMultiGraph{ 1, "green" }, * Gxs[j]->getGraph() );
			}
		}
	}
	
	// adding vertices and edges involving ui, with i \neq p
	
	for(unsigned int i = 0;i < Sis; i++){
	    vector <MyNode *> nodes_Sis = Sis.getNodes();
	    vector <MyNode *> nodes_T2 = trees[1].getNodes();
        for(unsigned int i = 0;i < nodes_T2.size(); i++){
            nodes_T2[i]->getInfos().setIsVisisted(false); //get isVisisted false for all node of T2
        }
        for(unsigned int i = 0;i < nodes_Sis.size(); i++){
            MyNode * z = nodes_T2[i]->getNodeWithStId(Sis[i]->getInfos().getStId());
            int cp_father_z_in_T2 = nodes_T2[i]->getNodeWithStId(Sis[i]->getFather()->getInfos().getStId())->getInfos().getCentroidPathNumber();
            while((z->getInfos().getCentroidPathNumber() != cp_father_z_in_T2) && corrNodesRootCentroidPaths[z->getInfos().getStId()]->hasFather()){
                z= corrNodesRootCentroidPaths[z->getInfos().getStId()]->getFather();
                z->getInfos().setIsVisisted(true);
            }   
        }
        
        for(unsigned int i = 0;i < nodes_T2.size(); i++){
            if(nodes_T2[i]->getInfos().getIsVisisted()){
            	pair< int, bool > l(ui[i]->getInfos().getStId(), false); //ui is to add to a Lx
            	int graph_to_add_edge_to = corrNodesRootCentroidPaths[nodes_T2[i]->getInfos().getStId()]->getInfos().getCentroidPathNumber() //get the Gx to which add ui
			    auto v_ui = boost::add_vertex(l, * Gxs[graph_to_add_edge_to]->getGraph());
			    pair< int, bool > r(nodes_T2[i]->getInfos().getStId(), true); //get the node in Rx corresponding to nodes_T2[i]
			    auto v_y = boost::vertex_by_label(l, * Gxs[graph_to_add_edge_to]->getGraph()); //y is already in Gx
                boost::edge_by_label(v_ui,v_y, EdgeMultiGraph{ -1, "white" }, * Gxs[graph_to_add_edge_to]->getGraph()); //weights not set yet
                boost::edge_by_label(v_ui,v_y, EdgeMultiGraph{ -1, "red" }, * Gxs[graph_to_add_edge_to]->getGraph()); //weights not set yet
                boost::edge_by_label(v_ui,v_y, EdgeMultiGraph{ -1, "green" }, * Gxs[graph_to_add_edge_to]->getGraph()); //weights not set yet

            }
        }
            
            
            edge_by_label
        }    	    
	}

	

	
	delete [] corrMiT2;

 }

#endif /*MAST_H_*/
