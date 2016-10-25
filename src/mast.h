// Created by: Celine Scornavacca

#ifndef MAST_H_
#define MAST_H_

#include <cassert>
#include "MyTree.h"
#include "MultiGraph.h"

using namespace bpp;

void agreementMatching(MultiGraph & Gx){
    #warning TODO
};


unsigned* mast(MyTree* T1, MyNode* root_of_T2);
unsigned* mast(MyNode* root_of_T1, MyTree* T2);


//! this computes the mast between T1 and each subtree of T2 in an array, mapping StId of a vertex v of T2 to the size of a mast between T1 & v
//NOTE: for StIds that do not occur in T2, the mapping is undefined
//NOTE: the caller should free the allocated space!
unsigned* mast(MyTree* T1, MyTree* T2)
{
	// ============= STEP 1 ==================
	//preprocessing for T1 and T2
  T2->sync_leaf_stids(*T1);
	for(MyTree* t: {T1, T2}){
		t->setDepthAndNumberOfDescendents();
		t->getCentroidDecompostion();
		t->lca_preprocess();
#warning another preoder traversal, eliminate one!
	  t->setLeafPreOrder();
	}

	// ============= STEP 2 ==================
	vector<MyNode*> leavesT2 = T2->getLeaves();
	StId corrMi[leavesT2.size()];
	StId orderedIds[leavesT2.size()] ; //leaf ids ordered in preorder in T2

	for(const MyNode* leaf: leavesT2)
		orderedIds[leaf->getInfos().getPreOrder()] = stid(leaf);
	
	vector<MyNode*> Mi_roots; // the roots of the subtrees that are hanging from the centroid path of the root
	MyNode* currentNode = T1->getRootNode();
	unsigned p = 0;
	vector<MyNode*> ui;
	
	while(!currentNode->isLeaf()){
    ui.push_back(currentNode); //the ui are set 
	  const bool Son0next = (currentNode->getSon(0)->getInfos().getCentroidPathNumber() == 0);
    currentNode = currentNode->getSon( Son0next ? 0 : 1);
    Mi_roots.push_back( currentNode->getSon( Son0next ? 1 : 0));

    // save the StIds of the leaves of Mi in corrMi
    for(const MyNode* Mi_leaf: TreeTemplateTools::getLeaves(*Mi_roots.back()))
      corrMi[stid(Mi_leaf)] = p; //this leaf is in subtree Mi
    
    ++p; //number of Mi
	}
  ui.push_back(currentNode); // finally, push u_p
	const StId up_id = stid(currentNode);  // the stid of up


  // for each i, store the StIds of the leaves in Mi in order that they appear in T2
	vector<StId>* corrMiT2 = new vector<StId>[p];
	for(const StId leaf_id: orderedIds)
		corrMiT2[corrMi[leaf_id]].push_back(leaf_id); // to comment

  vector<MyTree*> Si(Mi_roots.size());
  vector<unsigned*> MiSi_mast(Mi_roots.size()); // store the MAST sizes between Mi and each vertex of Si
	for(unsigned i=0; i < Mi_roots.size(); ++i){
		Si[i] = T2->induced_subtree(corrMiT2[i]);// subtrees Mi in T2
    MiSi_mast[i] = mast(Mi_roots[i], Si[i]);
	}
	delete [] corrMiT2;



	// ============= STEP 3 ==================
  const unsigned CPs_in_T2 = T2->num_centroid_paths();

  const vector<MyNode*> nodes_T2 = T2->getNodes();
  // get the maximum StId occuring in T2 in order to know array bounds in the future
  unsigned max_StId = 0;
  for(const MyNode* v: nodes_T2) max_StId = std::max(stid(v), max_StId);
  //creating all the muligraphs we will need
  //NOTE: this will automatically setup the sets R(x), with the nodes ordered in preorder (each CP from root to leaf)
  vector<MultiGraph*> Gx(CPs_in_T2);
  for(unsigned i = 0; i < Gx.size(); ++i)
    Gx[i] = new MultiGraph(T2->get_centroid_path(i), max_StId);
    
  // step 3.1: add vertices and edges involving u_p to each G_x
  if(! T2->isPreprocessed()) T2->setClades();
  const MyNode& up = *ui.back();
  const StId up_StId = stid(&up);
  const MyNode* twin_up_in_T2 = T2->getNodeWithStId(up_StId);
  const Clade& clade_twin_up_in_T2 = twin_up_in_T2->getInfos().getClade();
  const Label v(up_StId, false);
  for(unsigned j = 0; j < CPs_in_T2; j++){
    const MyNode* rootOfCP = T2->root_of_centroid_path(j); //this is true because of the way RxVertices are added
    const Clade& clade_rootOfCP = rootOfCP->getInfos().getClade();
    if(clade_rootOfCP.contains(clade_twin_up_in_T2)){
      Gx[j]->addLxVertex(up_StId);
      const MyNode*  vq = T2->leaf_of_centroid_path(j);
      const MyNode* lca = T2->getLCA(vq, twin_up_in_T2);      
      assert(lca->getInfos().getCentroidPathNumber() == j); //the LCA should be in the same centroid path
      Gx[j].add_edge(up_StId, stid(lca), 1);
    }// if
  }// for
  

  // NOTE: this allocates quadratic space, but leaves most of it uninitialized
  StId map[Mi_roots.size()][T2->getCorrespondanceLenghtId()]; //this map will be useful later on

  // step 3.2: add vertices and edges involving u_i to each G_x
  for(ssize_t i = p - 1; i >= 0; --i){
    const vector<MyNode*> nodes_Si = Si[i]->getNodes();
    const StId ui_id = stid(ui[i]); 
     
    for(unsigned j = 0;j < nodes_Si.size(); j++){
      const MyNode* z = nodes_Si[j];
      if(z->hasFather()){
        const StId id_z = stid(z);
        const MyNode* parent_z_in_Si = z->getFather();
        const unsigned cp_parent_z_in_Si = parent_z_in_Si->getInfos().getCentroidPathNumber();
        MyNode* vj = T2->getNodeWithStId(id_z);
        while((vj->getInfos().getCentroidPathNumber() != cp_parent_z_in_Si) && T2->root_of_centroid_path_of(vj)->hasFather()){
          const MyNode* const root_Nj = vj;
          vj = T2->root_of_centroid_path_of(vj)->getFather();
          const StId id_vj = stid(vj);
          const StId y_id = (Si[i]->getNodeWithStId(id_vj) != NULL) ? id_vj : id_z;

          map[i][id_vj] = y_id;
          // if map(i,j) == v_j, the white weight is the mast between Mi & Si rooted at the child of v_j,
          // if map(i,j) != v_j, the white weight is the mast between Mi & Si rooted at z
          const int white_weight = mast_MiSi[i][ (y_id == id_vj) ? stid(root_Nj) : id_z ];
          // the green weight is the mast between Mi and map(i,j)
          const int green_weight = mast_MiSi[i][y_id];
          // computing the red weight needs the agreement matching of G_y
//          const int red_weight = -1;
          const MultiGraph& Gy = *Gx[root_Nj->getInfos()->getCentroidPathNumber()];
          const int red_weight = Gy.agreement_matching_below(y_id, LEFT);

          // add u_i to G_j and add the corresponding edge
          const unsigned cp_of_vj = vj->getInfos().getCentroidPathNumber();
          MultiGraph& G = *Gx[cp_of_vj]; //get the Gx to which add ui
          G.addLxVertex(ui_id);
          G.add_edge(ui_id, id_vj, white_weight, green_weight, red_weight);
        }// while
      }// if
    }// for each vertex z of S_i
  }// for all Si
 
  // compute array to return
  unsigned* result = new[nodes_T2.size()];
  for(const MyNode*& v: nodes_T2){
    const unsigned x = v->getInfos().getCentroidPath();
    const StId v_id = stid(v);
    result[v_id] = Gx[x]->agreement_matching_below(v_id, RIGHT);
  }
	return result;
}

unsigned* mast(MyTree* T1, MyNode* root_of_T2)
{
  MyTree* T2 = new MyTree(root_of_T2);
  // copy the leaf-StIds from T1 using the names
  T2->sync_leaf_stids(*T1);
  unsigned* result = mast(T1, T2);
  delete T2;
  return result;
}

unsigned* mast(MyNode* root_of_T1, MyTree* T2)
{
  MyTree* T1 = new MyTree(root_of_T1);
  // copy the leaf-StIds from T1 using the names
  T1->sync_leaf_stids(*T2);
  unsigned* result = mast(T1, T2);
  delete T1;
  return result;
}

#endif /*MAST_H_*/
