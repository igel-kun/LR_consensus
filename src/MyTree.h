// Created by: Celine Scornavacca

#ifndef MYTREE_H_
#define MYTREE_H_


#include "rmq/lca.hpp"
#include "nearest_lower.hpp"

#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeExceptions.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

#include "NodeInfos.h"
#include "MatrixTriplets.h"
#include <iterator>
using namespace bpp;

using namespace std;

typedef NodeTemplate<NodeInfos> MyNode;

// forward declaration of MyTree
class MyTree;


// a class used to iterate over the children of a node, until Bio++ manages to return a const vector<NodeTemplate<NodeInfos>*>& for us...
template<class NodeType>
struct _ChildIterator{
  NodeType& node;
  unsigned index;

  _ChildIterator(NodeType& _node, const unsigned _index = 0):
    node(_node), index(_index)
  {}

  //! increment operator
  _ChildIterator& operator++() { ++index; return *this; }
  //! post-increment
  _ChildIterator& operator++(int) { unsigned i = index; ++(*this); return _ChildIterator(i); }
  //! dereferece
  NodeType& operator*() { return *node[index]; }
  //! comparison
  bool operator!=(const _ChildIterator& it) { return index != it.index; }
};
using ChildIterator = _ChildIterator<MyNode>;
using ChildConstIterator = _ChildIterator<const MyNode>;
// once Bio++ comes to its sences, replace this by a const vector<NodeTemplate<NodeInfos>*>
template<class NodeType>
struct _Children{
  NodeType& node;
  const size_t _size;

  _Children(NodeType& _node, const size_t _max): node(_node), _size(_max)   {}
  ChildConstIterator cbegin() const { return ChildConstIterator(node); }
  ChildIterator begin() const { return ChildIterator(node); }
  ChildConstIterator cend() const { return ChildConstIterator(node, _size); }
  ChildIterator end() const { return ChildIterator(node, _size); }
  size_t size() const {return _size; }
};
using Children = _Children<MyNode>;
using ConstChildren = _Children<const MyNode>;


// a comparator comparing nodes according to their depth DepthCompare<false>(u,v) <=> u has a smaller depth than v
template<bool invert = false>
struct DepthCompare{
  const MyTree& tree; //! holding the tree that helps map StIds to MyNodes to Depths
  DepthCompare(const MyTree& _tree): tree(_tree) {}

  // compare depths of nodes given by indices
  bool operator()(unsigned u_id, unsigned v_id) const;
};





class MyTree: public TreeTemplate<MyNode> {
public:

	vector<MyNode*> correspondanceId;
	vector<MyNode*> leavesPreordered;

	//vector< unsigned > centroidPaths ;

 	MyTree(): TreeTemplate<MyNode>(){}
 	MyTree(MyNode& root): TreeTemplate<MyNode>(&root) {}
 	MyTree(MyNode& root, unsigned dim): TreeTemplate<MyNode>(& root) {}

 	virtual ~MyTree(){
    if(lca_oracle) delete lca_oracle;
 	}

  bool is_root(const MyNode& node) const
  {
    return TreeTemplateTools::isRoot(node);
  }

  // get iterable children (this job should have been done by NodeTemplate<>...)
  ConstChildren get_children(const MyNode& node) const
  {
    return ConstChildren(node, node.getNumberOfSons());
  }
  Children get_children(MyNode& node) const
  {
    return Children(node, node.getNumberOfSons());
  }

  
  void setCorrespondanceLenghtId(unsigned dim){
 		correspondanceId.resize(dim);
 	}

  void resetCorrespondanceLenghtId(unsigned dim) {
    correspondanceId.clear();
 	}

	void setNodeWithStId(MyNode* node, unsigned stId) {
	 	correspondanceId[stId] = node;
	}

  void setup_StIds_subtree(MyNode* subtree_root, unsigned& offset){
    subtree_root->getInfos().setStId(offset);
    correspondanceId[offset] = subtree_root;
    std::cout << "setting StId "<<offset<<" to node @"<<subtree_root<<std::endl;
    for(auto& child: get_children(*subtree_root))
      setup_StIds_subtree(&child, ++offset);
  }
  // give stIds to the nodes of the tree starting with stId offset
  // if avoid_leaves is set, then do not give stIds to leaves
  void setup_StIds() {
    correspondanceId.resize(getNumberOfNodes());
    unsigned zero = 0;
    setup_StIds_subtree(getRootNode(), zero);
	}

	MyNode* getNodeWithStId(unsigned stId) const{
	 	return correspondanceId[stId];
	}

	void setDepth(MyNode& node) const {
	 	if(is_root(node))
	 		node.getInfos().setDepth(0); //root has depth 0
	 	else
	 		node.getInfos().setDepth(node.getFather()->getInfos().getDepth() + 1); //otherwise depth of father +1

    for(auto& child: get_children(node)) setDepth(child); // recursive calls for sons
	}

	void setNumberOfDescendents(MyNode& node) const{
    unsigned count = 0;
    for(auto& child: get_children(node)){
      setNumberOfDescendents(child); // recursive calls for sons
      count += child.getInfos().getNumberOfDescendents() + 1; //adding up all the descendants of sons + sons
    }
    node.getInfos().setNumberOfDescendents(count);
  }

	void setDepthAndNumberOfDescendents() {
    setDepth(*getRootNode());
    setNumberOfDescendents(*getRootNode());
  }

  //construct the centroid decomposition of a tree
	void getCentroidDecompostion(MyNode& node, unsigned partition_number) const {
    node.getInfos().setCentroidPathNumber(partition_number);
    unsigned max_desc = 0;
    unsigned other_partition_numbers = partition_number;
    MyNode* chosen_child = NULL;
    // get the child with the highest number of descendants
    for(auto& child: get_children(node)){
      unsigned num_desc = child.getInfos().getNumberOfDescendents();
      if(num_desc > max_desc){
        if(chosen_child)
          getCentroidDecompostion(*chosen_child, ++other_partition_numbers); // the other start new centroid paths
        chosen_child = &child;
        max_desc = num_desc;
      } else getCentroidDecompostion(child, ++other_partition_numbers); // the other start new centroid paths
    }
    //the child with the higher number of sons stay in the same centroid path than the father (ties broken arbitrarly)
    getCentroidDecompostion(*chosen_child, partition_number);
  }

  void getCentroidDecompostion() {
    getCentroidDecompostion(*getRootNode(), 0);
  }

	void setPreOrder(MyNode& node, unsigned& preorder_number) const {
    for(auto& child: get_children(node))
      setPreOrder(child, preorder_number);
		
		if(node.isLeaf()) node.getInfos().setPreOrder(++preorder_number);
  }

  void setPreOrder() {
    unsigned zero =0;
	 	setPreOrder(*getRootNode(), zero);
  }

	//not used for now
	 const vector<string> getNames(MyNode& node) const {
	   	vector<string> stId;
	  	if(node.isLeaf()) stId.push_back(node.getName());
	  	for(auto& child: get_children(node)){
	     	const vector<string> subStId = getNames(child);
        stId.insert(stId.end(), subStId.begin(), subStId.end());
      }
	   	return stId;
	 }

	vector<unsigned> getLeavesStId(MyNode& node) const {
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
  // preprocess can either be called before getLCA or getLCA will call it if it had not been called
  void lca_preprocess()
  {
    lca_oracle = new LCA_Oracle(*getRootNode(), correspondanceId.size(), access_Id(), access_children());
  }

  const MyNode* getLCA(const MyNode& u, const MyNode& v) const
  {
    if(!lca_oracle) throw TreeException("LCA query without preprocessing", this);
    return lca_oracle->query(u, v);
  }
  const MyNode* getLCA(const MyNode * u, const MyNode * v) const
  {
    if(!lca_oracle) throw TreeException("LCA query without preprocessing", this);
    return lca_oracle->query(*u, *v);
  }


  //======================= computing subtree induced by leaves ===============
  // restricts the given tree to the given leaves (given by StId)
  // NOTE: assumptions:
  //  1. idsLeaves is ordered with respect to a pre-order
  //  2. all depths and StIds have been set up
  MyTree* induced_subtree(const vector<unsigned>& idsLeaves) const {
    assert(!idsLeaves.empty());
    

    // Step 1: get the IDs of all nodes in the restricted tree
    // get all LCA in tree above the leaves in idsLeaves and construct the new nodes
    // NOTE: by the assumed ordering of idsLeaves, no LCA is missed
    vector<unsigned> ids; //holds the vertex IDs of all vertices in the induced subtree
    ids.reserve(2 * idsLeaves.size());
    for(unsigned y = 0; y < idsLeaves.size() - 1; y++){
      const unsigned y_id = idsLeaves[y];
      const MyNode* const y_node(getNodeWithStId(y_id));
      ids.push_back(y_id);

      const MyNode* const nextNode(getNodeWithStId(idsLeaves[y + 1]));
      unsigned LCA_id = getLCA(y_node, nextNode)->getInfos().getStId(); //internal nodes
      ids.push_back(LCA_id);
    }
    ids.push_back(*idsLeaves.rbegin());


    // Step 2: construct the restricted tree and reserve space for the correspondance
    MyTree* restrictedTree = new MyTree();
    // reserve enough space in the corresponding-ID vector for all ids
    restrictedTree->setCorrespondanceLenghtId(*std::max_element(ids.begin(), ids.end()));


    // Step 3: create all vertices according to the ids and set their correspondance in restrcitedTree
    vector<MyNode*> nodes;
    nodes.reserve(2 * idsLeaves.size());
    for(unsigned i = 0; i < ids.size(); ++i){
      const unsigned y_id = ids[i];
      const MyNode* const y_node(getNodeWithStId(y_id));
      MyNode* newNode = new MyNode();
      newNode->getInfos().setStId(y_id);
      restrictedTree->setNodeWithStId(newNode, y_id);
      if(y_node->hasName())
        newNode->setName(y_node->getName());
      nodes.push_back(newNode);
    }

    
    // Step 4: defining the edge set
    // use 2 instances of nearest_lower<> to find the vertex directly above y for each y with ID in ids
    nearest_lower<vector<unsigned>, DIR_LEFT,  DepthCompare<false> > nearest_left_index_above(ids, DepthCompare<false>(*this));
    nearest_lower<vector<unsigned>, DIR_RIGHT, DepthCompare<false> > nearest_right_index_above(ids, DepthCompare<false>(*this));
    for(unsigned y = 0; y < nodes.size(); y++){
      // find node above y_node using 2 instances of nearest_lower
      const ssize_t left_id_index = nearest_left_index_above.query(y);
      MyNode* v_left = (left_id_index == -1) ? NULL : nodes[left_id_index];
      
      const ssize_t right_id_index = nearest_right_index_above.query(y);
      MyNode* v_right = (right_id_index == -1) ? NULL : nodes[right_id_index];

      if(!v_left && v_right){
        v_right->addSon(nodes[y]);
      } else if(v_left && !v_right){
        v_left->addSon(nodes[y]);
      } else if(!v_left && !v_right){
        restrictedTree->setRootNode(nodes[y]);
        //restrictedTree->resetNodesId();
      } else{
        const unsigned left_depth = getNodeWithStId(ids[left_id_index])->getInfos().getDepth();
        const unsigned right_depth = getNodeWithStId(ids[right_id_index])->getInfos().getDepth();
        if(right_depth < left_depth){
          v_left->addSon(nodes[y]);
        } else{
          v_right->addSon(nodes[y]);
        }
      }
    }

    restrictedTree->resetNodesId();
    return restrictedTree;

  }

  /*
  // ====================== output ======================
  void pretty_print(const MyNode& node, std::ostream& os = std::cout) const
  {
    const unsigned num_desc = node.getInfos().getNumberOfDescendents();
    char to_print = ' ';
    const auto children = get_children(node);
    for(ChildConstIterator child_it = children.cbegin(); child_it != children.cend(); ++child_it){
      const unsigned num_child_desc = child_it->getInfos().getNumberOfDescendants();
      for(unsigned i = 0; i < num_child_desc; ++i) os << to_print;

//      to_print = 
    }
  }
  */

  void pretty_print(std::ostream& os = std::cout) const
  {
  }
  
  //clades stuff
  

  void setClades(MyNode & currentNode){
  	if(currentNode.isLeaf()){
  		pair <int,int> ids(leavesPreordered.size(),leavesPreordered.size()); // preprocess the clade of a leaf for later
  	  	currentNode.getInfos().setClade(ids);
  		leavesPreordered.push_back(& currentNode);
  	}	
  	else{
  		int startClade = leavesPreordered.size();
		for(auto& child: get_children(currentNode))
      		setClades(child); // recursive calls for sons
       int endClade = leavesPreordered.size()-1;
	   pair <int,int> ids(startClade,endClade);
       currentNode.getInfos().setClade(ids);
    }  		
  }
  
  void setClades(){
  	setClades(* getRootNode());
  }
  
  
  
  //triplets stuff
  

	
  void setTriplets(MatrixTriplets & Triplets, MyNode & currentNode){   
		
		for(auto& child: get_children(currentNode)){
      		setTriplets(Triplets, child); // recursive calls for sons
		}
		
		pair<int,int> cladeCN= currentNode.getInfos().getClade();

		for(auto childitA = get_children(currentNode).begin();childitA !=get_children(currentNode).end();childitA++){
			pair<int,int> cladeA= (* childitA).getInfos().getClade();
			for(auto childitB = next(childitA);childitB !=get_children(currentNode).end();childitB++){
				pair<int,int> cladeB= (* childitB).getInfos().getClade();
				
				for(unsigned int i = cladeA.first;i <= cladeA.second; i++){
					for(unsigned int j = cladeB.first;j < cladeB.second; j++){
						for(unsigned int z = 0;z < leavesPreordered.size(); z++){
							if(z==cladeCN.first)
								z=cladeCN.second;
							else
								Triplets.add(* leavesPreordered[i],* leavesPreordered[j], * leavesPreordered[z]);
						}	
					}
				}
			}
		}
	}
	
	MatrixTriplets * setTriplets(int dim ){   
		MatrixTriplets * Triplets =new MatrixTriplets();
		Triplets->setDim(dim);
		setTriplets(* Triplets, * getRootNode());
		return Triplets;
 }	

};


template<bool invert>
bool DepthCompare<invert>::operator()(unsigned u_id, unsigned v_id) const
{
  const MyNode* const u = tree.getNodeWithStId(u_id);
  const MyNode* const v = tree.getNodeWithStId(v_id);
  return (u->getInfos().getDepth() < v->getInfos().getDepth()) != invert;
}


#endif /*MYTREE_H_*/
