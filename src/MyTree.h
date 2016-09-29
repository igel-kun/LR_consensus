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

using namespace bpp;

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
  // NOTE: this assumes that idsLeaves is ordered with respect to a pre-order
#warning TODO: reorder the given leaves into a preorder at the beginning
  MyTree* induced_subtree(const vector<unsigned>& idsLeaves) const {
    assert(!idsLeaves.empty());
    vector<unsigned> ids; //holds the vertex IDs of all vertices in the induced subtree

    // get all LCA in tree above the leaves in idsLeaves
    // NOTE: by the assumed ordering of idsLeaves, no LCA is missed
    for(unsigned y = 0; y < idsLeaves.size() - 1; y++){
      const MyNode* yNode(getNodeWithStId(idsLeaves[y]));
      const MyNode* nextNode(getNodeWithStId(idsLeaves[y + 1]));
      ids.push_back(idsLeaves[y]);
      unsigned idNode = getLCA(yNode, nextNode)->getInfos().getStId(); //unsignedernal nodes
      ids.push_back(idNode);
    }
    ids.push_back(*idsLeaves.rbegin());

    MyTree* restrictedTree = new MyTree();
    // reserve enough space in the corresponding-ID vector for all ids
    restrictedTree->setCorrespondanceLenghtId(*std::max_element(ids.begin(), ids.end()));

    // use 2 instances of nearest_lower<> to find the vertex directly above y for each y with ID in ids
    nearest_lower<vector<unsigned>, DepthCompare<false>, DIR_LEFT > nearest_left_index_above(ids, DepthCompare<false>(*this));
    nearest_lower<vector<unsigned>, DepthCompare<false>, DIR_RIGHT> nearest_right_index_above(ids, DepthCompare<false>(*this));

    //defining the edge set
    for(unsigned y = 0; y < ids.size(); y++){
      // abbreviations:
      const unsigned y_id = ids[y];
      MyNode& y_node = *getNodeWithStId(y_id);

      // step1: construct new node to insert into the restrcited tree
      MyNode* newNode = new MyNode();
      newNode->getInfos().setStId(y_id);
      restrictedTree->setNodeWithStId(newNode, y_id);
      if(y_node.hasName()) newNode->setName(y_node.getName());

      // step2: find node above y_node using 2 instances of nearest_lower
      const ssize_t left_id_index = nearest_left_index_above.query(y);
      MyNode* v_left = (left_id_index == -1) ? NULL : restrictedTree->getNodeWithStId(ids[left_id_index]);
      
      const ssize_t right_id_index = nearest_right_index_above.query(y);
      MyNode* v_right = (right_id_index == -1) ? NULL : restrictedTree->getNodeWithStId(ids[right_id_index]);

      if(!v_left && v_right){
        v_right->addSon(newNode);
      } else if(v_left && !v_right){
        v_left->addSon(newNode);
      } else if(!v_left && !v_right){
        restrictedTree->setRootNode(newNode);
        //restrictedTree->resetNodesId();
      } else{
        if(v_right->getInfos().getDepth() < v_left->getInfos().getDepth())
          v_left->addSon(newNode);
        else
          v_right->addSon(newNode);
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

};


template<bool invert>
bool DepthCompare<invert>::operator()(unsigned u_id, unsigned v_id) const
{
  const MyNode* const u = tree.getNodeWithStId(u_id);
  const MyNode* const v = tree.getNodeWithStId(v_id);
  return (u->getInfos().getDepth() < v->getInfos().getDepth()) == invert; 
}


#endif /*MYTREE_H_*/
