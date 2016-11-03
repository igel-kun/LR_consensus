// Created by: Celine Scornavacca

#ifndef MYTREE_H_
#define MYTREE_H_

#include <iterator>
#include <map>
#include <list>

#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeExceptions.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

#include "utils.h"
#include "rmq/lca.hpp"
#include "nearest_lower.hpp"
#include "tree_iters.h"

#include "NodeInfos.h"
#include "MatrixTriplets.h"


using namespace bpp;

using namespace std;

typedef map<string, const MyNode*> LeafAssociation;

// forward declaration of MyTree
class MyTree;


template<class NodeType, class Iter>
struct _Traversal
{
  NodeType* root;

  _Traversal(NodeType* _r): root(_r) {}
  Iter begin() const { return Iter(root); }
  Iter end() const { return Iter(); }
};

using PreOrderTraversal = _Traversal<MyNode, PreOrderIterator>;
using PreOrderConstTraversal = const _Traversal<const MyNode, PreOrderConstIterator>;
using PostOrderTraversal = _Traversal<MyNode, PostOrderIterator>;
using PostOrderConstTraversal = const _Traversal<const MyNode, PostOrderConstIterator>;




// a comparator comparing nodes according to their depth DepthCompare<false>(u,v) <=> u has a smaller depth than v
template<bool invert = false>
struct DepthCompare{
  const MyTree& tree; //! holding the tree that helps map StIds to MyNodes to Depths
  DepthCompare(const MyTree& _tree): tree(_tree) {}

  // compare depths of nodes given by indices
  bool operator()(unsigned u_id, unsigned v_id) const;
};



class MyTree: public TreeTemplate<MyNode> 
{
protected:
	vector<MyNode*> StId_to_node; //! conversion array for StId -> node
	vector<MyNode*> leaves_po;    //! conversion array for leaf postorder numbers -> leaves
  vector<MyNode*> nodes_po;     //! conversion array for node postorder numbers -> nodes
  MatrixTriplets triplets;      //! a tree knows its triplets

  //vector<CentroidPath> centroid_paths;

	//vector< unsigned > centroidPaths ;

  //! copy the StIds of the nodes in T, assuming that we agree with T on our common leaves
  /** Here, "supertree" is in the phylogenetic sense (that is, modulo degree-2 nodes)
   * The return value is the new StId of "root", unless exactly one subtree of "root" in our tree is represented in T,
   *    in which case the StId of the root of this subtree is returned. This allows skipping over indeg-1 & outdeg-1 nodes
   * name_to_leaf should associate leaf names to the leaves in T **/
  StId sync_stids_from_leaf_names(const MyTree& T, MyNode& root, const LeafAssociation& name_to_leaf)
  {
    StId my_stid;
    if(root.isLeaf()){
      const auto node_iter = name_to_leaf.find(root.getName());
      if(node_iter == name_to_leaf.end()){
        // if this leaf does not exist in T, create a new StId at the end of StId_to_node
        return setNewStId(&root);
      } else my_stid = stid(node_iter->second);
    } else{
      // get the children of root that are represented in T
      vector<const MyNode*> childs_in_T;
      childs_in_T.reserve(root.getNumberOfSons());
      for(auto& child: get_children(root)){
        const StId child_id = sync_stids_from_leaf_names(T, child, name_to_leaf);
        // if the child is represented in T
        if(child_id < T.StId_to_node.size()) {
          childs_in_T.push_back(T[child_id]);
          assert(stid(childs_in_T.back()) == child_id);
        }
      }
      switch(childs_in_T.size()){
        case 0:
          // if no child is represnted in T, then set a new StId and return it
          return setNewStId(&root);
        case 1:
          // if only one child is represented in T, then set a new StId, but return this childs StId; this skips over degree-2-paths
          setNewStId(&root);
          return stid(childs_in_T.front());
        default:
          // if at least two children are represented in T, copy the StId of their father in T
          my_stid = stid(childs_in_T.front()->getFather());
      }
    }
    stid(&root) = my_stid;
    StId_to_node[my_stid] = &root;
    return my_stid;
  }

public:
 	MyTree(): TreeTemplate<MyNode>(){}

  // copy the subtree rooted at _root from somewhere, StIds are NOT copied
  MyTree(MyNode* _root): TreeTemplate<MyNode>(TreeTemplateTools::cloneSubtree<MyNode>(*_root)) {}

  //! copy the given tree and synchronize StIds, setting-up StId_to_node
  MyTree(const MyTree& T): TreeTemplate<MyNode>(T)
  {
    // we like to use Bpp's TreeTemplate cppy constructor, but it doesn't copy infos because... reasons
    // so we have to associate the leaves using their names (which ARE copied) and then copy StIds bottom-up
    sync_stids_from(T);
  }

 	virtual ~MyTree(){ if(lca_oracle) delete lca_oracle; }

  bool is_root(const MyNode& node) const
  {
    return TreeTemplateTools::isRoot(node);
  }

  //! return the number of nodes, assuming the node infos are set up
  size_t num_nodes() const { assert(node_infos_set_up()); return getRootNode()->getInfos().subtree_size; }
  //! return the number of leaves, assuming the node infos are set up
  size_t num_leaves() const { return leaves_po.size(); }
  //! return the size of the StId storage
  size_t num_stids() const { return StId_to_node.size(); }


  //! compute an association of leaf-names to leaf-pointers
  void compute_leaf_association(LeafAssociation& name_to_leaf) const
  {
    for(const auto& leaf: getLeaves()) name_to_leaf.emplace(leaf->getName(), leaf);
  }

  //! copy StIds from T, assuming that we agree on T on our common leaves (that is, our subtrees induced by the common leaves are isomorphic)
  /** NOTE: if the topoly of T differs from our own, the behavior is undefined **/
  void sync_stids_from(const MyTree& T, LeafAssociation* name_to_leaf = NULL)
  {
    const bool delete_asso = (name_to_leaf == NULL);
    if(delete_asso) {
      // step 1: associate leaves
      name_to_leaf = new LeafAssociation();
      T.compute_leaf_association(*name_to_leaf);
    }

    // copy the StIds bottom-up from the corresponding leaves
    StId_to_node.resize(T.StId_to_node.size());
    sync_stids_from_leaf_names(T, *getRootNode(), *name_to_leaf);

    if(delete_asso) delete name_to_leaf;
  }

  //! give a specific StId to a specific node. This might mean that we have to swap with someone
  void assign_StId(MyNode* const leaf, const StId stid_want)
  {
    const StId stid_have = stid(leaf);
    if(stid_want != stid_have){
      //cout << "swapping stid("<<leaf->getName()<<")="<<stid_have<<" with "<<stid_want<<" (parent stid "<<stid(name_to_leaf->at(leaf->getName())->getFather())<<"), storage size is "<<StId_to_node.size()<<endl;
      // make sure there's enough room in the StId_to_node vector
      if(StId_to_node.size() > stid_want) {
        MyNode* swap_target = StId_to_node[stid_want];
        if(swap_target && (stid(swap_target) == stid_want)){
          stid(swap_target) = stid_have;
          StId_to_node[stid_have] = swap_target;
        }
      } else {
        StId_to_node.resize(stid_want + 1);
        StId_to_node[stid_have] = NULL;
      }
      StId_to_node[stid_want] = leaf;
      stid(leaf) = stid_want;
    }
  }

  //! setup StIds such that all leaves with the same label have the same StIds; all other StIds are (more or less) random
  void sync_leaf_stids(const MyTree& T, LeafAssociation* name_to_leaf = NULL)
  {
    const bool delete_association = (name_to_leaf == NULL);
    if(delete_association){
      name_to_leaf = new LeafAssociation();
      T.compute_leaf_association(*name_to_leaf);
    }

    // for each leaf, switch the StId's of that leaf and the node with the StId that we want
    for(MyNode* leaf: getLeaves()) assign_StId(leaf, stid(name_to_leaf->at(leaf->getName())));

    if(delete_association) delete name_to_leaf;
  }

  //! consolidate StIds such that the tree uses consequtive StIds
  //NOTE: if dont_touch_leaves is set, the leaf StIds will remain in tact. In this case, the StIds might not be consolidated
  void consolidate_StIds(const bool dont_touch_leaves = true)
  {
    bool seen_leaf = false;
    StId first = 0;
    while(StId_to_node.back() == NULL) StId_to_node.pop_back();
    StId last = StId_to_node.size() - 1;
    while(first < last){
      while((StId_to_node[first] != NULL) && (first < StId_to_node.size())) ++first;
      if(first < StId_to_node.size()) {
        // at this point, first points to the first NULL element, and last to the last non-NULL element
        MyNode* const u = StId_to_node[last];
        if(dont_touch_leaves && u->isLeaf()){
          seen_leaf = true;
          --last;
        } else assign_StId(u, first);
      } else return; // if all other StIds are given, we can savely return
      if(!seen_leaf)
        while(StId_to_node.back() == NULL) { StId_to_node.pop_back(); --last; }
      else
        while((StId_to_node[last] == NULL) && (last > 0)) --last;
    }// while first < last
  }

  
  ConstChildren get_children_by_StId(const StId& id) const
  {
    return get_children((const MyNode&)*StId_to_node[id]);
  }
  // get iterable traversals
  PreOrderConstTraversal preorder_traversal() const
  {
    return PreOrderConstTraversal(getRootNode());
  }
  PreOrderTraversal preordero_traversal()
  {
    return PreOrderTraversal(getRootNode());
  }
  PostOrderConstTraversal postorder_traversal() const
  {
    return PostOrderConstTraversal(getRootNode());
  }
  PostOrderTraversal postorder_traversal()
  {
    return PostOrderTraversal(getRootNode());
  }

  // return a copy of us that misses the vertex with the given StId
  MyTree* operator-(const StId x) const 
  {
    MyTree* result = new MyTree(*this);
    *result -= x;
    return result;
  }

  // return a copy of us that misses the vertex with the given StId
  MyTree* operator-(const MyNode* x) const 
  {
    MyTree* result = new MyTree(*this);
    *result -= x;
    return result;
  }

  // delete node from the tree
  void operator-=(const StId x)
  {
    if(x >= StId_to_node.size()) throw NodeNotFoundException("operator-=(): could't find node by StId: ", x);
    MyNode* x_node = StId_to_node[x];
    StId_to_node[x] = NULL;
    if(x_node->hasFather()) {
      MyNode* parent = x_node->getFather();

      // for some reason, Bio++ decided that removing x_node from the parent child-list causes x_node to become a leaf...
      const bool x_is_leaf = x_node->isLeaf();

      parent->removeSon(x_node);
      if(!x_is_leaf){
        for(auto& child: get_children(*x_node)){ 
          child.removeFather();
          parent->addSon(&child);
        }
      } else operator-=(parent);
    } else {
      // x_node is the root
      if(x_node->getNumberOfSons() != 1) throw Exception("can only remove the root if it has exactly one child");
      MyNode* child = x_node->getSon(0);
      child->removeFather();
      x_node->removeSon(child);
      setRootNode(child);
    }
    delete x_node;
  }

  void operator-=(const MyNode* const node)
  {
    operator-=(stid(node));
  }

  //! take the leaf x and regraft it above v, keeping StIds in tact
  void regraft_leaf_above(MyNode* const x, MyNode* const v)
  {
    assert(x->hasFather());
    MyNode* const x_parent = x->getFather();
    assert(v != x_parent); // if v is x_parent then grafing x above v makes v deg-2 and we end up with the same tree

    //cout << "regrafting "<<x->getName()<<" onto "<<stid(v)<<" in "<<endl;
    //pretty_print(cout, true);
    if(x_parent->hasFather()){
      
      MyNode* const x_grand_parent = x_parent->getFather();
      // step 1: remove x_parent from its parent
      x_grand_parent->removeSon((Node*)x_parent);
      // step 1: move all children of x_parent onto their grand parent, except x
      for(ssize_t i = x_parent->getNumberOfSons() - 1; i >= 0; --i){
        MyNode* const child = x_parent->getSon(i);
        x_parent->removeSon(i);
        if(child != x)
          x_grand_parent->addSon(child);
      }
      // now, graft x_parent and x onto v
      if(!v->hasFather()){
        // if v is the root
        x_parent->addSon(v);
        x_parent->addSon(x);
        setRootNode(x_parent);
      } else {
        MyNode* const v_parent = v->getFather();
        v_parent->removeSon((Node*)v);
        v_parent->addSon(x_parent);
        x_parent->addSon(v);
        x_parent->addSon(x);
      }
    } else {
      assert(v->hasFather());
      MyNode* const v_parent = v->getFather();
      //NOTE: if x and v have a common parent, then the input tree is equal to the output tree, so don't do anything
      if(v_parent == x_parent) return;
      // x_parent is the root
      assert(x_parent->getNumberOfSons() == 2);
      const unsigned new_root_son = (x_parent->getSon(0) == x) ? 1 : 0;
      MyNode* const new_root = x_parent->getSon(new_root_son);

      //cout << stid(new_root)<<" is going to be the new root"<<endl;
      x_parent->removeSon(new_root_son);
      //NOTE: v cannot be the new root, as in this case v->getFather() == x_parent
      v_parent->removeSon((Node*)v);
      v_parent->addSon(x_parent);
      x_parent->addSon(v);
      setRootNode(new_root);
    }
    //cout << "done regrafting: "<<endl;
    //pretty_print();
  }

  void regraft_leaf_above(const StId x, const StId v)
  {
    regraft_leaf_above(StId_to_node[x], StId_to_node[v]);
  }

  //! setup traversal numbers, subtree sizes, and StIds
  // if requested, returns the vertices in post order; in any case, return our tree to enable chaining preprocessings
  //NOTE: unset change_StIds to keep StIds in tact
  MyTree* setup_node_infos(const bool change_StIds = true)
  {
    PostOrderTraversal po = postorder_traversal();
    leaves_po.clear();
    nodes_po.clear();
    if(change_StIds) StId_to_node.clear();

    for(auto u_iter = po.begin(); u_iter.is_valid(); ++u_iter){
      MyNode& u = *u_iter;
      NodeInfos& u_info = u.getInfos();
      u_info.subtree_size = 1;
      u_info.depth = u_iter.current_depth();
      u_info.po_num = nodes_po.size();
      nodes_po.push_back(&u);

      if(change_StIds) setNewStId(&u);
      if(u.isLeaf()){
        u_info.clade = Clade(leaves_po.size());
  		  leaves_po.push_back(&u);
      } else {
        u_info.clade = u.getSon(0)->getInfos().clade;
        for(const MyNode& v: get_children(u)){
          u_info.clade = get_spanning_clade(u_info.clade, v.getInfos().clade);
          u_info.subtree_size += v.getInfos().subtree_size;
        }
      }
    }
    //cout << "done setting up "<<nodes_po.size()<<" nodes"<<endl;
    return this;
  }

  bool node_infos_set_up() const
  {
    //cout << leaves_po.size() << " leaves "<<nodes_po.size()<<" nodes"<<endl;
    if(leaves_po.empty()) return false;
    //cout << "root-po: "<<getRootNode()->getInfos().po_num<<" root-subtree: "<<getRootNode()->getInfos().subtree_size<<endl;
    return (getRootNode()->getInfos().po_num + 1 == getRootNode()->getInfos().subtree_size);
  }


  //! give a previously unused StId (at the end of StId_to_node) to the node
  StId setNewStId(MyNode* node){
    const StId my_stid = StId_to_node.size();
    stid(node) = my_stid;
    StId_to_node.push_back(node);
    return my_stid;
  }

  //! translate StId to node using []
  MyNode* operator[](const StId id) const
  {
    return StId_to_node[id];
  }

  MyNode* leaf_by_po_num(const unsigned po_num) const
  {
    return leaves_po[po_num];
  }

  MyNode* node_by_po_num(const unsigned po_num) const
  {
    return nodes_po[po_num];
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
	     	stId.push_back(stid(&node));
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
    unsigned operator[](const MyNode& u) const { return stid(&u); }
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
  //NOTE: we'll return ourself in order to allow chaining preprocessings
  MyTree* lca_preprocess(const bool force = true)
  {
    if(!lca_oracle || force) {
      if(lca_oracle) delete lca_oracle;
      lca_oracle = new LCA_Oracle(*getRootNode(), StId_to_node.size(), access_Id(), access_children());
    }
    return this;
  }
  bool is_lca_preprocessed() const
  {
    return (lca_oracle != NULL);
  }

  const MyNode* getLCA(const MyNode& u, const MyNode& v) const
  {
    if(!lca_oracle) throw TreeException("LCA query without preprocessing", this);
    return lca_oracle->query(u, v);
  }
  const MyNode* getLCA(const MyNode* u, const MyNode * v) const
  {
    if(!lca_oracle) throw TreeException("LCA query without preprocessing", this);
    return lca_oracle->query(*u, *v);
  }

  //======================= computing subtree induced by leaves ===============
  // restricts the given tree to the given leaves (given by StId)
  // NOTE: assumptions:
  //  1. idsLeaves is ordered with respect to a pre- or post-order
  //  2. all depths and StIds have been set up
  template<class IdContainer = vector<unsigned> > // we also accept a list of stids
  MyTree* induced_subtree(const IdContainer& idsLeaves) const 
  {
    assert(!idsLeaves.empty());
    
    // Step 1: get the IDs of all nodes in the restricted tree
    // get all LCA in tree above the leaves in idsLeaves and construct the new nodes
    // NOTE: by the assumed ordering of idsLeaves, no LCA is missed
    vector<unsigned> ids; //holds the vertex IDs of all vertices in the induced subtree
    ids.reserve(2 * idsLeaves.size());
    for(auto y_iter = idsLeaves.begin();;){
      const unsigned y_id = *y_iter;
      const MyNode* const y_node(StId_to_node[y_id]);
      ids.push_back(y_id);

      if(++y_iter != idsLeaves.end()) {
        const MyNode* const nextNode(StId_to_node[*y_iter]);
        const StId lca_id = stid(getLCA(*y_node, *nextNode));
        ids.push_back(lca_id);
      } else break;
    }
    //cout << "ids in the induced subtree: "<<ids<<endl;

    // Step 2: construct the restricted tree and reserve space for the correspondance
    MyTree* restrictedTree = new MyTree();
    // reserve enough space in the corresponding-ID vector for all ids
    restrictedTree->StId_to_node.resize(*std::max_element(ids.begin(), ids.end()) + 1);


    // Step 3: create all vertices according to the ids and set their correspondance in restrcitedTree
    vector<MyNode*> nodes;
    nodes.reserve(2 * idsLeaves.size());
    for(unsigned i = 0; i < ids.size(); ++i){
      const unsigned y_id = ids[i];
      const MyNode* const y_node(StId_to_node[y_id]);
      MyNode* newNode = new MyNode(y_node->getId());
      stid(newNode) = y_id;
      restrictedTree->StId_to_node[y_id] = newNode;
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
        const unsigned left_depth = StId_to_node[ids[left_id_index]]->getInfos().depth;
        const unsigned right_depth = StId_to_node[ids[right_id_index]]->getInfos().depth;
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

  // ====================== output ======================
  void pretty_print(const string prefix, const MyNode* const root, const bool print_id = false, std::ostream& os = std::cout) const
  {
    if(!root->isLeaf()){
      const string rootid = to_string(stid(root)) + (string)(print_id ? " [" + to_string(root->getId()) + "]" : "");
      os << "-" << rootid;
      string new_prefix = prefix + string(rootid.length(), ' ');
      for(unsigned i = 0; i < root->getNumberOfSons() - 1; ++i){
        pretty_print(new_prefix + '|', root->getSon(i), print_id, os);
        os << new_prefix + '+';
      }
      pretty_print(new_prefix + ' ', root->getSon(root->getNumberOfSons()-1), print_id, os);
    } else os << "-" << stid(root)<<" (" << (string)(root->hasName() ? root->getName() : "") << ")" + (string)(print_id ? " [" + to_string(root->getId()) + "]" : "")<<std::endl;
    
  }

  void pretty_print(std::ostream& os = std::cout, const bool print_id = false) const
  {
    pretty_print("", getRootNode(), print_id, os);
  }
  
  // ====================== conflict triples ======================
	
  void setup_triplets(const MyNode& current_node) 
  {
		for(auto& child: get_children(current_node)) setup_triplets(child); // recursive calls
		
		const Clade& cladeCN = current_node.getInfos().clade;
    const ConstChildren childs = get_children(current_node);
    
		for(auto childitA = childs.begin(); childitA != childs.end(); ++childitA){
			const Clade& cladeA = childitA->getInfos().clade;
			for(auto childitB = childitA + 1; childitB != childs.end(); ++childitB){
				const Clade& cladeB = childitB->getInfos().clade;
				for(unsigned i = cladeA.first;i <= cladeA.second; i++){
					for(unsigned j = cladeB.first;j <= cladeB.second; j++){
						for(unsigned z = 0; z < cladeCN.first; ++z)
								add_triple(triplets, i, j, stid(leaves_po[z]));
						for(unsigned z = cladeCN.second + 1; z < num_leaves(); ++z)
								add_triple(triplets, i, j, stid(leaves_po[z]));
					}// for all leaves j in the clade of B
				}// for all leaves i in the clade of A
			}// for each child B that is right of A
		}// for each child A of the current node
	}

  //! gets all conflicting triples
  /** NOTE: this assumes that the clades have been set up*/
	void setup_triplets() {
    assert(!leaves_po.empty());

		triplets.resize(num_leaves(), num_leaves());
		setup_triplets(*getRootNode());
  }

  bool triplets_set_up() const
  {
    return triplets.size() == make_pair(num_leaves(), num_leaves());
  }

  const MatrixTriplets& get_triplets() const
  {
    return triplets;
  }

  //! after regrafting a leaf x, this will update all triplets involving x
  void update_triplets(const unsigned old_po_num, const unsigned new_po_num)
  {
#warning TODO: write me
  }

/*
  //============================ centroid paths =============================
  //construct the centroid decomposition of a tree
	void getCentroidDecompostion(MyNode& node, unsigned partition_number, const bool is_root_of_path = true) {
    // make centroid_paths fit
    if(centroid_paths.size() <= partition_number) centroid_paths.resize(partition_number + 1);
    // set the centroid path of node and, if it's a root, remember that
    if(is_root_of_path) centroid_paths[partition_number].set_root(&node);
    node.getInfos().cp_num = partition_number;
    // recurse unless node is a leaf
    if(!node.isLeaf()){
      unsigned max_desc = 0;
      MyNode* chosen_child = NULL;
      // get the child with the highest number of descendants
      for(auto& child: get_children(node)){
        const unsigned num_desc = child.getInfos().subtree_size;
        if(num_desc > max_desc){
          // if we find a new chosen child, recurse into the old chosen child first
          if(chosen_child){
            getCentroidDecompostion(*chosen_child, centroid_paths.size()); // the other start new centroid paths
          }
          chosen_child = &child;
          max_desc = num_desc;
        } else getCentroidDecompostion(child, centroid_paths.size()); // the other start new centroid paths
      }
      //the child with the higher number of sons stay in the same centroid path than the father (ties broken arbitrarly)
      getCentroidDecompostion(*chosen_child, partition_number, false);
    } else centroid_paths[partition_number].set_leaf(&node);
  }

  void getCentroidDecompostion() {
    centroid_paths.clear();
    // NOTE: there are as many centroid paths in a tree as there are leaves
    centroid_paths.reserve(StId_to_node.size() / 2); // at most 1/2 of the vertices of the tree are leaves
    getCentroidDecompostion(*getRootNode(), 0);
  }
  unsigned num_centroid_paths() const { return centroid_paths.size(); }
  const CentroidPath& get_centroid_path(const unsigned i) const { return centroid_paths[i]; }

#define root_of_centroid_path(x) get_centroid_path(x).get_root()
#define leaf_of_centroid_path(x) get_centroid_path(x).get_leaf()

  MyNode* root_of_centroid_path_of(const MyNode* node) const 
  { 
    return centroid_paths[node->getInfos().cp_num].get_root();
  }
*/

};


template<bool invert>
bool DepthCompare<invert>::operator()(unsigned u_id, unsigned v_id) const
{
  const MyNode* const u = tree[u_id];
  const MyNode* const v = tree[v_id];
  return (u->getInfos().depth < v->getInfos().depth) != invert;
}


/* this function reads a list of trees written in a file (path) in a newich format, separed by semicolons and returns a list of MyTree*/
vector < MyTree *>  readTrees(const string & path) throw (Exception) {
    // Checking the existence of specified file
    
    ifstream file(path.c_str(), ios::in);
    //if (! file) { throw IOException ("\nError reading file.\nInvalid options!\nUsage:\n ./physic -s sourceTreeFile -t threshold.\nwhere:\n - sourceTreeFile contains a set of rooted trees in newick format with bootstrap values and possibly edge lengths.\n - threshold indicates bootstrap values under which clades are not considered for building the supertree\n(typically a threshold of 70 can be used when source trees where obtained from 100 bootstrap replicates).\n"); }
    if (! file.good()) { throw IOException ("\nError reading file.\n"); }
    
    vector<MyTree*> trees;
    string temp, description;
    while (file.good()) {
        temp = FileTools::getNextLine(file);
        if(temp.size() != 0){
            string::size_type index = temp.find(";");
            if(index == string::npos) throw Exception("readTrees(). Bad format: no semi-colon found.");
            if(index < temp.size()) {
                description += temp.substr(0, index + 1);   
                TreeTemplate<Node> * tree = TreeTemplateTools::parenthesisToTree(description,true);
                trees.push_back(new MyTree((MyNode*)(tree->getRootNode())));
                delete tree;
                description = temp.substr(index + 1);   
            } else description += temp;
        }
    }
    file.close();
    return trees;   
};

void collapseEdge(MyTree & tree, MyNode * on) {
  MyNode * temp = (* on).getFather();
  MyNode * newNode;
  unsigned i_max = (* on).getNumberOfSons();  
  for (unsigned i=0; i< i_max; i++){
  
      if(((* on).getSon(0))->hasDistanceToFather() && ( on)->hasDistanceToFather())
          (* (* on).getSon(0)).setDistanceToFather(((* on).getSon(0))->getDistanceToFather() + ( on)->getDistanceToFather());
          
      newNode=(* on).getSon(0);   // we take always the first son...it's always a different one
      
      (* on).removeSon((* on).getSon(0));
      temp->addSon( newNode);
  } 
  (* temp).removeSon( on);
  if((* temp).getNumberOfSons()==1 && (((* temp).hasFather())))   {  //degree ==2... not so good! we have to collapse again
      collapseEdge(tree, temp);
  }
}



#endif /*MYTREE_H_*/
