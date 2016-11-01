

#pragma once

#include "NodeInfos.h" 


// a class used to iterate over the children of a node, until Bio++ manages to return a const vector<NodeTemplate<NodeInfos>*>& for us...
template<class NodeType>
struct _ChildIterator
{
  NodeType& node;
  unsigned index;

  _ChildIterator(NodeType& _node, const unsigned _index = 0):
    node(_node), index(_index)
  {}

  //! increment operator
  _ChildIterator& operator++() { ++index; return *this; }
  //! post-increment
  _ChildIterator operator++(int) { unsigned i = index; ++(*this); return _ChildIterator(node, i); }
  //! addition
  _ChildIterator operator+(const int& j) const { return _ChildIterator(node, index + j); }
  //! dereferece
  NodeType& operator*() const { return *node[index]; }
  NodeType* operator->() const { return node[index]; }
  //! comparison
  bool operator!=(const _ChildIterator& it) const { return index != it.index; }
};
using ChildIterator = _ChildIterator<MyNode>;
using ChildConstIterator = _ChildIterator<const MyNode>;
// once Bio++ comes to its sences, replace this by a const vector<NodeTemplate<NodeInfos>*>
template<class NodeType>
struct _Children{
  typedef _ChildIterator<NodeType> Iter;

  NodeType& node;
  const size_t _size;

  _Children(NodeType& _node, const size_t _max): node(_node), _size(_max)   {}
  Iter begin() const { return Iter(node); }
  Iter end() const { return Iter(node, _size); }
  size_t size() const {return _size; }
};
using Children = _Children<MyNode>;
using ConstChildren = _Children<const MyNode>;


// get iterable children (this job should have been done by NodeTemplate<>...)
ConstChildren get_children(const MyNode& node)
{
  return ConstChildren(node, node.getNumberOfSons());
}
Children get_children(MyNode& node)
{
  return Children(node, node.getNumberOfSons());
}





//! an iterator visiting the tree in preorder
template<class NodeType>
class _PreOrderIterator
{
private:
  NodeType* node;
  stack<unsigned> child_index;

public:

  //! to have an "end" iterator
  _PreOrderIterator() {}

  _PreOrderIterator(NodeType* _node): node(_node), child_index() 
  {
    child_index.push(0);
  }

  //! increment operator
  _PreOrderIterator& operator++() {
    if(!child_index.empty()){
      if(node->isLeaf()){
        while(++child_index.top() >= node->getNumberOfSons()){
          child_index.pop();
          if(child_index.empty()) return *this;
          node = node->getFather();
        }
        node = node->getSon(child_index.top());
      } else {
        node = node->getSon(child_index.top());
        child_index.push(node->getNumberOfSons());
      }
    }// if the stack is non-empty
    return *this; 
  }
  //! dereferece
  NodeType& operator*() { return *node; }
  NodeType* operator->() { return node; }
  //! comparison
  bool operator==(const _PreOrderIterator& it) const
  {
    return (child_index.size() == it.child_index.size()) &&
           (node == it.node) &&
           (child_index == it.child_index);
  }
  bool operator!=(const _PreOrderIterator& it) const
  {
    return !operator==(it);
  }

};


//! an iterator visiting the tree in post order
template<class NodeType>
class _PostOrderIterator
{
private:
  NodeType* node;
  stack<unsigned> child_index;

  //! whenever we see the root of a new subtree, call this to go down to the first leaf
  void reached_new_subtree()
  {
    while(!node->isLeaf()){
      node = node->getSon(0);
      child_index.push(0);
    }
  }
public:

  //! to have an "end" iterator
  _PostOrderIterator() {}

  _PostOrderIterator(NodeType* _node): node(_node), child_index() 
  {
    child_index.push(0);
    reached_new_subtree();
  }

  bool is_valid() const
  {
    return !child_index.empty();
  }

  //! get the depth of the current vertex
  unsigned current_depth() const
  {
    return child_index.size() - 1;
  }

  //! increment operator
  _PostOrderIterator& operator++() 
  {
    if(!child_index.empty()){
      // take the current nodes child index off the stack and go up to its parent
      child_index.pop();
      if(child_index.empty()) return *this;
      node = node->getFather();
      // if there are unexplored children, go there
      if(++child_index.top() < node->getNumberOfSons()){
        node = node->getSon(child_index.top());
        child_index.push(0);
        reached_new_subtree();
      }
    }// if the stack is non-empty
    return *this; 
  }
  //! dereferece
  NodeType& operator*() { return *node; }
  NodeType* operator->() { return node; }
  //! comparison
  bool operator==(const _PostOrderIterator& it) const
  {
    // to compare correctly to end(), we consider two iterators the same if they both have an empty child stack
    if(child_index.empty() && it.child_index.empty())
      return true;
    else
      return (child_index.size() == it.child_index.size()) &&
             (node == it.node) &&
             (child_index == it.child_index);
  }
  bool operator!=(const _PostOrderIterator& it) const
  {
    return !operator==(it);
  }

};


using PreOrderIterator = _PreOrderIterator<MyNode>;
using PreOrderConstIterator = _PreOrderIterator<const MyNode>;
using PostOrderIterator = _PostOrderIterator<MyNode>;
using PostOrderConstIterator = _PostOrderIterator<const MyNode>;


template<class NodeType, class Iter>
struct _Traversal
{
  NodeType* const root;

  _Traversal(NodeType* const _r): root(_r) {}
  Iter begin() const { return Iter(root); }
  Iter end() const { return Iter(); }
};

using PreOrderTraversal = _Traversal<MyNode, PreOrderIterator>;
using PreOrderConstTraversal = const _Traversal<const MyNode, PreOrderConstIterator>;
using PostOrderTraversal = _Traversal<MyNode, PostOrderIterator>;
using PostOrderConstTraversal = const _Traversal<const MyNode, PostOrderConstIterator>;





