

#ifndef NODEINFOS_H_
#define NODEINFOS_H_


// From the STL:
#include <string>
#include <vector>


#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/NodeTemplate.h>
#include <Bpp/Phyl/TreeExceptions.h>

using namespace std;
using namespace bpp;

// abbreviations to increase readability
#define clade(x) (x)->getInfos().clade
#define stid(x) (x)->getInfos().stId
#define same_centroid_path(x,y) ((x)->getInfos().cp_num == (y)->getInfos().cp_num)
#define is_ancestor(x,y) ((x)->getInfos().clade.contains((y)->getInfos().clade))

typedef unsigned StId;

class NodeInfos;
typedef NodeTemplate<NodeInfos> MyNode;



struct Clade: public pair<size_t, size_t>
{
  using pair<size_t, size_t>::pair;
  Clade(): pair() {}
  Clade(const size_t x): pair(x,x) {}
  // return whether the given clade is a subclade of us
  bool contains(const Clade& _clade) const
  {
    return (first <= _clade.first) && (_clade.second <= second);
  }
  // return the number of leaves in the clade
  size_t size() const
  {
    assert(second >= first);
    return (second - first) + 1;
  }
};

ostream& operator<<(ostream& os, const Clade& c)
{
  //return os << reinterpret_cast<const pair<unsigned, unsigned>&>(c); // why does this segfault?
  return os << "(" << c.first<<", "<<c.second<<")";
}

// return a minimum clade spanning the given clades
Clade get_spanning_clade(const Clade& c1, const Clade& c2)
{
  return Clade(std::min(c1.first, c2.first), std::max(c1.second, c2.second));
}


struct NodeInfos
{
  StId stId;
  unsigned cp_num;        // the centroid path number
  unsigned depth; 
  unsigned subtree_size;  // the size of the subtree rooted at this node
  unsigned po_num;        // the nodes preorder number
  Clade clade;            // the clade below the node

  // return the leaf_po number, in case we are a leaf
  unsigned leaf_po_num() const
  {
    assert(clade.first == clade.second);
    return clade.first;
  }

};


// get the child of v that is not in the same centroid path as v
const MyNode* get_non_cp_child(const MyNode* const v)
{
  assert(!v->isLeaf());
  if(v->getInfos().cp_num != v->getSon(0)->getInfos().cp_num)
    return v->getSon(0);
  else
    return v->getSon(1);
}

// compute the size of the subtree rooted at the child of v that is not in the same centroid path as v
unsigned compute_nj(const MyNode* v)
{
  return v->isLeaf() ? 1 : get_non_cp_child(v)->getInfos().subtree_size;
}


#endif /*NODEINFOS_H_*/

