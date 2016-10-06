

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


typedef unsigned StId;


struct Clade: public pair<unsigned, unsigned>
{
  using pair<unsigned, unsigned>::pair;
  Clade(): pair() {}
  Clade(const unsigned x): pair(x,x) {}
  // return whether the given clade is a subclade of us
  bool contains(const Clade& _clade) const
  {
    return (first <= _clade.first) && (_clade.second <= second);
  }
  operator unsigned() const
  {
    assert(first == second);
    return first;
  }
  unsigned size() const
  {
    assert(second >= first);
    return (second - first) + 1;
  }
};

class NodeInfos  {
protected:
		
  unsigned centroidPathNumber;
  unsigned depth; 
  unsigned numberOfDescendents;
  StId stId;
  bool isVisited;
  unsigned preorder;
  Clade clade;

public:

      void setIsVisisted(bool vis){
    isVisited = vis;
  }
  
  bool getIsVisisted() const{
    return isVisited ;
  }


      void setCentroidPathNumber(int id){
    centroidPathNumber = id;
  }
  
  unsigned getCentroidPathNumber() const{
    return centroidPathNumber ;
  }

  void setPreOrder(int id){
    preorder = id;
  }
  
  int getPreOrder() const{
    return preorder ;
  }
  
    
  void setStId(int id){
    stId = id;
  }
  
  int getStId() const {
    return stId ;
  }
  

      void setDepth(int id){
    depth = id;
  }
  
  int getDepth() const{
    return depth ;
  }
  
  void setNumberOfDescendents(int id){
    numberOfDescendents = id;
  }
  
  int getNumberOfDescendents() const{
    return numberOfDescendents;
  }
  
  void setClade(const Clade& _clade){
    clade = _clade;
  }
  
  const Clade& getClade() const{
    return clade;
  }
  
  
};

#endif /*NODEINFOS_H_*/

