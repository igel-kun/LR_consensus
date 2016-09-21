

#ifndef NODEINFOS_H_
#define NODEINFOS_H_


// From the STL:
#include <string>
#include <vector>


#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/NodeTemplate.h>
#include <Bpp/Phyl/TreeExceptions.h>

using namespace std;

class NodeInfos  {
  typedef NodeTemplate<NodeInfos> MyNode;
		
  protected:
  vector<int> clade;
  vector<string> cladeNames;
  vector<MyNode*>  cladeNode;
  vector<MyNode*>  cladeNodeAbove;
  vector<MyNode*> cladeNodeAboveDuplications;
  int counter;
  bool square;	
  bool visited;
  int stId;

  public:

	void copyInClade(const vector<int>& leaves)
  {
		clade = leaves;
	}

  void copyInCladeNames(const vector<string>& leaves)
  {
    cladeNames = leaves;
  }
	
  void copyInCladeNode(const vector<MyNode*>& leaves)
  {
    cladeNode = leaves;
	}
		
	void copyInCladeNodeAbove(const vector<MyNode*>& leaves)
  {
    cladeNodeAbove = leaves;
  }	
		
	void copyInCladeNodeAboveDuplications(const vector<MyNode*>& leaves)
  {
    cladeNodeAboveDuplications = leaves;
	}	
		
	void copyInCladeNodeAboveNoClearDuplications(const vector<MyNode*>& leaves)
  {
    copy(leaves.begin(), leaves.end(), back_inserter(cladeNodeAboveDuplications));
	}	
			
	void copyInCladeNodeAboveNoClear(const vector<MyNode*>& leaves)
  {
    copy(leaves.begin(), leaves.end(), back_inserter(cladeNodeAbove));
  }	
		
	void setCladeNodeAbove(const vector<MyNode*>& leaves)
  {
			cladeNodeAbove = leaves;
	}
		
	void setCladeNodeAboveDuplications(const vector<MyNode*>& leaves){
	  cladeNodeAboveDuplications = leaves;
	}
		
	const vector<MyNode*>& getCladeNode() const
  {
	  return cladeNode;
	}

  const vector<MyNode*>& getCladeNodeAbove() const
  {
	  return cladeNodeAbove;
	}

	const vector<MyNode*>& getCladeNodeAboveDuplications() const
  {
	  return cladeNodeAboveDuplications;
	} 
		
	const vector<int>& getClade() const
  {
	  return clade;
	}
		
	const vector<string>& getCladeNames() const
  {
	  return cladeNames;
	}
	
	void setStid(int id)
  {
	  stId = id;
	}
		
	int getStid()
  {
	  return stId ;
	}
		
		
	void eraseCladeNode()
  {
	  cladeNode.clear();
	}
		
	void eraseCladeNodeAbove()
  {
    cladeNodeAbove.clear();
  }
		
	void eraseCladeNodeAboveDuplications()
  {
    cladeNodeAboveDuplications.clear();
  }
					
  void setVisited(bool a) {visited = a;};
  bool getVisited(){return visited;};
  
  void setSquare(bool a) { square = a;};
  bool getSquare(){return square;};
  
  void setCounter(int a) {counter= a;};
  void addCounter() {counter ++;};
  int getCounter(){return counter;};	
  
  void copyInCladeNodeLinear(MyNode *leaves){
    //(cladeNode).clear();
    cladeNode.push_back(  leaves );
  }
  
  void copyInCladeNodeAboveLinear(MyNode *leaves){
    //(cladeNodeAbove).clear();
    cladeNodeAbove.push_back(  leaves );		
  }	
		
};

#endif /*NODEINFOS_H_*/


