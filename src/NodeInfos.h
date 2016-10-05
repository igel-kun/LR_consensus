

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

class NodeInfos  {
	 typedef NodeTemplate<NodeInfos> MyNode;
		protected:
		
		int centroidPathNumber;
		int depth; 
		int numberOfDescendents;
		int stId;
		bool isVisited;
		int preorder;
		pair<int, int> clade;

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
		
		int getCentroidPathNumber() const{
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
		
		void setClade(pair <int, int> inds){
			clade=inds;
		}
		
		pair <int, int> getClade(){
			return clade;
		}
		
		
};

#endif /*NODEINFOS_H_*/

