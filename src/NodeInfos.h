

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
		
		int centroidPathNumber;
		int depth; 
		int numberOfDescendents;
		int stId;

		//not used for now

	    vector<int> clade;
		vector<string> cladeNames;
		vector <MyNode*>  cladeNode;
		vector <MyNode*>  cladeNodeAbove;
		vector <MyNode*> cladeNodeAboveDuplications;
		int counter;
		bool square;	
		bool visited;	

	public:

        void setCentroidPathNumber(int id){
			centroidPathNumber = id;
		}
		
		int getCentroidPathNumber(){
			return centroidPathNumber ;
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
		
		int getDepth(){
			return depth ;
		}
		
		void setNumberOfDescendents(int id){
			numberOfDescendents = id;
		}
		
		int getNumberOfDescendents(){
			return numberOfDescendents;
		}
	
		//not used for now		

		void copyInClade(vector <int> leaves){
			(clade).clear(); 
			for( unsigned int i = 0; i <  leaves.size(); i++ ){
     				(clade).push_back( leaves[i] );
			}
		}
		
		void copyInCladeNames(vector <string> leaves){
			(cladeNames).clear(); 
			for( unsigned int i = 0; i <  leaves.size(); i++ ){
     				(cladeNames).push_back( leaves[i] );
			}
		}
		
		
		void copyInCladeNode(vector <MyNode *> leaves){
			(cladeNode).clear();
			for( unsigned int i = 0; i <  leaves.size(); i++ ){
     				(cladeNode).push_back(  leaves[i] );
			}
		}
		
		void copyInCladeNodeAbove(vector <MyNode *> leaves){
			(cladeNodeAbove).clear();
			for( unsigned int i = 0; i <  leaves.size(); i++ ){
     				(cladeNodeAbove).push_back(  leaves[i] );
			}
		}	
		
		void copyInCladeNodeAboveDuplications(vector <MyNode *> leaves){
			(cladeNodeAboveDuplications).clear();
			for( unsigned int i = 0; i <  leaves.size(); i++ ){
     				(cladeNodeAboveDuplications).push_back(  leaves[i] );
			}
		}	
		
		void copyInCladeNodeAboveNoClearDuplications(vector <MyNode *> leaves){
			for( unsigned int i = 0; i <  leaves.size(); i++ )
     			(cladeNodeAboveDuplications).push_back(  leaves[i] );
		}	
			
		void copyInCladeNodeAboveNoClear(vector <MyNode *> leaves){
			for( unsigned int i = 0; i <  leaves.size(); i++ )
     			(cladeNodeAbove).push_back(  leaves[i] );
		}	
		
		void setCladeNodeAbove(vector <MyNode *> leaves){
			cladeNodeAbove = leaves;
		}
		
		void setCladeNodeAboveDuplications(vector <MyNode *> leaves){
			cladeNodeAboveDuplications = leaves;
		}
		
		vector <MyNode * >getCladeNode(){
			return cladeNode;
		} 
		vector <MyNode * >getCladeNodeAbove(){
			return cladeNodeAbove;
		} 
		
		vector <MyNode * >getCladeNodeAboveDuplications(){
			return cladeNodeAboveDuplications;
		} 
		
		vector<int> getClade(){
			return clade;
		};
		
		vector<string> getCladeNames(){
			return cladeNames;
		};
	
		
		void eraseCladeNode(){
			int j= cladeNode.size();
			for(int i = 0; i <  j; i++ ){
     				(cladeNode).pop_back();
			}
		}
		
		void eraseCladeNodeAbove(){	
			int j = cladeNodeAbove.size();
			for(int i = 0; i <  j; i++ ){
     				(cladeNodeAbove).pop_back();
			}
		}
		
		void eraseCladeNodeAboveDuplications(){	
			int j = cladeNodeAboveDuplications.size();
			for(int i = 0; i <  j; i++ ){
     				(cladeNodeAboveDuplications).pop_back();
			}
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
			(cladeNode).push_back(  leaves );
		}
		
		void copyInCladeNodeAboveLinear(MyNode *leaves){
			//(cladeNodeAbove).clear();
			(cladeNodeAbove).push_back(  leaves );		
		}	

		
		
};

#endif /*NODEINFOS_H_*/

