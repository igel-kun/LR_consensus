#ifndef MATRIXTRIPLETS__H_
#define MATRIXTRIPLETS__H_

using namespace bpp;
#include <unordered_set>
typedef NodeTemplate<NodeInfos> MyNode;

#define  maxDim 1000
class MatrixTriplets
{	
	private:
	int dim;
	unordered_set< int > ** GroupOfTriplets ;
	
	public:
	


	void setDim(int dimTrees){
		dim= dimTrees;
		GroupOfTriplets = new unordered_set< int >*[dimTrees];
		for (int i=0;i< dimTrees;i++){
			GroupOfTriplets[i] = new unordered_set< int >[dimTrees];
		}
	}
	
	void eraseValues(){
		for (int i=0;i< dim;i++){
			for (int j=0;j< dim;j++){
				GroupOfTriplets[i][j].clear();		
			}
		}
	}
		 
	void deleteMatrix(){
		for (int i=0;i< dim;i++){		
			delete [] (GroupOfTriplets[i]);		
		}	
		delete [] (GroupOfTriplets);			
	}
	
	int getDim(){
		return dim;
	}
	
	bool getValue (int i,int j,int z){
		if(i<=j){
			return GroupOfTriplets[i][j].find(z)!=GroupOfTriplets[i][j].end();
		}
		else
			return GroupOfTriplets[j][i].find(z)!=GroupOfTriplets[j][i].end();
	}
	
	void add(MyNode  & I,MyNode  & J,MyNode  & Z) {
		//cout << I.getName()  << " "<< J.getName()  << " " << Z.getName() << " " << endl;
		int i= (I.getInfos()).getStId();
		int j=  (J.getInfos()).getStId();
		int z= ( Z.getInfos()).getStId();
		//cout << i << " " << j << " " << z << endl;
		if (i<j) 
			GroupOfTriplets[i][j].insert(z);
		else 
			GroupOfTriplets[j][i].insert(z);
	}
	
// 	void print(map<string,int> association, char* out ){
// 		ofstream list(out, 1 ? (ios::out) : (ios::out|ios::app));
// 		map<string,int>::iterator iter;
// 		iter = association.begin();
// 		for (int i=0;i< dim;i++){
// 			for (int j=i;j< dim;j++){
// 				for (int z=0;z< dim;z++){
// 					if(GroupOfTriplets[i][j][z] !=0){
// 						iter = association.begin();
// 						for(int h=0;h< i;h++)
// 							(iter++);
// 						list << (iter++)->first << ",";	
// 						iter = association.begin();	
// 						for(int h=0;h< j;h++)
// 							(iter++);
// 						list << (iter++)->first << "|";
// 						iter = association.begin();	
// 						for(int h=0;h< z;h++)
// 							(iter++);	
// 						list << (iter++)->first ;
// 						list << " = " << GroupOfTriplets[i][j][z] << endl;	
// 					}
// 				}
// 			}
// 		}
// 	}
// 	
		
};



#endif /*MATRIXTRIPLETS_H_*/
