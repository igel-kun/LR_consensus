#ifndef MATRIXTRIPLETS__H_
#define MATRIXTRIPLETS__H_

#include <unordered_set>

#include "vector2d.h"

#include "NodeInfos.h"

using namespace bpp;

//! a 2d vector assigning a pair of PREORDER NUMBERS to a set of STIDs such that the three leaves form a conflict
typedef symmetric_vector2d<unordered_set<unsigned> > MatrixTriplets;

void add_triple(MatrixTriplets& trip, const unsigned x, const unsigned y, const unsigned z)
{
  trip[{x, y}].insert(z);
}

bool is_not_triple(const MatrixTriplets& trip, const unsigned x, const unsigned y, const unsigned z)
{
  const auto& conflict_set = trip[{x ,y}];
  return conflict_set.find(z) == conflict_set.end();
}

/*

#define  maxDim 1000
class MatrixTriplets
{	
	unsigned dim;
	unordered_set< StId > ** GroupOfTriplets ;
	
	public:

	void setDim(unsigned dimTrees){
		dim= dimTrees;
		GroupOfTriplets = new unordered_set< StId >*[dimTrees];
		for (unsigned i=0;i< dimTrees;i++){
			GroupOfTriplets[i] = new unordered_set< StId >[dimTrees];
		}
	}
	
	void eraseValues(){
		for (unsigned i=0;i< dim;i++){
			for (unsigned j=0;j< dim;j++){
				GroupOfTriplets[i][j].clear();		
			}
		}
	}
		 
	void deleteMatrix(){
		for (unsigned i=0;i< dim;i++){		
			delete [] (GroupOfTriplets[i]);		
		}	
		delete [] (GroupOfTriplets);			
	}
	
	unsigned getDim(){
		return dim;
	}

	unordered_set< StId > & getSet (unsigned i, unsigned j){
		if(i<j)
			return GroupOfTriplets[i][j];
		else
			return GroupOfTriplets[j][i];
	}
		
	bool isTriplet(StId i, StId j, StId z){
		if(i<=j){
			return GroupOfTriplets[i][j].find(z)!=GroupOfTriplets[i][j].end();
		}
		else
			return GroupOfTriplets[j][i].find(z)!=GroupOfTriplets[j][i].end();
	}
	
	void add(MyNode& I, MyNode& J, MyNode& Z) {
		//cout << I.getName()  << " "<< J.getName()  << " " << Z.getName() << " " << endl;
		StId i= stid(&I);
		StId j= stid(&J);
		StId z= stid(&Z);
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

*/

#endif /*MATRIXTRIPLETS_H_*/
