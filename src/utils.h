

// File: utils.h
// Created by: Celine Scornavacca 
// Created on: Mon Oct  21 11:00 2013


#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <queue>

/**
 * @brief This function collapses the edge having the node on as end.
 * @param tree The tree to modify.
 * @param on The node at the lower estremity of the edge to collapse.
 */

#define update(x,y) x = std::max(x, y)

//output freakin pairs
template<class A, class B>
std::ostream& operator<<(std::ostream& os, const std::pair<A,B>& e)
{
  return os << "(" << e.first << ", "<< e.second << ")";
}

//output freakin vectors
template<class E>
std::ostream& operator<<(std::ostream& os, const std::vector<E>& e)
{
  os << "(";
  for(const auto& item: e) os << item << " ";
  return os << ")";
}

//output freakin lists
template<class E>
std::ostream& operator<<(std::ostream& os, const std::list<E>& e)
{
  os << "(";
  for(const auto& item: e) os << item << " ";
  return os << ")";
}

//output freakin sets
template<class E>
std::ostream& operator<<(std::ostream& os, const std::set<E>& e)
{
  os << "(";
  for(const auto& item: e) os << item << " ";
  return os << ")";
}





    /**
     * @brief This function returns the union of two ordered vectors.
     * @param a First ordered vector.
     * @param b Second ordered vector.
     * @return A vector which is the ordered union of a and b.
     */
    
    template<class T>
    std::vector<T> unionVector(const std::vector<T> & a, const std::vector<T> & b){
        std::vector<T> unione;
        unsigned int i=0;
        unsigned int j=0;
        unsigned int A =  a.size();
        unsigned int B =  b.size();
        while((i != A)&&(j != B)){
            if(a[i]>b[j]){
                unione.push_back(b[j]);
                j++;
            }
            else if(a[i]<b[j]){
                unione.push_back(a[i]);
                i++;
            }
            else if(a[i]==b[j]){
                unione.push_back(a[i]);
                i++;
                j++;
            }
        }
        if(i==A){
            for (unsigned int k=j;k<b.size();k++){
                unione.push_back(b[k]);
            }
        }
        else if(j==B){
            for (unsigned int k=i;k<a.size();k++){
                unione.push_back(a[k]);
            }
        }
        return unione;      
}


    /**
     * @brief This function returns the union of all ordered std::vectors given in input.
     * @param a An array of ordered std::vectors.
     * @param b The size of the array.
     * @return A std::vector which is the ordered union of all ordered std::vectors given in input.
     */
     
    /* remove, because the version below is more efficient
    template<class T>
    std::vector<T> allLeaves(const std::vector<T> a[], unsigned size){
      if(size == 1) return a[0];
      std::vector<T> all = a[0];
      for(int i=1; i< size; ++i){
        const unsigned all_size = all.size();
        copy(a[i].begin(), a[i].end(), back_inserter(all));
        inplace_merge(all.begin(), all.begin() + all_size, all.end());
      }
      return all;
    }
    */
    template<class Container>
    struct iter_cmp{
      typedef typename Container::const_iterator it;
      bool operator()(const std::pair<it,it>& x, const std::pair<it,it>& y) const { return *x.first > *y.first; }
    };
    
    template<class T>
    std::vector<T> allLeaves(const std::vector<T> a[], unsigned size){
      // in the std::queue, we map the smallest element of each std::vector to its std::vector
      typedef typename std::vector<T>::const_iterator T_iter;
      typedef std::pair<T_iter, T_iter> T_iterpair;
      std::priority_queue<T_iterpair, std::vector<T_iterpair>, iter_cmp<std::vector<T> > > smallest_items;

      unsigned total_size = 0;
      for(unsigned i = 0; i < size; ++i){
        total_size += a[i].size();
        smallest_items.emplace(a[i].begin(), a[i].end());
      }

      std::vector<T> result;
      result.reserve(total_size);
      while(!smallest_items.empty()){
        T_iterpair smallest = smallest_items.top();
        smallest_items.pop();
        // add the item to the result avoiding duplicates
        if(*smallest.first != *result.rbegin())
          result.push_back(*smallest.first);
        // see if we're at the end of the corresponding std::vector
        if(++smallest.first != smallest.second)
          smallest_items.emplace(smallest);
      }
      return result;
    }


        /**
     * @brief This function returns the intersection of two ordered std::vectors.
     * @param a First ordered std::vector.
     * @param b Second ordered std::vector.
     * @return A std::vector which is the ordered intersection of a and b.
     */
  template<class T>
  vector<T> intersection(vector<T> a, vector<T> b){
		vector<T> intersection;
		unsigned int i=0;
		unsigned int j=0;
		unsigned int A =  a.size();
		unsigned int B =  b.size();
		while((i != A)&&(j != B)){
			if(a[i]>b[j]){
				j++;
			}
			else if(a[i]<b[j]){
				i++;
			}
			else if(a[i]==b[j]){
				intersection.push_back(a[i]);
				i++;
				j++;
			}
		}
		return intersection;		
	}



    template<class T>
    bool binarySearchOrdered(const std::vector < T> & a,T b){    
        int fine = a.size();
        int p,u,m;
        p = 0;
        u = fine-1;
        while(p<=u) {
            m = (p+u)/2;
            if ( a[m]==b) return true; 
            if( a[m]<b)
            p = m+1;
            else
            u = m-1;
        }
        return false;
    } ;


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

};

/* this function reads a list of trees written in a file (path) in a newich format, separed by semicolons and returns a list of MyTree*/
vector < MyTree *>  readTrees(const string & path) throw (Exception) {
    // Checking the existence of specified file
    
    ifstream file(path.c_str(), ios::in);
    //if (! file) { throw IOException ("\nError reading file.\nInvalid options!\nUsage:\n ./physic -s sourceTreeFile -t threshold.\nwhere:\n - sourceTreeFile contains a set of rooted trees in newick format with bootstrap values and possibly edge lengths.\n - threshold indicates bootstrap values under which clades are not considered for building the supertree\n(typically a threshold of 70 can be used when source trees where obtained from 100 bootstrap replicates).\n"); }
    if (! file) { throw IOException ("\nError reading file.\n"); }
    
    vector<MyTree*> trees;
    string temp, description;
    while (! file.eof()) {
        temp = FileTools::getNextLine(file);
        if(temp.size() != 0){
            string::size_type index = temp.find(";");
            if(index == string::npos) throw Exception("readTrees(). Bad format: no semi-colon found.");
            if(index < temp.size()) {
                description += temp.substr(0, index + 1);   
                TreeTemplate<Node> * tree = TreeTemplateTools::parenthesisToTree(description,true);
                trees.push_back(new MyTree((MyNode*)(tree->getRootNode())));
                delete tree;
                description = temp.substr(index);   
            } else description += temp;
        }
    }
    file.close();
    return trees;   
};



#endif /*UTILS_H_*/


