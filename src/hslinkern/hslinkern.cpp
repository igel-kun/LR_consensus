
#include <iostream>
#include "hslinkern.hpp"

using namespace std;

Hypergraph kernelize(const Hypergraph &G, const Graphstats &s, unsigned int k) {
  Hypergraph H;

//  cout << endl << "... computing kernel for k = " << k << " ..." << endl;
  
  initialize_db(G, s);
  
  for(Hypergraph::const_iterator i = G.begin(); i != G.end(); ++i) {
    petals_less_than lt(k);
    foreach_intersection(G, *i, s, lt);
    if(lt.yes) {
      add_petal addp(*i);
      H.push_back(*i);
      foreach_intersection(G, *i, s, addp);
      //      foreach_intersection(G, *i, s, unforbid);
      // cout << "petals " << k+1 << " for ";
      // copy(i->begin(), i->end(), ostream_iterator<Vertex>(cout, " "));
      // cout << endl;
      lookup(*i)->petals = k+1;
    // } else {
    //   cout << "threw away ";
    //   copy(i->begin(), i->end(), ostream_iterator<Vertex>(cout, " "));
    //   cout << endl;
    }

    if(lookup(Edge())->petals > k) {
      cout << "No-instance" << endl;
      break;
    }
  }
  return H;
}

void print_stats(const Graphstats &stats) {
  cout << "----------------" <<endl;
  cout << "Edges              : " << stats.edges << endl;
  cout << "Maximum edge size  : " << stats.edgesize << endl;
  cout << "Smallest huge edge : " << stats.hugesize << endl;
  cout << "Huge edges         : " << stats.hugeedges << endl;
}

