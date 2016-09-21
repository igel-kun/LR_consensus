#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <cassert>
#include "Hypergraph.hpp"
#include "bitops.hpp"

using namespace std;

Hypergraph read_hypergraph(const char *file) {
  Hypergraph G;
  ifstream input(file, ifstream::in);
  string line;

  getline(input, line);
  while(input.good()) {
    istringstream linestream(line);
    Edge e;
    Vertex x;
    Vertex y = 0;
    bool sorted = 1;
    while(linestream >> x) {
      e.push_back(x);
      if(x < y) sorted = 0;
      y = x;
    }
    if(!sorted) sort(e.begin(), e.end());
    G.push_back(e);
    getline(input, line);
  }
  return G;
}

Graphstats get_stats(const Hypergraph &G) {
  Graphstats stats;
  stats.edges = G.size();
  stats.edgesize = 0;
  stats.edgesums = 0;
  stats.vertices = 0;
  stats.hugesize = 0;
  stats.hugeedges = 0;
  
  for(Hypergraph::const_iterator i = G.begin(); i != G.end(); ++i) {
    stats.edgesize = max(stats.edgesize, i->size());
    stats.edgesums += i->size();
    if((!ls(i->size())) || (ls(i->size()) > stats.edges)) {
      stats.hugeedges++;
      if(stats.hugesize == 0) stats.hugesize = i->size();
      else stats.hugesize = min(stats.hugesize, i->size());
    }
    for(Edge::const_iterator j = i->begin(); j != i->end(); ++j)
      stats.vertices = max(stats.vertices, *j);
  }
  if(stats.hugesize==0) stats.hugesize=stats.edgesize + 1;
  return stats;
}

void bar(int have, int max, char c) {
  int barsize = (70*have)/max;
  for(int i = 1; i <= barsize+1; ++i)
    cout << c;
  cout << endl;
}

// sort edges of a hypergraph by size, return edge statistics
void edges_by_size(Hypergraph &G, const Graphstats &s) {
  vector<Hypergraph> edges(s.edgesize + 1);

  Hypergraph::iterator i = G.begin();
  while(i != G.end()) {
    edges[i->size()].push_back(*i);
    i = G.erase(i);
  }

  assert(G.empty());

  cout << "Percentage of edges of size:  (\"=\" small edges, \"+\" large edges)" << endl;
  for(Hypergraph::size_type i = 0; i < edges.size(); ++i) {
    unsigned size(edges[i].size());
    if(size > 0) {
      cout << setw(3) << i << ": ";
      if(i < s.hugesize) bar(size, s.edges, '=');
      else bar(size, s.edges, '+');
    }
    
    Hypergraph::iterator j = edges[i].begin();
    while(j != edges[i].end()) {
      G.push_back(*j);
      j = edges[i].erase(j);
    }
  }
}

Hypergraph approx_hs(const Hypergraph& G, const Graphstats &s) {
  vector<int> dead(s.vertices + 1);
  Hypergraph apx;
  
  for(Hypergraph::const_iterator j = G.begin(); j != G.end(); ++j) {
    bool skip = 0;
    for(Edge::const_iterator k = j->begin(); k != j->end(); ++k) if(dead[*k]) skip = 1;
    if(!skip) {
      apx.push_back(*j);
      for(Edge::const_iterator k = j->begin(); k != j->end(); ++k) dead[*k] = 1;
    }
  }
  return apx;
}

class contains_vertex {
private:
  Vertex v;
public:
  contains_vertex(Vertex _v) : v(_v) {}
  bool operator() (const Edge &e) {
    return (e.end() != find(e.begin(), e.end(), v));
  }
};

int greedy_hs(const Hypergraph &G, const Graphstats &s) {
  int res(0);
  vector<int> deg(s.vertices+1);
  vector<bool> hit(s.vertices+1);

  for(Hypergraph::const_iterator i = G.begin(); i != G.end(); ++i)
    for(Edge::const_iterator j = i->begin(); j != i->end(); ++j)
      deg[*j]++;
  
  for(Hypergraph::const_iterator i = G.begin(); i != G.end(); ++i){
    Vertex maxdeg(0);
    bool ishit = false;
    for(Edge::const_iterator j = i->begin(); j != i->end(); ++j) {
      if(hit[*j]) { ishit = true; break; }
      if(deg[*j] > deg[maxdeg]) maxdeg = *j;
    }
    if(!ishit) { 
      hit[maxdeg] = true;
      res++;
    }
  }
  
  return res;
}

int searchtree_hs(const Hypergraph &G, const Graphstats &s, vector<bool>& hit, int k) {
  if(k<0) return -1; // solution is larger than k
  for(Hypergraph::const_iterator i = G.begin(); i != G.end(); ++i){
    bool ishit = false;
    int minsol = s.vertices + 1;
    for(Edge::const_iterator j = i->begin(); j != i->end(); ++j)
      if(hit[*j]) { ishit = true; break; }
    
    if(!ishit) {
      for(Edge::const_iterator j = i->begin(); j != i->end(); ++j) {
	hit[*j] = true;
	int solsize = searchtree_hs(G,s,hit,k-1);
	hit[*j] = false;
	if(solsize < minsol) minsol = solsize;
      }
      return minsol + 1;
    }
  }
  return(0); // here we land if everything was hit
}


Edge get_subedge(const Edge& e, unsigned long mask) {
  Edge res;
  res.reserve(e.size());
  for(Edge::size_type i = 0; i < e.size(); ++i)
    if(ls(i) & mask) res.push_back(e[i]);
  return res;
}

Edge edge_intersection(const Edge& e1, const Edge& e2) {
  // requires that edges are already sorted!
  int length = min(e1.size(), e2.size());
  Edge cap;
  cap.reserve(length);
  int i1 = 0;
  int i2 = 0;
  while(i1 < length && i2 < length) {
    if(e1[i1] < e2[i2]) i1++;
    else if(e1[i1] > e2[i2]) i2++;
    else {
      cap.push_back(e1[i1]);
      i1++;
      i2++;
    }
  }
  return cap;
}

Edge edge_subtraction(const Edge& e1, const Edge& e2) {
  if(e2.empty()) return e1;
  
  Edge sub;
  sub.reserve(e1.size());
  Edge::size_type i1 = 0;
  Edge::size_type i2 = 0;
  while(i1 < e1.size() && i2 < e2.size()) {
    if(e1[i1] < e2[i2]) {
      sub.push_back(e1[i1]);
      i1++;
    } else {
      if(e1[i1] > e2[i2]) i2++;
      else {
	i1++;
	i2++;
      }
    }
  }

  if(i2 >= e2.size()) {
    while(i1 < e1.size())
      sub.push_back(e1[i1++]);
  }
  return sub;
}

void print_edge(const Edge& e) {
  copy(e.begin(), e.end(), ostream_iterator<Vertex>(cout, " "));
  cout <<endl;
}
