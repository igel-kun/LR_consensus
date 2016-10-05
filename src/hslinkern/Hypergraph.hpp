#ifndef _HYPERGRAPH_H
#define _HYPERGRAPH_H

#include <vector>
#include <list>
#include <cassert>
#include "bitops.hpp"
using namespace std;

typedef unsigned Vertex;
typedef vector<Vertex> Edge;
typedef list<Edge> Hypergraph;

struct Graphstats {
  unsigned vertices;
  unsigned edgesums;
  unsigned hugeedges;
  Edge::size_type hugesize;
  Hypergraph::size_type edges;
  Edge::size_type edgesize;
};

Hypergraph read_hypergraph(const char*);
Graphstats get_stats(const Hypergraph&);
void edges_by_size(Hypergraph &, const Graphstats&);
Hypergraph approx_hs(const Hypergraph&, const Graphstats&);
Edge edge_intersection(const Edge&, const Edge&);
Edge edge_subtraction(const Edge&, const Edge&);
void print_edge(const Edge&);
Edge get_subedge(const Edge& e, unsigned long mask);
int greedy_hs(const Hypergraph &, const Graphstats &);
int searchtree_hs(const Hypergraph &G, const Graphstats &s, vector<bool>& hit, int k);

template<class F>
void foreach_intersection(const Hypergraph &G, const Edge &e, const Graphstats &s, F &f) {
  if(e.size() >= s.hugesize) {
    for(Hypergraph::const_iterator j = G.begin(); j != G.end(); ++j)
      f(edge_intersection(e, *j));
  } else {
    unsigned long max = ls(e.size());
    assert(max);
    for(unsigned int j = 0; j < max; ++j)
      f(get_subedge(e, j));
  }
}

#endif
