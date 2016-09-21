
#ifndef _HSLINKERN_HPP
#define _HSLINKERN_HPP

#include "db.hpp"
#include "Hypergraph.hpp"

struct petals_less_than {
  bool yes;
  unsigned int k;

  petals_less_than(int _k) : k(_k) { yes = true; }
  void operator()(const Edge &c) {
    if(lookup(c)->petals > k) yes = false;
  }
};

struct add_petal {
  const Edge &e;

  add_petal(const Edge &_e) : e(_e) {}

  void operator()(const Edge &c) {
    Leaf *l(lookup(c));
    Edge enotc(edge_subtraction(e, c));
    
    for(Edge::const_iterator i = enotc.begin(); i != enotc.end(); ++i)
      if(l->used(*i)) return; // l->forbidden = true;

    //    if(l->forbidden) return; // XXX what is this for?

    l->petals++;

    // cout << "petals " << l->petals << " for ";
    // copy(c.begin(), c.end(), ostream_iterator<Vertex>(cout, " "));
    // cout << "(added ";
    // copy(e.begin(), e.end(), ostream_iterator<Vertex>(cout, " "));
    // cout << ") and enotc is ";
    // copy(enotc.begin(), enotc.end(), ostream_iterator<Vertex>(cout, " "));
    // cout << endl;

    //    l->forbidden = true;

    for(Edge::const_iterator i = enotc.begin(); i != enotc.end(); ++i)
      l->setused(*i);
  }
};

Hypergraph kernelize(const Hypergraph &G, const Graphstats &s, unsigned int k);
void print_stats(const Graphstats &stats);

#endif
