#ifndef _TRIE_HPP
#define _TRIE_HPP

#include <Hypergraph.hpp>

class Leaf {
public:
  Leaf() : petals(0) {}
  unsigned int petals;
  virtual bool used(Vertex n) = 0;
  virtual void setused(Vertex n) = 0;
};

void initialize_db(const Hypergraph &, const Graphstats &);
Leaf* lookup(const Edge &);

#endif
