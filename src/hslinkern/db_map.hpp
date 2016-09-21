#ifndef _DB_MAP_HPP
#define _DB_MAP_HPP

#include <set>
#include <Hypergraph.hpp>
#include "db.hpp"

using namespace std;

class MLeaf : public Leaf {
public:
  virtual bool used(Vertex n);
  virtual void setused(Vertex n);
private:
  set<Vertex> usedb;
};

void initialize_db(const Hypergraph &, const Graphstats &);
Leaf* lookup(const Edge &);

#endif
